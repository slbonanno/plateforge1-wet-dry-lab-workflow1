"""True column alignment from ANARCI IMGT numbering.

`sequence_alignment_aa` is NOT a fixed-width alignment. Measured on Briney et
al. 2019, IGHV3-23 has 54 distinct lengths across 16,497 sequences, because
reads cover different spans of the domain. Laying those out side by side by
string position puts every residue in the wrong column.

OAS carries `ANARCI_numbering`: a per-residue IMGT position label, grouped by
region. That is an alignment — each residue already knows which column it
belongs in, so sequences of any length drop into shared columns with no
aligner and no guessing.

The stored value is a Python dict literal:

    {'fwh1': {'15 ': 'P', '16 ': 'G', ...},
     'cdrh1': {'27 ': 'G', ...},
     ...
     'cdrh3': {'111 ': 'G', '111A': 'I', '112A': 'D', '112 ': 'R', ...},
     'fwh4': {'118 ': 'W', ...}}

Position keys are IMGT numbers, space-padded, with a letter for insertions.
"""
from __future__ import annotations

import ast
import re
from dataclasses import dataclass

import pandas as pd

# Region keys as ANARCI emits them, in order along the chain. Heavy and light
# differ only in the letter, so both are accepted.
REGION_ORDER_H = ["fwh1", "cdrh1", "fwh2", "cdrh2", "fwh3", "cdrh3", "fwh4"]
REGION_ORDER_L = ["fwl1", "cdrl1", "fwl2", "cdrl2", "fwl3", "cdrl3", "fwl4"]
REGION_LABELS = {
    "fwh1": "FR1", "cdrh1": "CDR1", "fwh2": "FR2", "cdrh2": "CDR2",
    "fwh3": "FR3", "cdrh3": "CDR3", "fwh4": "FR4",
    "fwl1": "FR1", "cdrl1": "CDR1", "fwl2": "FR2", "cdrl2": "CDR2",
    "fwl3": "FR3", "cdrl3": "CDR3", "fwl4": "FR4",
}

_POS = re.compile(r"^\s*(\d+)\s*([A-Z]*)\s*$")


def parse(value) -> dict[str, dict[str, str]] | None:
    """Parse a stored ANARCI_numbering value. Returns None if unusable."""
    if isinstance(value, dict):
        return value
    if not isinstance(value, str) or not value.strip().startswith("{"):
        return None
    try:
        parsed = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return None
    return parsed if isinstance(parsed, dict) else None


def position_key(pos: str) -> tuple[int, int, str]:
    """Sort key placing an IMGT position, insertions included, in column order.

    IMGT inserts into CDR3 outward from the middle: 111, 111A, 111B, then
    112C, 112B, 112A, 112. So insertions after 111 ascend and insertions
    before 112 descend, which is why this is not a plain string sort.
    """
    m = _POS.match(str(pos))
    if not m:
        return (10_000, 0, str(pos))
    num, ins = int(m.group(1)), m.group(2)
    if not ins:
        return (num, 0, "")
    rank = sum((ord(c) - 64) * (27 ** (len(ins) - i - 1)) for i, c in enumerate(ins))
    return (num, -rank, ins) if num == 112 else (num, rank, ins)


@dataclass
class Alignment:
    """Residues laid out in shared IMGT columns."""
    matrix: pd.DataFrame            # rows = seq_id, columns = IMGT position labels
    regions: list[tuple[str, int, int]]   # (label, start_col, end_col) half-open

    @property
    def width(self) -> int:
        return self.matrix.shape[1]

    def as_strings(self, gap: str = ".") -> dict[str, str]:
        return {idx: "".join(row.fillna(gap))
                for idx, row in self.matrix.iterrows()}

    def occupancy(self) -> pd.Series:
        """Fraction of sequences with a residue at each column. Columns near
        zero are insertions only a few sequences carry."""
        return self.matrix.notna().mean()


def build(rows: pd.DataFrame, numbering_column: str = "anarci_numbering",
          id_column: str = "seq_id", min_occupancy: float = 0.0) -> Alignment | None:
    """Lay sequences into shared IMGT columns. None if no numbering is usable.

    `min_occupancy` drops columns held by fewer than that fraction of
    sequences, which trims rare CDR3 insertions that would otherwise stretch
    the figure with near-empty columns.
    """
    if numbering_column not in rows.columns:
        return None

    per_seq: dict[str, dict[str, str]] = {}
    col_region: dict[str, str] = {}
    for _, row in rows.iterrows():
        parsed = parse(row[numbering_column])
        if not parsed:
            continue
        flat: dict[str, str] = {}
        for region, residues in parsed.items():
            if not isinstance(residues, dict):
                continue
            label = REGION_LABELS.get(region, region)
            for pos, aa in residues.items():
                key = str(pos).strip()
                flat[key] = aa
                col_region.setdefault(key, label)
        if flat:
            per_seq[str(row[id_column])] = flat

    if not per_seq:
        return None

    columns = sorted({c for flat in per_seq.values() for c in flat}, key=position_key)
    matrix = pd.DataFrame.from_dict(per_seq, orient="index").reindex(columns=columns)

    if min_occupancy > 0:
        keep = matrix.notna().mean() >= min_occupancy
        matrix = matrix.loc[:, keep]
        columns = list(matrix.columns)

    regions: list[tuple[str, int, int]] = []
    for i, col in enumerate(columns):
        label = col_region.get(col, "?")
        if regions and regions[-1][0] == label:
            regions[-1] = (label, regions[-1][1], i + 1)
        else:
            regions.append((label, i, i + 1))
    return Alignment(matrix=matrix, regions=regions)


def germline_by_column(row: pd.Series, numbering_column: str = "anarci_numbering",
                       seq_column: str = "aa_gapped",
                       germline_column: str = "germline_aa") -> dict[str, str] | None:
    """Place one sequence's germline residues into IMGT columns.

    `germline_alignment_aa` is aligned to `sequence_alignment_aa` position by
    position, and the numbering covers that sequence's non-gap residues in
    order. Zipping the two puts the germline into the same columns. Returns
    None if the counts disagree, rather than producing a shifted reference --
    a germline row one residue out would make every position look mutated.
    """
    parsed = parse(row.get(numbering_column))
    seq, germ = row.get(seq_column), row.get(germline_column)
    if not parsed or not isinstance(seq, str) or not isinstance(germ, str):
        return None
    if len(seq) != len(germ):
        return None

    keys = sorted({str(p).strip() for residues in parsed.values()
                   if isinstance(residues, dict) for p in residues},
                  key=position_key)
    pairs = [(s, g) for s, g in zip(seq, germ) if s not in ".-"]
    if len(pairs) != len(keys):
        return None
    return {key: g for key, (_s, g) in zip(keys, pairs)}


def germline_row(rows: pd.DataFrame, columns: list[str], **kwargs) -> tuple[dict[str, str], int]:
    """A germline reference per column, taken as the modal residue across rows.

    Returns (mapping, how many rows contributed). Rows whose numbering and
    germline disagree are skipped rather than guessed at.
    """
    from collections import Counter, defaultdict

    votes: dict[str, Counter] = defaultdict(Counter)
    used = 0
    for _, row in rows.iterrows():
        mapped = germline_by_column(row, **kwargs)
        if not mapped:
            continue
        used += 1
        for col, aa in mapped.items():
            votes[col][aa] += 1
    return ({col: votes[col].most_common(1)[0][0]
             for col in columns if votes.get(col)}, used)
