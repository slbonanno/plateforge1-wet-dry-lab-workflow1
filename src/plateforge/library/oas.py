"""Reading Observed Antibody Space data units.

An OAS data unit is a gzipped CSV with a one-line JSON metadata header:

    line 0   {"Species": "human", "Chain": "Heavy", ...}   (CSV-quoted JSON)
    line 1   AIRR column names
    line 2+  one sequence per row

Full OAS is ~1.1 TB, so nothing here downloads the database. A unit is read
straight from its URL or from a file already on disk, filtered while streaming,
and only the surviving rows are kept.

Format notes and a fixture live in docs/formats/oas.md.
"""
from __future__ import annotations

import csv
import gzip
import os
import hashlib
import io
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import pandas as pd

from ..core import ids
from . import germlines

BASE_URL = "https://opig.stats.ox.ac.uk/webapps/ngsdb"

# OPIG's server rejects the default Python-urllib User-Agent with 403, so a
# browser-like one is sent for URL reads. Override with the PLATEFORGE_UA
# environment variable if it ever needs to change.
USER_AGENT = os.environ.get(
    "PLATEFORGE_UA",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/125.0 Safari/537.36")


def _is_url(src) -> bool:
    return str(src).startswith(("http://", "https://"))


def _storage_options(src) -> dict | None:
    """pandas passes these as HTTP headers for URL reads, and rejects them
    for local paths, so they are only supplied for URLs."""
    return {"User-Agent": USER_AGENT} if _is_url(src) else None


# IMGT regions, in order along the sequence. OAS ships each as its own column.
REGION_COLUMNS = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa"]

# Columns we keep in the SQLite index. Everything else stays in parquet.
INDEX_COLUMNS = [
    "seq_id", "source", "source_ref", "chain", "v_call", "v_gene", "v_family",
    "j_call", "j_gene", "cdr3_aa", "cdr3_len", "aa_seq", "redundancy", "n_ambiguous",
    "anarci_status", "has_liability", "species", "btype", "disease", "isotype",
]


def unit_url(study: str, filename: str, paired: bool = False,
             scheme: str = "https") -> str:
    """Build a data unit URL. `filename` is as listed by the OAS web search."""
    kind = "paired" if paired else "unpaired"
    if not filename.endswith(".csv.gz"):
        filename = f"{filename}.csv.gz"
    base = BASE_URL if scheme == "https" else BASE_URL.replace("https://", "http://", 1)
    return f"{base}/{kind}/{study}/csv/{filename}"


def _open_text(src: str | Path):
    """Open a local .csv or .csv.gz as text. URLs are handed to pandas instead."""
    p = Path(src)
    if str(src).endswith(".gz"):
        return gzip.open(p, "rt", newline="")
    return open(p, "rt", newline="")


def read_metadata(src: str | Path) -> dict:
    """The JSON blob on line 0 of a data unit."""
    if _is_url(src):
        # pandas handles the URL and the gzip; the metadata arrives as headers.
        cols = pd.read_csv(src, nrows=0, storage_options=_storage_options(src)).columns
        return json.loads(",".join(cols))
    with _open_text(src) as fh:
        first = fh.readline()
    row = next(csv.reader(io.StringIO(first)))
    return json.loads(row[0]) if len(row) == 1 else json.loads(",".join(row))


def iter_unit(src: str | Path, chunksize: int = 50_000) -> Iterator[pd.DataFrame]:
    """Stream a data unit in chunks, skipping the metadata line."""
    yield from pd.read_csv(src, header=1, chunksize=chunksize, low_memory=False,
                           storage_options=_storage_options(src))


def read_unit(src: str | Path) -> pd.DataFrame:
    """Read a whole data unit at once. Prefer iter_unit for real ones."""
    return pd.read_csv(src, header=1, low_memory=False,
                       storage_options=_storage_options(src))


@dataclass
class Filter:
    """What to keep from a data unit. Every field is optional.

    Filtering is a funnel, and a funnel that silently empties is useless. Each
    step is named, so `explain()` can show exactly where rows died instead of
    leaving "kept 0" as the only clue.
    """
    genes: list[str] | None = None          # e.g. ["IGHV3-23", "IGHV1-69"]
    families: list[str] | None = None       # e.g. ["IGHV3"]
    chain: str | None = None                # "H" or "L"
    cdr3_len_range: tuple[int, int] | None = None
    min_redundancy: int | None = None
    exclude_liabilities: bool = True
    productive_only: bool = True
    drop_duplicate_cdr3: bool = False
    max_ambiguous: int | None = 0   # X and other non-standard residues allowed

    def _steps(self):
        """(name, function) pairs, applied in order."""
        steps = []
        if self.productive_only:
            def productive(df):
                # Studies differ: 'T'/'F', True/False, 1/0, or not populated.
                # A null is unknown, not unproductive -- dropping those would
                # silently empty a whole study, so they are kept.
                if "productive" not in df.columns:
                    return df
                col = df["productive"]
                if col.dtype == bool:
                    return df[col]
                text = col.astype(str).str.strip().str.upper()
                unknown = col.isna() | text.isin(["NAN", "NONE", ""])
                return df[unknown | text.isin(["T", "TRUE", "YES", "1"])]
            steps.append(("productive", productive))
        if self.chain:
            steps.append((f"chain == {self.chain}",
                          lambda df: df[df["chain"] == self.chain]))
        if self.genes:
            steps.append((f"v_gene in {len(self.genes)} panel genes",
                          lambda df: df[df["v_gene"].isin(self.genes)]))
        if self.families:
            steps.append((f"v_family in {self.families}",
                          lambda df: df[df["v_family"].isin(self.families)]))
        if self.cdr3_len_range:
            lo, hi = self.cdr3_len_range
            steps.append((f"cdr3_len {lo}-{hi}",
                          lambda df: df[(df["cdr3_len"] >= lo) & (df["cdr3_len"] <= hi)]))
        if self.min_redundancy is not None:
            steps.append((f"redundancy >= {self.min_redundancy}",
                          lambda df: df[df["redundancy"].fillna(1) >= self.min_redundancy]
                          if "redundancy" in df.columns else df))
        if self.exclude_liabilities:
            steps.append(("no ANARCI liability",
                          lambda df: df[~df["has_liability"].fillna(False)]
                          if "has_liability" in df.columns else df))
        if self.max_ambiguous is not None:
            cap = self.max_ambiguous
            steps.append((f"ambiguous residues <= {cap}",
                          lambda df: df[df["n_ambiguous"].fillna(0) <= cap]
                          if "n_ambiguous" in df.columns else df))
        if self.drop_duplicate_cdr3:
            steps.append(("unique cdr3", lambda df: df.drop_duplicates(subset=["cdr3_aa"])))
        return steps

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df
        for _name, fn in self._steps():
            out = fn(out)
        return out.reset_index(drop=True)

    def explain(self, df: pd.DataFrame) -> pd.DataFrame:
        """Rows surviving each step, in order. Run this when a filter keeps
        nothing: the first step that drops everything is the culprit."""
        rows = [{"step": "input", "rows": len(df), "dropped": 0}]
        out = df
        for name, fn in self._steps():
            before = len(out)
            out = fn(out)
            rows.append({"step": name, "rows": len(out), "dropped": before - len(out)})
        return pd.DataFrame(rows)


def _seq_id(source_ref: str, aa_seq: str) -> str:
    """Deterministic so re-ingesting the same unit does not duplicate rows.

    The body is a content hash, which makes ingest idempotent. It is still
    minted exactly once per distinct sequence and never recomputed downstream.
    """
    h = hashlib.blake2b(f"{source_ref}|{aa_seq}".encode(), digest_size=6).hexdigest()
    return ids.mint("SEQ", h)


# ANARCI_status is pipe-delimited flags, and most of them describe the READ
# rather than the molecule. Measured on Briney et al. 2019 (19,864 sequences
# after gene and length filters):
#
#   Shorter than IMGT defined: fw1                       17,592
#   Deletions: 1..15, 73                                 11,471
#   Deletions: 1..14, 73                                  3,453
#   Shorter than IMGT defined: fw1, fw4                   2,272
#   Missing Conserved Cysteine: 104                          480
#
# Treating every flag as a liability discarded all 125,743 rows, which is how
# this was found. "Shorter than IMGT defined" and low-numbered "Deletions" are
# both the same artifact: the amplicon starts inside FR1, so the first dozen
# IMGT positions are simply not covered. Calling that a liability would reject
# an entire repertoire on the basis of primer design.
#
# A judgment call worth recording: "Deletions" is treated as coverage even
# though a mid-domain deletion would be real. In this dataset position 73 is
# deleted in ~90% of sequences, which is systematic, not 18,000 real events.
# If a dataset ever shows sparse isolated deletions, revisit this.
#
# Missing Conserved Cysteine is a genuine structural problem and stays a
# liability. So does anything unrecognised: unknown flags are liabilities until
# someone looks at them.
COVERAGE_PATTERNS = (
    "shorter than imgt defined",
    "deletions:",
)


def status_flags(status) -> list[str]:
    """The individual flags in an ANARCI_status string."""
    if not isinstance(status, str):
        return []
    return [f.strip() for f in status.split("|") if f.strip()]


def is_coverage_flag(flag: str, patterns: tuple[str, ...] = COVERAGE_PATTERNS) -> bool:
    """True when a flag describes incomplete read coverage, not the molecule."""
    low = flag.lower()
    return any(p in low for p in patterns)


def _has_liability(status, benign: tuple[str, ...] = COVERAGE_PATTERNS) -> bool:
    """True when a sequence carries a flag that is about the molecule."""
    return any(not is_coverage_flag(f, benign) for f in status_flags(status))


def flag_counts(df: pd.DataFrame) -> pd.Series:
    """How often each ANARCI flag appears. Run this on real data before
    trusting any liability filter."""
    from collections import Counter
    counts: Counter = Counter()
    for status in df.get("anarci_status", pd.Series(dtype=object)):
        counts.update(status_flags(status))
    return pd.Series(counts).sort_values(ascending=False)


def normalize(df: pd.DataFrame, meta: dict, source_ref: str) -> pd.DataFrame:
    """Turn raw OAS rows into the canonical pool schema.

    Unrecognised OAS columns are preserved; only the canonical ones are added.
    """
    out = df.copy()
    aa = out.get("sequence_alignment_aa")
    if aa is None:
        raise KeyError("data unit has no sequence_alignment_aa column")
    # Keep BOTH forms. The gapped one is IMGT-numbered, so sequences from one
    # germline are already column-aligned -- that is what makes an alignment
    # view possible without running an aligner.
    out["aa_gapped"] = aa.astype(str)
    out["aa_seq"] = out["aa_gapped"].str.replace("-", "", regex=False).str.replace(".", "", regex=False)

    # OAS ships IgBlast's own germline reference per sequence, already aligned
    # column-for-column with sequence_alignment_aa. That is the authoritative
    # reference for a divergence view -- no germline table to maintain, and no
    # need to fall back on the consensus of whatever happens to be in the pool.
    germ = out.get("germline_alignment_aa")
    if germ is None:
        germ = out.get("v_germline_alignment_aa")
    out["germline_aa"] = germ.astype(str) if germ is not None else None

    # Per-region amino acid strings, used to draw FR/CDR boundaries on an
    # alignment. IMGT regions, as IgBlast reports them.
    for col in REGION_COLUMNS:
        if col not in out.columns:
            out[col] = None

    out["v_gene"] = out["v_call"].map(germlines.gene)
    out["v_family"] = out["v_call"].map(germlines.family)
    out["j_gene"] = out.get("j_call", pd.Series(dtype=object)).map(germlines.gene)
    out["chain"] = out["v_call"].map(germlines.chain_of)

    cdr3 = out.get("cdr3_aa")
    if cdr3 is None:
        raise KeyError("data unit has no cdr3_aa column")
    out["cdr3_aa"] = cdr3.astype(str)
    out["cdr3_len"] = out["cdr3_aa"].str.len()

    status = out.get("ANARCI_status")
    out["anarci_status"] = status if status is not None else ""
    out["has_liability"] = out["anarci_status"].map(_has_liability)

    # X marks a position the basecaller could not resolve. Real repertoires
    # carry plenty; they are not biology and must never reach a picked clone.
    out["n_ambiguous"] = out["aa_seq"].str.count(r"[^ACDEFGHIKLMNPQRSTVWY]")

    out["redundancy"] = pd.to_numeric(out.get("Redundancy", 1), errors="coerce").fillna(1).astype(int)
    out["isotype"] = out.get("Isotype", meta.get("Isotype", ""))

    out["source"] = "OAS"
    out["source_ref"] = source_ref
    out["species"] = meta.get("Species", "")
    out["btype"] = meta.get("BType", "")
    out["disease"] = meta.get("Disease", "")
    out["study"] = meta.get("Author", "") or meta.get("Study", "")

    out["seq_id"] = [_seq_id(source_ref, s) for s in out["aa_seq"]]
    return out


def load(src: str | Path, filt: Filter | None = None, chunksize: int = 50_000,
         limit: int | None = None) -> tuple[pd.DataFrame, dict]:
    """Read, normalize and filter a data unit. Returns (rows, metadata).

    Works the same for a URL and a local file, so the download step can happen
    once and everything after it runs offline.
    """
    meta = read_metadata(src)
    source_ref = Path(str(src)).name
    filt = filt or Filter()
    kept: list[pd.DataFrame] = []
    total = 0
    funnel = None
    for chunk in iter_unit(src, chunksize=chunksize):
        total += len(chunk)
        normalized = normalize(chunk, meta, source_ref)
        if funnel is None:
            funnel = filt.explain(normalized)
        part = filt.apply(normalized)
        if len(part):
            kept.append(part)
        if limit and sum(len(k) for k in kept) >= limit:
            break
    df = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=INDEX_COLUMNS)
    if limit:
        df = df.head(limit)
    df = df.drop_duplicates(subset=["seq_id"]).reset_index(drop=True)
    meta = dict(meta)
    meta["_rows_scanned"] = total
    meta["_rows_kept"] = len(df)
    meta["_source_ref"] = source_ref
    meta["_funnel"] = funnel
    return df, meta
