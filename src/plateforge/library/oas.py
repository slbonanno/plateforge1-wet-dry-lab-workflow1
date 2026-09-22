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

# Columns we keep in the SQLite index. Everything else stays in parquet.
INDEX_COLUMNS = [
    "seq_id", "source", "source_ref", "chain", "v_call", "v_gene", "v_family",
    "j_call", "j_gene", "cdr3_aa", "cdr3_len", "aa_seq", "redundancy",
    "anarci_status", "has_liability", "species", "btype", "disease", "isotype",
]


def unit_url(study: str, filename: str, paired: bool = False) -> str:
    """Build a data unit URL. `filename` is as listed by the OAS web search."""
    kind = "paired" if paired else "unpaired"
    if not filename.endswith(".csv.gz"):
        filename = f"{filename}.csv.gz"
    return f"{BASE_URL}/{kind}/{study}/csv/{filename}"


def _open_text(src: str | Path):
    """Open a local .csv or .csv.gz as text. URLs are handed to pandas instead."""
    p = Path(src)
    if str(src).endswith(".gz"):
        return gzip.open(p, "rt", newline="")
    return open(p, "rt", newline="")


def read_metadata(src: str | Path) -> dict:
    """The JSON blob on line 0 of a data unit."""
    if str(src).startswith(("http://", "https://")):
        # pandas handles the URL and the gzip; the metadata arrives as headers.
        cols = pd.read_csv(src, nrows=0).columns
        return json.loads(",".join(cols))
    with _open_text(src) as fh:
        first = fh.readline()
    row = next(csv.reader(io.StringIO(first)))
    return json.loads(row[0]) if len(row) == 1 else json.loads(",".join(row))


def iter_unit(src: str | Path, chunksize: int = 50_000) -> Iterator[pd.DataFrame]:
    """Stream a data unit in chunks, skipping the metadata line."""
    yield from pd.read_csv(src, header=1, chunksize=chunksize, low_memory=False)


def read_unit(src: str | Path) -> pd.DataFrame:
    """Read a whole data unit at once. Prefer iter_unit for real ones."""
    return pd.read_csv(src, header=1, low_memory=False)


@dataclass
class Filter:
    """What to keep from a data unit. Every field is optional."""
    genes: list[str] | None = None          # e.g. ["IGHV3-23", "IGHV1-69"]
    families: list[str] | None = None       # e.g. ["IGHV3"]
    chain: str | None = None                # "H" or "L"
    cdr3_len_range: tuple[int, int] | None = None
    min_redundancy: int | None = None
    exclude_liabilities: bool = True
    productive_only: bool = True
    drop_duplicate_cdr3: bool = False

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df
        if self.productive_only and "productive" in out.columns:
            prod = out["productive"].astype(str).str.upper()
            out = out[prod.isin(["T", "TRUE", "YES"])]
        if self.chain:
            out = out[out["chain"] == self.chain]
        if self.genes:
            out = out[out["v_gene"].isin(self.genes)]
        if self.families:
            out = out[out["v_family"].isin(self.families)]
        if self.cdr3_len_range:
            lo, hi = self.cdr3_len_range
            out = out[(out["cdr3_len"] >= lo) & (out["cdr3_len"] <= hi)]
        if self.min_redundancy is not None and "redundancy" in out.columns:
            out = out[out["redundancy"].fillna(1) >= self.min_redundancy]
        if self.exclude_liabilities and "has_liability" in out.columns:
            out = out[~out["has_liability"].fillna(False)]
        if self.drop_duplicate_cdr3:
            out = out.drop_duplicates(subset=["cdr3_aa"])
        return out.reset_index(drop=True)


def _seq_id(source_ref: str, aa_seq: str) -> str:
    """Deterministic so re-ingesting the same unit does not duplicate rows.

    The body is a content hash, which makes ingest idempotent. It is still
    minted exactly once per distinct sequence and never recomputed downstream.
    """
    h = hashlib.blake2b(f"{source_ref}|{aa_seq}".encode(), digest_size=6).hexdigest()
    return ids.mint("SEQ", h)


def _has_liability(status) -> bool:
    """ANARCI_status is '|...|' with flags between the bars; empty means clean."""
    if not isinstance(status, str):
        return False
    return bool(status.strip("| ").strip())


def normalize(df: pd.DataFrame, meta: dict, source_ref: str) -> pd.DataFrame:
    """Turn raw OAS rows into the canonical pool schema.

    Unrecognised OAS columns are preserved; only the canonical ones are added.
    """
    out = df.copy()
    aa = out.get("sequence_alignment_aa")
    if aa is None:
        raise KeyError("data unit has no sequence_alignment_aa column")
    out["aa_seq"] = aa.astype(str).str.replace("-", "", regex=False).str.replace(".", "", regex=False)

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
    for chunk in iter_unit(src, chunksize=chunksize):
        total += len(chunk)
        part = filt.apply(normalize(chunk, meta, source_ref))
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
    return df, meta
