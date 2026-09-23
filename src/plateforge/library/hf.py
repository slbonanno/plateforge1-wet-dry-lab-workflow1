"""Reading OAS from the HuggingFace parquet mirror.

OPIG's own download paths return 403 (directory listings are readable, the
.csv.gz files are not), and OAS has no API. The mirror is the working route:

    ConvergeBio/oas-unpaired
      data/unpaired_heavy/<Study>/train-000NN.parquet
      data/unpaired_light/<Study>/train-000NN.parquet

It keeps the AIRR columns, including germline_alignment_aa and the per-region
columns the alignment figure needs. Two differences from a raw OAS data unit:

  * metadata arrives as meta_* COLUMNS, not a JSON header line
  * one folder holds every species, so a species filter is not optional --
    'unpaired_heavy' includes rabbit, mouse and rhesus alongside human

Parquet is read over range requests, so a column projection and a row-group
limit fetch megabytes rather than the whole 2.45 TB.
"""
from __future__ import annotations

import re
from typing import Iterator

import pandas as pd

from . import oas

REPO = "ConvergeBio/oas-unpaired"
HEAVY = "data/unpaired_heavy"
LIGHT = "data/unpaired_light"

# meta_<Key> in the mirror -> the key an OAS data unit puts in its JSON header.
META_PREFIX = "meta_"

# Columns worth pulling. Anything absent is skipped rather than erroring, so a
# schema change costs a missing feature instead of a crash.
WANTED = [
    "sequence_alignment_aa", "germline_alignment_aa", "v_germline_alignment_aa",
    "v_call", "j_call", "d_call", "productive",
    "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa",
    "junction_aa", "Redundancy", "ANARCI_status", "ANARCI_numbering",
    "v_identity", "j_identity",
]


def _require_hf():
    try:
        import huggingface_hub  # noqa: F401
    except ImportError as exc:                                     # pragma: no cover
        raise ImportError(
            "reading the mirror needs huggingface_hub and fsspec:\n"
            "    pip install huggingface_hub fsspec"
        ) from exc


def list_files(chain: str = "heavy") -> list[str]:
    """Every parquet path in the mirror for one chain."""
    _require_hf()
    from huggingface_hub import HfApi

    root = HEAVY if chain == "heavy" else LIGHT
    return sorted(f for f in HfApi().list_repo_files(REPO, repo_type="dataset")
                  if f.startswith(root) and f.endswith(".parquet"))


def list_studies(chain: str = "heavy") -> list[str]:
    """Study folder names, as the mirror spells them (e.g. 'Briney et al., 2019')."""
    root = HEAVY if chain == "heavy" else LIGHT
    seen = {f[len(root) + 1:].split("/")[0] for f in list_files(chain)}
    return sorted(seen)


def study_files(study: str, chain: str = "heavy") -> list[str]:
    """Shards belonging to one study. Matching is case- and punctuation-insensitive."""
    def norm(s: str) -> str:
        return re.sub(r"[^a-z0-9]", "", s.lower())

    root = HEAVY if chain == "heavy" else LIGHT
    target = norm(study)
    hits = [f for f in list_files(chain)
            if norm(f[len(root) + 1:].split("/")[0]) == target]
    if not hits:
        raise KeyError(
            f"no study matching {study!r}. Try list_studies() -- the mirror "
            "spells them like 'Briney et al., 2019'.")
    return hits


def _url(path: str) -> str:
    return f"hf://datasets/{REPO}/{path}"


def read_shard(path: str, columns: list[str] | None = None,
               row_groups: int | None = None) -> pd.DataFrame:
    """One parquet shard, optionally only the first N row groups."""
    _require_hf()
    import huggingface_hub  # noqa: F401  (registers the hf:// filesystem)
    import pyarrow.parquet as pq

    pf = pq.ParquetFile(_url(path))
    available = set(pf.schema_arrow.names)
    cols = [c for c in (columns or WANTED + meta_columns(pf)) if c in available]
    if row_groups is None:
        return pf.read(columns=cols).to_pandas()
    frames = [pf.read_row_group(i, columns=cols).to_pandas()
              for i in range(min(row_groups, pf.metadata.num_row_groups))]
    return pd.concat(frames, ignore_index=True)


def meta_columns(parquet_file) -> list[str]:
    return [c for c in parquet_file.schema_arrow.names if c.startswith(META_PREFIX)]


def extract_meta(df: pd.DataFrame) -> dict:
    """Turn the meta_* columns into the dict an OAS JSON header would carry."""
    out = {}
    for col in df.columns:
        if col.startswith(META_PREFIX) and len(df):
            vals = df[col].dropna().unique()
            if len(vals) == 1:
                out[col[len(META_PREFIX):]] = vals[0]
            elif len(vals) > 1:
                out[col[len(META_PREFIX):]] = f"<{len(vals)} values>"
    return out


def load_study(study: str, filt: oas.Filter | None = None, chain: str = "heavy",
               species: str | None = "human", limit: int | None = 20000,
               max_shards: int | None = None,
               row_groups_per_shard: int | None = 1) -> tuple[pd.DataFrame, dict]:
    """Load one study from the mirror, normalized into the pool schema.

    `species` is applied before anything else and defaults to human, because a
    chain folder in the mirror mixes species and a rabbit IGHV1S45 sequence
    would otherwise sail into a human panel unnoticed. Pass None to keep all.
    """
    filt = filt or oas.Filter()
    shards = study_files(study, chain)
    if max_shards:
        shards = shards[:max_shards]

    kept: list[pd.DataFrame] = []
    scanned = 0
    used_meta: dict = {}
    funnel = None
    species_dropped = 0
    for path in shards:
        raw = read_shard(path, row_groups=row_groups_per_shard)
        scanned += len(raw)
        meta = extract_meta(raw)
        used_meta = used_meta or meta

        if species is not None and "meta_Species" in raw.columns:
            before = len(raw)
            raw = raw[raw["meta_Species"].astype(str).str.lower() == species.lower()]
            species_dropped += before - len(raw)
        if not len(raw):
            continue

        normalized = oas.normalize(raw, meta,
                                   source_ref=f"hf:{study}:{path.split('/')[-1]}")
        if funnel is None:
            funnel = filt.explain(normalized)
        part = filt.apply(normalized)
        if len(part):
            kept.append(part)
        if limit and sum(len(k) for k in kept) >= limit:
            break

    df = (pd.concat(kept, ignore_index=True) if kept
          else pd.DataFrame(columns=oas.INDEX_COLUMNS))
    if limit:
        df = df.head(limit)
    df = df.drop_duplicates(subset=["seq_id"]).reset_index(drop=True)

    used_meta = dict(used_meta)
    used_meta.update({
        "_rows_scanned": scanned,
        "_rows_kept": len(df),
        "_source_ref": f"hf:{REPO}:{study}",
        "_shards_read": len(shards),
        "_species_dropped": species_dropped,
        "_funnel": funnel,
    })
    return df, used_meta


def alignment_report(df: pd.DataFrame) -> pd.DataFrame:
    """Per germline: are the sequences the same length, i.e. column-aligned?

    The alignment figure assumes IMGT-gapped sequences of one germline share
    columns. That holds for OAS's own files; whether it survives the mirror is
    a question about real data, so it gets measured rather than assumed.
    """
    rows = []
    for gene, sub in df.groupby("v_gene"):
        lens = sub["aa_gapped"].str.len()
        rows.append({
            "v_gene": gene,
            "n": len(sub),
            "distinct_lengths": int(lens.nunique()),
            "min_len": int(lens.min()),
            "max_len": int(lens.max()),
            "column_aligned": bool(lens.nunique() == 1),
        })
    return pd.DataFrame(rows).sort_values("n", ascending=False).reset_index(drop=True)
