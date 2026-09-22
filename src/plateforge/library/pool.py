"""The sequence pool: a growing store of sequences to sample from later.

Split per decisions/0006: full rows in parquet ("oas_pool"), filter columns in
the SQLite `sequences` table. Ingest is idempotent because seq_id is a content
hash, so re-running the same unit adds nothing.
"""
from __future__ import annotations

import json
from typing import Any

import pandas as pd

from ..core import artifacts, bulk, ids, stores
from . import oas

POOL_TABLE = "oas_pool"


def ingest(df: pd.DataFrame, meta: dict, produced_by: str = "library.oas") -> str:
    """Add normalized rows to the pool. Returns the ingest RUN artifact id."""
    run_id = ids.mint_stamped("RUN", "ingest")
    if len(df):
        existing = set(known_ids())
        fresh = df[~df["seq_id"].isin(existing)].reset_index(drop=True)
    else:
        fresh = df

    if len(fresh):
        bulk.write(POOL_TABLE, fresh, append=True)
        conn = stores.connect("library")
        cols = [c for c in oas.INDEX_COLUMNS if c in fresh.columns]
        rows = fresh[cols].copy()
        rows["has_liability"] = rows.get("has_liability", False).astype(int)
        rows["meta_json"] = "{}"
        rows["notes"] = None
        rows["created_at"] = ids.now_iso()
        rows.to_sql("sequences", conn, if_exists="append", index=False)
        conn.commit()

    artifacts.register(
        run_id, "sequence_ingest", produced_by,
        label=meta.get("_source_ref"),
        store="library", path=str(bulk.paths.bulk_path(POOL_TABLE)),
        params={"source_ref": meta.get("_source_ref")},
        meta={k: v for k, v in meta.items() if not k.startswith("_")}
             | {"rows_scanned": meta.get("_rows_scanned"),
                "rows_kept": meta.get("_rows_kept"),
                "rows_added": int(len(fresh))},
    )
    return run_id


def known_ids() -> set[str]:
    conn = stores.connect("library")
    return {r[0] for r in conn.execute("SELECT seq_id FROM sequences")}


def index(where: str | None = None, params: tuple = ()) -> pd.DataFrame:
    """Query the SQLite index. `where` is raw SQL without the WHERE keyword."""
    conn = stores.connect("library")
    sql = "SELECT * FROM sequences" + (f" WHERE {where}" if where else "")
    return pd.read_sql(sql, conn, params=params)


def fetch(seq_ids: list[str], columns: list[str] | None = None) -> pd.DataFrame:
    """Full parquet rows for a set of seq_ids."""
    df = bulk.read(POOL_TABLE, columns=columns)
    return df[df["seq_id"].isin(seq_ids)].reset_index(drop=True)


def stats() -> dict[str, Any]:
    df = index()
    if df.empty:
        return {"total": 0, "by_gene": {}, "by_source": {}}
    return {
        "total": int(len(df)),
        "by_gene": df["v_gene"].value_counts().to_dict(),
        "by_source": df["source_ref"].value_counts().to_dict(),
        "cdr3_len": {
            "min": int(df["cdr3_len"].min()),
            "median": float(df["cdr3_len"].median()),
            "max": int(df["cdr3_len"].max()),
        },
    }


def make_library(name: str, seq_ids: list[str], *, produced_by: str,
                 params: dict | None = None, parents: dict[str, str] | None = None,
                 meta: dict | None = None) -> str:
    """Freeze a set of sequences as a named LIB artifact."""
    lib_id = ids.mint_stamped("LIB", name)
    conn = stores.connect("library")
    conn.executemany(
        "INSERT OR REPLACE INTO library_members (lib_id, seq_id, rank) VALUES (?,?,?)",
        [(lib_id, s, i) for i, s in enumerate(seq_ids)],
    )
    conn.commit()
    artifacts.register(
        lib_id, "sequence_set", produced_by, label=name, store="library",
        params=params or {}, parents=parents or {},
        meta=(meta or {}) | {"n": len(seq_ids)},
    )
    return lib_id


def library_members(lib_id: str) -> pd.DataFrame:
    conn = stores.connect("library")
    return pd.read_sql(
        "SELECT s.*, m.rank FROM library_members m JOIN sequences s USING (seq_id)"
        " WHERE m.lib_id = ? ORDER BY m.rank", conn, params=(lib_id,))
