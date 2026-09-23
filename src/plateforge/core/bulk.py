"""Parquet side of the library.

SQLite holds what you query and join; parquet holds what you scan in bulk.
The sequence pool will outgrow SQLite long before anything else does, so the
split exists from the start rather than being retrofitted.

Convention: a bulk table has a SQLite index table holding IDs and the few
columns you filter on, and a parquet file holding the full rows.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from . import paths


def write(name: str, df: pd.DataFrame, append: bool = False,
          key: str | None = None) -> Path:
    """Write (or append to) a bulk table. Returns the parquet path.

    `key` makes an append idempotent on its own terms. Without it, a parquet
    table can accumulate a second copy of everything whenever the SQLite index
    that normally guards against re-ingest is missing or was deleted -- the two
    halves of a store drift, and nothing notices until a count looks wrong.
    """
    path = paths.bulk_path(name)
    if append and path.exists():
        df = pd.concat([pd.read_parquet(path), df], ignore_index=True)
    if key and key in df.columns:
        df = df.drop_duplicates(subset=[key], keep="first").reset_index(drop=True)
    df.to_parquet(path, index=False)
    return path


def read(name: str, columns: list[str] | None = None, where: dict[str, Any] | None = None) -> pd.DataFrame:
    """Read a bulk table, optionally projecting columns and filtering on equality."""
    path = paths.bulk_path(name)
    if not path.exists():
        raise FileNotFoundError(f"no bulk table {name!r} at {path}")
    filters = [(k, "==", v) for k, v in (where or {}).items()] or None
    return pd.read_parquet(path, columns=columns, filters=filters)


def exists(name: str) -> bool:
    return paths.bulk_path(name).exists()


def info(name: str) -> dict[str, Any]:
    path = paths.bulk_path(name)
    if not path.exists():
        return {"name": name, "exists": False}
    df = pd.read_parquet(path, columns=None)
    return {"name": name, "exists": True, "rows": len(df),
            "columns": list(df.columns), "bytes": path.stat().st_size}
