"""Where files live. Nothing constructs a data path by hand.

    $PLATEFORGE_DATA/
      stores/            named SQLite databases
      bulk/              large tables as parquet (too big for SQLite)
      raw/<source>/      downloaded external files, never modified after landing
      outputs/<obj_id>/  everything a module produced for one artifact

Point PLATEFORGE_DATA at an external drive and the whole tree moves with it.
"""
from __future__ import annotations

import os
from pathlib import Path


def data_root() -> Path:
    root = Path(os.environ.get("PLATEFORGE_DATA", Path.cwd() / "data")).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    return root


def _sub(name: str) -> Path:
    p = data_root() / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def store_path(name: str) -> Path:
    return _sub("stores") / f"{name}.sqlite"


def bulk_path(name: str) -> Path:
    """Parquet file for a table too large to keep in SQLite."""
    return _sub("bulk") / f"{name}.parquet"


def raw_dir(source: str) -> Path:
    """Landing zone for downloads. Treat contents as read-only once written."""
    d = _sub("raw") / source
    d.mkdir(parents=True, exist_ok=True)
    return d


def output_dir(obj_id: str) -> Path:
    """One directory per artifact. The artifact's `path` field points here."""
    d = _sub("outputs") / obj_id
    d.mkdir(parents=True, exist_ok=True)
    return d
