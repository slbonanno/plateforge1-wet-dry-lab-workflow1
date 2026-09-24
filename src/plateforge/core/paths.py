"""Where files live. Nothing constructs a data path by hand.

    <repo>/plateforge-data/        (the default; gitignored)
      stores/            named SQLite databases
      bulk/              large tables as parquet (too big for SQLite)
      raw/<source>/      downloaded external files, never modified after landing
      outputs/<obj_id>/  everything a module produced for one artifact

Point PLATEFORGE_DATA at an external drive and the whole tree moves with it.

The default sits inside the checkout rather than in the home directory, so a
project is one directory: clone it, run it, delete it and nothing is left
behind. `plateforge-data/` is gitignored -- a pool of a few hundred thousand
sequences has no business in git.

The repo root is found from this file, not from the working directory. A cwd
default means the same command writes to a different store depending on which
directory you happened to be in, and two half-populated pools is a bad
afternoon.
"""
from __future__ import annotations

import os
from pathlib import Path


DEFAULT_DIRNAME = "plateforge-data"


def repo_root() -> Path | None:
    """The checkout this package was installed from, if it is one."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pyproject.toml").exists():
            return parent
    return None


def default_root() -> Path:
    """Inside the checkout when there is one, else beside the caller."""
    repo = repo_root()
    return (repo or Path.cwd()) / DEFAULT_DIRNAME


def data_root() -> Path:
    root = Path(os.environ.get("PLATEFORGE_DATA") or default_root()).expanduser()
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
