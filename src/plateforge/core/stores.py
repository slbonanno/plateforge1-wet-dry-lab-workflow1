"""Named SQLite stores. No module ever hardcodes a database path.

Adding a store later = one register_store() call plus a .sql schema file.
"""
from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

from . import paths

SCHEMA_DIR = Path(__file__).parent / "schema"


@dataclass(frozen=True)
class StoreSpec:
    name: str
    schema: str
    description: str = ""

    @property
    def path(self) -> Path:
        return paths.store_path(self.name)


_STORES: dict[str, StoreSpec] = {}
_CONNS: dict[str, sqlite3.Connection] = {}


def register_store(name: str, schema: str, description: str = "") -> StoreSpec:
    if name in _STORES:
        raise KeyError(f"store {name!r} already registered")
    spec = StoreSpec(name, schema, description)
    _STORES[name] = spec
    return spec


def known_stores() -> dict[str, StoreSpec]:
    return dict(_STORES)


def connect(name: str) -> sqlite3.Connection:
    """Open (creating and migrating if needed) a named store."""
    if name in _CONNS:
        return _CONNS[name]
    if name not in _STORES:
        raise KeyError(f"unknown store {name!r}; known: {sorted(_STORES)}")
    spec = _STORES[name]
    conn = sqlite3.connect(spec.path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA journal_mode = WAL")
    conn.executescript((SCHEMA_DIR / spec.schema).read_text())
    conn.commit()
    _CONNS[name] = conn
    return conn


def attach(conn: sqlite3.Connection, name: str) -> None:
    """Make another store queryable in the same SQL statement, as `name.table`."""
    if name not in _STORES:
        raise KeyError(f"unknown store {name!r}")
    conn.execute("ATTACH DATABASE ? AS ?", (str(_STORES[name].path), name))


def close_all() -> None:
    for conn in _CONNS.values():
        conn.close()
    _CONNS.clear()


LEDGER = register_store("ledger", "ledger.sql", "artifact provenance across all modules")
LIBRARY = register_store("library", "library.sql", "in-house sequence database")
ASSAY = register_store("assay", "assay.sql", "plates, wells, reads")
REAGENTS = register_store("reagents", "reagents.sql", "reagent catalog, lots, caveats")
