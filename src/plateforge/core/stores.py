"""Named SQLite stores. No module ever hardcodes a database path.

Adding a store later = one register_store() call plus a .sql schema file.
"""
from __future__ import annotations

import re
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


_CREATE_TABLE = re.compile(
    r"CREATE\s+TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?([A-Za-z_][\w]*)\s*\((.*?)\)\s*;",
    re.IGNORECASE | re.DOTALL)

# Constraints SQLite cannot attach with ALTER TABLE ADD COLUMN.
_UNADDABLE = re.compile(r"\b(PRIMARY\s+KEY|UNIQUE)\b", re.IGNORECASE)
_TABLE_CONSTRAINT = re.compile(
    r"^\s*(PRIMARY\s+KEY|FOREIGN\s+KEY|UNIQUE|CHECK|CONSTRAINT)\b", re.IGNORECASE)


def _split_columns(body: str) -> list[str]:
    """Split a CREATE TABLE body on top-level commas."""
    parts, depth, current = [], 0, []
    for ch in body:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    if current:
        parts.append("".join(current))
    return [p.strip() for p in parts if p.strip()]


def strip_comments(sql: str) -> str:
    """Remove -- comments, preserving line structure.

    This has to happen BEFORE splitting a table body on commas, not after.
    A comment containing a comma -- and a comment enumerating the allowed
    values of a column almost always does -- otherwise splits one column
    definition into two pieces. The piece before the comma ends inside the
    comment and is discarded as empty, so the NEXT column silently vanishes
    from the migration; the piece after begins mid-sentence and is read as a
    column definition, which is how a table ends up being asked for a column
    called `not`.

    Silently skipping a column is the bad half: the schema declares it, the
    database never gains it, and the failure surfaces later as a missing
    column at write time.
    """
    return "\n".join(line.split("--", 1)[0] for line in sql.splitlines())


def declared_columns(schema_sql: str) -> dict[str, dict[str, str]]:
    """{table: {column: full definition}} as the schema file declares them."""
    out: dict[str, dict[str, str]] = {}
    for table, body in _CREATE_TABLE.findall(strip_comments(schema_sql)):
        cols: dict[str, str] = {}
        for piece in _split_columns(body):
            line = " ".join(piece.split()).strip()
            if not line or _TABLE_CONSTRAINT.match(line):
                continue
            cols[line.split()[0]] = line
        out[table] = cols
    return out


def migrate(conn: sqlite3.Connection, schema_sql: str) -> list[str]:
    """Add columns the schema declares but the database lacks.

    Schemas are expected to gain columns (CLAUDE.md rule 7), and
    CREATE TABLE IF NOT EXISTS does not alter an existing table -- so without
    this, adding one column breaks every store that already holds data.

    Additive only, deliberately. Renames, type changes and drops are not
    guessed at: they need a considered migration, and silently rewriting a
    table is how data disappears.
    """
    applied: list[str] = []
    for table, cols in declared_columns(schema_sql).items():
        existing_rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
        if not existing_rows:
            continue                       # table is new; the script created it
        existing = {r[1] for r in existing_rows}
        for name, definition in cols.items():
            if name in existing:
                continue
            addable = _UNADDABLE.sub("", definition).strip()
            if re.search(r"\bNOT\s+NULL\b", addable, re.IGNORECASE) and \
               not re.search(r"\bDEFAULT\b", addable, re.IGNORECASE):
                # SQLite rejects NOT NULL without a default on a populated table.
                addable = re.sub(r"\bNOT\s+NULL\b", "", addable, flags=re.IGNORECASE).strip()
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {addable}")
            applied.append(f"{table}.{name}")
    if applied:
        conn.commit()
    return applied


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
    schema_sql = (SCHEMA_DIR / spec.schema).read_text()
    # Tables first, then migrate, then the whole script. An index declared on a
    # column an older database lacks would fail if the script ran in one go.
    tables_only = "\n".join(
        f"CREATE TABLE IF NOT EXISTS {t} (\n{body}\n);"
        for t, body in _CREATE_TABLE.findall(schema_sql))
    conn.executescript(tables_only)
    migrate(conn, schema_sql)
    conn.executescript(schema_sql)
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
