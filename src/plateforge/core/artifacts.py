"""The contract between modules.

Modules do not import each other's functions. A module writes something,
registers it here, and returns an obj_id. The next module resolves that id.
Everything is traceable both directions through artifact_links.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from . import ids, stores, version


@dataclass
class Artifact:
    obj_id: str
    obj_type: str
    produced_by: str
    label: str | None = None
    store: str | None = None
    path: str | None = None
    code_version: str | None = None
    params: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)
    notes: str | None = None
    created_at: str = ""

    @classmethod
    def _from_row(cls, row) -> "Artifact":
        return cls(
            obj_id=row["obj_id"],
            obj_type=row["obj_type"],
            produced_by=row["produced_by"],
            label=row["label"],
            store=row["store"],
            path=row["path"],
            code_version=row["code_version"],
            params=json.loads(row["params_json"]),
            meta=json.loads(row["meta_json"]),
            notes=row["notes"],
            created_at=row["created_at"],
        )


def register(
    obj_id: str,
    obj_type: str,
    produced_by: str,
    *,
    label: str | None = None,
    store: str | None = None,
    path: str | None = None,
    code_version: str | None = None,
    params: dict | None = None,
    meta: dict | None = None,
    notes: str | None = None,
    parents: dict[str, str] | None = None,
) -> str:
    """Record an artifact. `parents` maps parent obj_id -> relationship."""
    ids.parse(obj_id)
    conn = stores.connect("ledger")
    conn.execute(
        "INSERT INTO artifacts (obj_id, obj_type, label, store, path, produced_by,"
        " code_version, params_json, meta_json, notes, created_at)"
        " VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        (
            obj_id, obj_type, label, store, path, produced_by,
            code_version if code_version is not None else version.code_version(),
            json.dumps(params or {}), json.dumps(meta or {}), notes, ids.now_iso(),
        ),
    )
    for parent_id, rel in (parents or {}).items():
        link(parent_id, obj_id, rel, _conn=conn)
    conn.commit()
    return obj_id


def link(parent_id: str, child_id: str, relationship: str, _conn=None) -> None:
    conn = _conn or stores.connect("ledger")
    conn.execute(
        "INSERT OR IGNORE INTO artifact_links (parent_id, child_id, relationship) VALUES (?,?,?)",
        (parent_id, child_id, relationship),
    )
    if _conn is None:
        conn.commit()


def get(obj_id: str) -> Artifact:
    conn = stores.connect("ledger")
    row = conn.execute("SELECT * FROM artifacts WHERE obj_id = ?", (obj_id,)).fetchone()
    if row is None:
        raise KeyError(f"no artifact {obj_id!r}")
    return Artifact._from_row(row)


def set_meta(obj_id: str, key: str, value: Any) -> None:
    """Attach a field nobody anticipated. Promote to a column if it recurs."""
    art = get(obj_id)
    art.meta[key] = value
    conn = stores.connect("ledger")
    conn.execute("UPDATE artifacts SET meta_json = ? WHERE obj_id = ?", (json.dumps(art.meta), obj_id))
    conn.commit()


def find(obj_type: str | None = None, produced_by: str | None = None, **meta_eq) -> list[Artifact]:
    sql = "SELECT * FROM artifacts WHERE 1=1"
    args: list[Any] = []
    if obj_type:
        sql += " AND obj_type = ?"
        args.append(obj_type)
    if produced_by:
        sql += " AND produced_by = ?"
        args.append(produced_by)
    for key, val in meta_eq.items():
        sql += " AND json_extract(meta_json, ?) = ?"
        args += [f"$.{key}", val]
    sql += " ORDER BY created_at"
    conn = stores.connect("ledger")
    return [Artifact._from_row(r) for r in conn.execute(sql, args)]


def parents(obj_id: str) -> list[tuple[str, str]]:
    conn = stores.connect("ledger")
    return [(r["parent_id"], r["relationship"])
            for r in conn.execute("SELECT parent_id, relationship FROM artifact_links WHERE child_id = ?", (obj_id,))]


def children(obj_id: str) -> list[tuple[str, str]]:
    conn = stores.connect("ledger")
    return [(r["child_id"], r["relationship"])
            for r in conn.execute("SELECT child_id, relationship FROM artifact_links WHERE parent_id = ?", (obj_id,))]


def lineage(obj_id: str, direction: str = "up") -> list[str]:
    """Walk the provenance graph. 'up' = what made this, 'down' = what came from it."""
    step = parents if direction == "up" else children
    seen: list[str] = []
    stack = [obj_id]
    while stack:
        node = stack.pop()
        for nxt, _rel in step(node):
            if nxt not in seen:
                seen.append(nxt)
                stack.append(nxt)
    return seen


def supersede(old_id: str, new_id: str) -> None:
    """Record that new_id replaces old_id. Neither artifact is modified."""
    link(old_id, new_id, "superseded_by")


def is_superseded(obj_id: str) -> bool:
    return any(rel == "superseded_by" for _cid, rel in children(obj_id))


def current(obj_type: str, **meta_eq) -> list[Artifact]:
    """Artifacts of a type that nothing has replaced."""
    return [a for a in find(obj_type=obj_type, **meta_eq) if not is_superseded(a.obj_id)]
