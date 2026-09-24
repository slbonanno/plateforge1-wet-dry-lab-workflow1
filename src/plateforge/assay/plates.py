"""Physical plates, and the wells in them.

The organising idea, from decision 0019: **the plate is not the sample.** A
`CLN` clone id is the through-line; plates are containers it passes through.
Every well records which clone it descends from and which well it came from.

That one choice is what makes steps addable and removable. A step is a
function from plates to plates, and nothing downstream depends on how many
steps preceded it — well D07 of an IgG prep knows it holds clone CLN-… and
came from well D07 of a supernatant plate, without either of them knowing
that a ProteinA capture happened to sit between them.

Per-well lineage rather than per-plate is deliberate. Most steps are 1:1 and
preserve the map, which is what gives "the same layout as the plate we
ordered" for free. But a rearray, a cherry-pick or a plate split is then the
same mechanism rather than a special case, and `trace()` works the same way
through all of them.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone

import pandas as pd

from ..core import artifacts, ids, stores, wells as wellmod

# What a plate holds. Not an enum: a lab will invent a content kind this list
# does not have, and a registry of strings that a step declares is better than
# a migration every time that happens.
CONTENT_KINDS = (
    "dna", "cells", "supernatant", "beads", "eluate", "igg_prep",
    "dilution", "assay",
)

WELL_ROLES = ("sample", "control", "blank", "empty")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Plate:
    """One plate, real or planned, and what is in each well."""
    plate_id: str
    barcode: str | None
    plate_format: int
    content_kind: str
    wells: pd.DataFrame
    role: str | None = None
    is_virtual: bool = False
    produced_by: str | None = None
    notes: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def occupied(self) -> pd.DataFrame:
        return self.wells[self.wells["role"] != "empty"]

    @property
    def clones(self) -> list[str]:
        got = self.wells["clone_id"].dropna()
        return list(dict.fromkeys(got))

    def well(self, position: str) -> pd.Series:
        position = wellmod.normalize(position)
        hit = self.wells[self.wells["well"] == position]
        if hit.empty:
            raise KeyError(f"{self.plate_id} has no well {position}")
        return hit.iloc[0]

    def map_of(self, column: str = "clone_id") -> dict[str, object]:
        return dict(zip(self.wells["well"], self.wells[column]))

    def describe(self) -> str:
        n = len(self.occupied)
        return (f"{self.plate_id} [{self.barcode or 'no barcode'}] "
                f"{self.content_kind}, {n}/{self.plate_format} wells"
                + (" (virtual)" if self.is_virtual else ""))


WELL_COLUMNS = ["plate_id", "well", "clone_id", "content", "role", "volume_ul",
                "concentration", "conc_units", "source_plate_id",
                "source_well", "flags", "meta_json"]


def _blank_wells(plate_id: str, plate_format: int) -> pd.DataFrame:
    rows = [{"plate_id": plate_id, "well": w, "clone_id": None, "content": None,
             "role": "empty", "volume_ul": None, "concentration": None,
             "conc_units": None, "source_plate_id": None, "source_well": None,
             "flags": None, "meta_json": "{}"}
            for w in wellmod.all_wells(plate_format, "column")]
    return pd.DataFrame(rows, columns=WELL_COLUMNS)


def create(content_kind: str, plate_format: int = 96, *,
           barcode: str | None = None, label: str = "plate",
           role: str | None = None, is_virtual: bool = False,
           produced_by: str | None = None, notes: str = "",
           parents: dict[str, str] | None = None,
           meta: dict | None = None) -> Plate:
    """Mint a new plate with every well present and empty."""
    plate_id = ids.mint_stamped("PLT", label)
    plate = Plate(plate_id=plate_id, barcode=barcode, plate_format=plate_format,
                  content_kind=content_kind, wells=_blank_wells(plate_id, plate_format),
                  role=role, is_virtual=is_virtual, produced_by=produced_by,
                  notes=notes, meta=meta or {})
    artifacts.register(
        plate_id, "plate", produced_by or "assay.plates", label=label,
        store="assay", params={"format": plate_format, "content": content_kind},
        parents=parents or {},
        meta=(meta or {}) | {"barcode": barcode, "content_kind": content_kind},
    )
    return plate


def fill(plate: Plate, rows: pd.DataFrame, *, clone_column: str = "clone_id",
         well_column: str = "well", role: str = "sample",
         volume_ul: float | None = None, content: str | None = None,
         source_plate_id: str | None = None,
         source_well_column: str | None = None) -> Plate:
    """Put clones into wells. Positions are normalised; nothing is guessed."""
    by_well = {}
    for _, row in rows.iterrows():
        position = wellmod.normalize(str(row[well_column]))
        by_well[position] = row

    updated = plate.wells.copy()
    for i, position in enumerate(updated["well"]):
        row = by_well.get(position)
        if row is None:
            continue
        updated.loc[i, "clone_id"] = row.get(clone_column)
        updated.loc[i, "role"] = role
        updated.loc[i, "content"] = content or plate.content_kind
        if volume_ul is not None:
            updated.loc[i, "volume_ul"] = volume_ul
        if source_plate_id:
            updated.loc[i, "source_plate_id"] = source_plate_id
            updated.loc[i, "source_well"] = (
                wellmod.normalize(str(row[source_well_column]))
                if source_well_column else position)
    plate.wells = updated
    return plate


def save(plate: Plate) -> str:
    """Persist a plate and its wells. Returns the plate id."""
    conn = stores.connect("assay")
    conn.execute(
        "INSERT OR REPLACE INTO plates (plate_id, barcode, format, role, "
        "content_kind, is_virtual, produced_by, meta_json, notes, created_at) "
        "VALUES (?,?,?,?,?,?,?,?,?,?)",
        (plate.plate_id, plate.barcode, plate.plate_format, plate.role,
         plate.content_kind, int(plate.is_virtual), plate.produced_by,
         json.dumps(plate.meta, sort_keys=True), plate.notes, _now()))
    conn.executemany(
        f"INSERT OR REPLACE INTO wells ({', '.join(WELL_COLUMNS)}) "
        f"VALUES ({', '.join('?' * len(WELL_COLUMNS))})",
        [tuple(None if pd.isna(v) else v for v in row)
         for row in plate.wells[WELL_COLUMNS].itertuples(index=False)])
    conn.commit()
    return plate.plate_id


def save_reads(plate_id: str, frame: pd.DataFrame, *, channel: str,
               units: str | None = None, value_column: str | None = None,
               run_id: str | None = None) -> int:
    """Store one channel of readings against a plate."""
    value_column = value_column or channel
    conn = stores.connect("assay")
    rows = [(plate_id, wellmod.normalize(str(r["well"])), channel,
             float(r[value_column]), units, run_id)
            for _, r in frame.iterrows()]
    conn.executemany(
        "INSERT OR REPLACE INTO reads (plate_id, well, channel, value, units, "
        "run_id) VALUES (?,?,?,?,?,?)", rows)
    conn.commit()
    return len(rows)


def read_values(plate_id: str, channel: str | None = None) -> pd.DataFrame:
    conn = stores.connect("assay")
    sql = ("SELECT r.well, r.channel, r.value, r.units, w.clone_id "
           "FROM reads r LEFT JOIN wells w USING (plate_id, well) "
           "WHERE r.plate_id = ?")
    params = [plate_id]
    if channel:
        sql += " AND r.channel = ?"
        params.append(channel)
    return pd.read_sql(sql, conn, params=tuple(params))


def load(plate_id: str) -> Plate:
    conn = stores.connect("assay")
    head = conn.execute(
        "SELECT barcode, format, role, content_kind, is_virtual, produced_by, "
        "meta_json, notes FROM plates WHERE plate_id = ?", (plate_id,)).fetchone()
    if head is None:
        raise KeyError(f"no plate {plate_id}")
    frame = pd.read_sql(
        f"SELECT {', '.join(WELL_COLUMNS)} FROM wells WHERE plate_id = ?",
        conn, params=(plate_id,))
    # Canonical order is column-major, the order a plate is filled and a
    # multichannel head visits it -- not the lexicographic order SQL would
    # give. A loaded plate that disagreed with the one in memory about well
    # order would make every round-trip comparison subtly wrong.
    order = {w: i for i, w in enumerate(wellmod.all_wells(head[1], "column"))}
    frame = (frame.assign(_order=frame["well"].map(order))
                  .sort_values("_order", kind="stable")
                  .drop(columns="_order")
                  .reset_index(drop=True))
    # SQLite gives back NaN where a plate built in memory holds None. Left
    # alone, every comparison downstream has to know which of the two it is
    # looking at, so the difference is removed here rather than everywhere.
    frame = frame.astype(object).where(pd.notna(frame), None)
    return Plate(plate_id=plate_id, barcode=head[0], plate_format=head[1],
                 role=head[2], content_kind=head[3], is_virtual=bool(head[4]),
                 produced_by=head[5], meta=json.loads(head[6] or "{}"),
                 notes=head[7] or "", wells=frame)


def find(content_kind: str | None = None, barcode: str | None = None) -> pd.DataFrame:
    conn = stores.connect("assay")
    sql = "SELECT plate_id, barcode, format, content_kind, produced_by, created_at FROM plates"
    clauses, params = [], []
    if content_kind:
        clauses.append("content_kind = ?")
        params.append(content_kind)
    if barcode:
        clauses.append("barcode = ?")
        params.append(barcode)
    if clauses:
        sql += " WHERE " + " AND ".join(clauses)
    return pd.read_sql(sql + " ORDER BY created_at DESC", conn, params=tuple(params))


# --- lineage ----------------------------------------------------------------

def trace(plate_id: str, well: str) -> list[dict]:
    """Every well this one descends from, newest first.

    Walks `source_plate_id`/`source_well` back to the origin. Works the same
    whether each step was 1:1 or a rearray, because the lineage is recorded
    per well rather than inferred from a plate-level mapping.
    """
    conn = stores.connect("assay")
    chain, seen = [], set()
    current = (plate_id, wellmod.normalize(well))
    while current and current not in seen:
        seen.add(current)
        row = conn.execute(
            "SELECT w.plate_id, w.well, w.clone_id, w.role, w.volume_ul, "
            "w.concentration, w.conc_units, w.source_plate_id, w.source_well, "
            "p.content_kind, p.barcode "
            "FROM wells w LEFT JOIN plates p USING (plate_id) "
            "WHERE w.plate_id = ? AND w.well = ?", current).fetchone()
        if row is None:
            break
        chain.append({
            "plate_id": row[0], "well": row[1], "clone_id": row[2],
            "role": row[3], "volume_ul": row[4], "concentration": row[5],
            "conc_units": row[6], "content_kind": row[9], "barcode": row[10],
        })
        current = (row[7], row[8]) if row[7] and row[8] else None
    return chain


def where_is(clone_id: str) -> pd.DataFrame:
    """Every well of every plate that currently holds this clone."""
    conn = stores.connect("assay")
    return pd.read_sql(
        "SELECT w.plate_id, p.barcode, p.content_kind, w.well, w.volume_ul, "
        "w.concentration, w.conc_units, p.created_at "
        "FROM wells w JOIN plates p USING (plate_id) "
        "WHERE w.clone_id = ? ORDER BY p.created_at",
        conn, params=(clone_id,))


def states_of(barcode: str) -> pd.DataFrame:
    """Every recorded state of one physical plate, oldest first.

    A barcode labels a container, not a state. Adding beads to a plate and
    then washing it produces two PLT artifacts with one barcode, because the
    artifact is immutable (rule 10) and the plastic is not.
    """
    conn = stores.connect("assay")
    return pd.read_sql(
        "SELECT plate_id, content_kind, produced_by, created_at FROM plates "
        "WHERE barcode = ? ORDER BY created_at", conn, params=(barcode,))


def layout_matches(a: Plate, b: Plate) -> bool:
    """True when both plates hold the same clone in the same well."""
    return a.map_of("clone_id") == b.map_of("clone_id")


def grid(plate: Plate, column: str = "clone_id", width: int = 12) -> str:
    """The plate as a human reads it."""
    n_rows, n_cols = wellmod.dims(plate.plate_format)
    values = {}
    for _, row in plate.wells.iterrows():
        value = row.get(column)
        text = "-" if value is None or pd.isna(value) else str(value)
        values[row["well"]] = text if len(text) <= width else text[-(width - 1):]

    out = ["     " + "".join(f"{c:>{width + 1}d}" for c in range(1, n_cols + 1))]
    for r in range(1, n_rows + 1):
        label = wellmod.row_label(r)
        cells = [f"{values.get(f'{label}{c:02d}', '-'):>{width + 1}}"
                 for c in range(1, n_cols + 1)]
        out.append(f"{label:>3}  " + "".join(cells))
    return "\n".join(out)
