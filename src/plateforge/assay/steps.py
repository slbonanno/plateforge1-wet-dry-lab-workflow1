"""Protocol steps: what happened to a plate, and what came out.

A step is a function from plates to plates, registered by name. Adding one is
a new file plus a decorator (rule 5); removing one is not calling it. Nothing
downstream depends on how many steps preceded it, because the lineage lives
per well (see `assay.plates`).

## plan and record are separate

    plan(inputs, **params)  -> Plan     no side effects, re-runnable
    record(plan, **observed) -> Run     what actually happened

They are separate because they happen hours apart. You plan the ProteinA
capture at 9am, the robot and a human do it over the morning, and what gets
recorded at noon includes the lot of beads that was actually opened and the
fact that column 7 was short. A single function would have to invent that
history at planning time.

They are also separate so that a **human step has the same shape as a robot
step**. `immunoprecipitate` plans by emitting instructions and records by
taking the scientist's word for it; `harvest_supernatant` plans by emitting a
worklist and records the same way. The protocol does not fork based on who
did the work, which is what lets a manual step be automated later without
touching anything around it.

## What a step does NOT do

It does not write instrument files. A step declares what transfers it needs;
turning those into a worklist for a particular machine is `emit`'s job, and
is gated on having a real example of that machine's format (rule 6).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone

import pandas as pd

from ..core import artifacts, ids, registry, stores
from . import plates as platemod

STEPS = registry.Registry("protocol step")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Transfer:
    """One liquid movement. The unit `emit` turns into a worklist line."""
    source_plate_id: str
    source_well: str
    dest_plate_id: str
    dest_well: str
    volume_ul: float

    def as_dict(self) -> dict:
        return {"source_plate_id": self.source_plate_id,
                "source_well": self.source_well,
                "dest_plate_id": self.dest_plate_id,
                "dest_well": self.dest_well,
                "volume_ul": self.volume_ul}


@dataclass
class Plan:
    """What a step intends to do. No side effects until `commit`."""
    run_id: str
    step_type: str
    inputs: list[str]
    outputs: list[platemod.Plate]
    transfers: list[Transfer] = field(default_factory=list)
    params: dict = field(default_factory=dict)
    instructions: list[str] = field(default_factory=list)
    instrument: str | None = None
    notes: str = ""

    @property
    def transfer_frame(self) -> pd.DataFrame:
        return pd.DataFrame([t.as_dict() for t in self.transfers])

    def describe(self) -> str:
        out = [f"{self.step_type}  {self.run_id}",
               f"  in : {', '.join(self.inputs) or '(none)'}",
               f"  out: {', '.join(p.describe() for p in self.outputs)}"]
        if self.transfers:
            volumes = {t.volume_ul for t in self.transfers}
            out.append(f"  {len(self.transfers)} transfers, "
                       + (f"{volumes.pop()} uL each" if len(volumes) == 1
                          else f"{min(volumes)}-{max(volumes)} uL"))
        out += [f"  - {line}" for line in self.instructions]
        return "\n".join(out)


@dataclass
class Run:
    """A step that happened."""
    run_id: str
    step_type: str
    status: str
    outputs: list[str]
    operator: str | None = None
    lots: list[str] = field(default_factory=list)
    deviations: str = ""


class Step:
    """Base class. A subclass implements `plan`; `record` is usually enough.

    `consumes` and `produces` are content kinds, and they are checked. A step
    that expects supernatant and is handed a plate of cells should say so
    before a robot moves anything.
    """
    name: str = ""
    consumes: str | None = None
    produces: str | None = None
    robot: bool = False

    def plan(self, source: platemod.Plate, **params) -> Plan:
        raise NotImplementedError

    # -- shared machinery --------------------------------------------------
    def check_input(self, source: platemod.Plate) -> None:
        if self.consumes and source.content_kind != self.consumes:
            raise ValueError(
                f"{self.name} consumes {self.consumes!r} but "
                f"{source.plate_id} holds {source.content_kind!r}")

    def new_plan(self, inputs: list[str], outputs: list[platemod.Plate],
                 **kw) -> Plan:
        return Plan(run_id=ids.mint_stamped("RUN", self.name),
                    step_type=self.name, inputs=inputs, outputs=outputs, **kw)

    def carry_over(self, source: platemod.Plate, dest: platemod.Plate,
                   volume_ul: float | None = None,
                   roles: tuple[str, ...] = ("sample", "control")) -> list[Transfer]:
        """1:1 into the same well, recording where each well came from.

        The common case, and the reason a downstream plate keeps the layout
        of the plate that was ordered.
        """
        transfers = []
        dest_wells = dest.wells.copy()
        index = {w: i for i, w in enumerate(dest_wells["well"])}
        for _, row in source.wells.iterrows():
            if row["role"] not in roles:
                continue
            i = index.get(row["well"])
            if i is None:
                continue
            dest_wells.loc[i, "clone_id"] = row["clone_id"]
            dest_wells.loc[i, "role"] = row["role"]
            dest_wells.loc[i, "content"] = dest.content_kind
            dest_wells.loc[i, "source_plate_id"] = source.plate_id
            dest_wells.loc[i, "source_well"] = row["well"]
            if volume_ul is not None:
                dest_wells.loc[i, "volume_ul"] = volume_ul
                transfers.append(Transfer(source.plate_id, row["well"],
                                          dest.plate_id, row["well"], volume_ul))
        dest.wells = dest_wells
        return transfers


def commit(plan: Plan, parents: dict[str, str] | None = None) -> Plan:
    """Persist the plan: its plates, its transfers, and a `planned` step row."""
    conn = stores.connect("assay")
    for plate in plan.outputs:
        plate.produced_by = plan.run_id
        platemod.save(plate)

    if plan.transfers:
        xfr_id = ids.mint_stamped("XFR", plan.step_type)
        conn.executemany(
            "INSERT OR REPLACE INTO transfers (xfr_id, run_id, source_plate_id, "
            "source_well, dest_plate_id, dest_well, volume_ul, ordinal) "
            "VALUES (?,?,?,?,?,?,?,?)",
            [(xfr_id, plan.run_id, t.source_plate_id, t.source_well,
              t.dest_plate_id, t.dest_well, t.volume_ul, i)
             for i, t in enumerate(plan.transfers)])
        plan.params = dict(plan.params, xfr_id=xfr_id)

    conn.execute(
        "INSERT OR REPLACE INTO steps (run_id, step_type, status, inputs_json, "
        "outputs_json, params_json, planned_at, instrument, meta_json) "
        "VALUES (?,?,?,?,?,?,?,?,?)",
        (plan.run_id, plan.step_type, "planned", json.dumps(plan.inputs),
         json.dumps([p.plate_id for p in plan.outputs]),
         json.dumps(plan.params, sort_keys=True, default=str), _now(),
         plan.instrument, json.dumps({"instructions": plan.instructions})))
    conn.commit()

    artifacts.register(
        plan.run_id, "protocol_step", f"assay.steps.{plan.step_type}",
        label=plan.step_type, store="assay", params=plan.params,
        parents=(parents or {}) | {p: "input_to" for p in plan.inputs},
        meta={"outputs": [p.plate_id for p in plan.outputs]})
    return plan


def record(plan: Plan, *, operator: str | None = None,
           lots: list[str] | None = None, deviations: str = "",
           status: str = "done", started_at: str | None = None,
           observed: dict | None = None) -> Run:
    """Close a step out with what actually happened.

    `deviations` is a human sentence and belongs in one: "column 7 short,
    ~350 uL" is worth more later than any field this could have instead.
    """
    conn = stores.connect("assay")
    conn.execute(
        "UPDATE steps SET status = ?, operator = ?, lots_json = ?, "
        "deviations = ?, started_at = COALESCE(?, started_at), "
        "finished_at = ?, meta_json = json_patch(meta_json, ?) "
        "WHERE run_id = ?",
        (status, operator, json.dumps(lots or []), deviations, started_at,
         _now(), json.dumps({"observed": observed or {}}, default=str),
         plan.run_id))
    conn.commit()
    artifacts.set_meta(plan.run_id, "status", status)
    if deviations:
        artifacts.set_meta(plan.run_id, "deviations", deviations)
    return Run(run_id=plan.run_id, step_type=plan.step_type, status=status,
               outputs=[p.plate_id for p in plan.outputs], operator=operator,
               lots=lots or [], deviations=deviations)


def history(plate_id: str | None = None) -> pd.DataFrame:
    """Steps, newest first. With a plate id, only steps that touched it."""
    conn = stores.connect("assay")
    frame = pd.read_sql(
        "SELECT run_id, step_type, status, inputs_json, outputs_json, "
        "operator, instrument, planned_at, finished_at, deviations "
        "FROM steps ORDER BY planned_at DESC", conn)
    if plate_id and len(frame):
        touched = frame.apply(
            lambda r: plate_id in json.loads(r["inputs_json"])
            or plate_id in json.loads(r["outputs_json"]), axis=1)
        frame = frame[touched]
    return frame


def available() -> pd.DataFrame:
    """Every registered step, and what it moves between."""
    return pd.DataFrame([
        {"step": name,
         "consumes": STEPS.meta(name).get("consumes"),
         "produces": STEPS.meta(name).get("produces"),
         "robot": STEPS.meta(name).get("robot", False)}
        for name in STEPS
    ])


def register(step_class):
    """Class decorator: register a Step subclass under its own `name`."""
    STEPS.add(step_class.name, step_class(), consumes=step_class.consumes,
              produces=step_class.produces, robot=step_class.robot)
    return step_class


def get(name: str) -> Step:
    return STEPS.get(name)
