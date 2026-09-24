"""What is actually in this lab.

A lab has the instruments it has. That is a property of the place, not of any
one run, so it belongs in a config file read once rather than an argument
threaded through every call. A run may override it -- someone borrows a
colleague's robot for an afternoon -- and the override is recorded on that
run, which is the honest way round: the exception is written down and the rule
is not restated ninety-six times.

The file is JSON with a schema, per rule 8, because eventually an agent writes
it. The agent's version of "which Integra do you have?" is a question asked
once and stored here.

    $PLATEFORGE_DATA/lab.json        (or PLATEFORGE_LAB, or a path)

Nothing here fails hard when the file is missing: an unconfigured lab gets
defaults that are honest about being defaults, so a fresh clone runs.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path

from . import paths

FILENAME = "lab.json"

DEFAULTS: dict = {
    "_comment": "What this lab has. Written once; a run may override it.",
    "name": "unconfigured",
    "liquid_handler": "generic",
    "liquid_handler_model": None,
    "has_single_channel_module": False,
    "plate_reader": None,
    "bli": None,
    "plate_format": 96,
    "default_vendor": "generic",
    "operator": None,
}


@dataclass
class Lab:
    """One lab's equipment, as declared."""
    name: str = "unconfigured"
    liquid_handler: str = "generic"
    liquid_handler_model: str | None = None
    has_single_channel_module: bool = False
    plate_reader: str | None = None
    bli: str | None = None
    plate_format: int = 96
    default_vendor: str = "generic"
    operator: str | None = None
    source: str = "defaults"
    extra: dict = field(default_factory=dict)

    @property
    def configured(self) -> bool:
        return self.source != "defaults"

    def as_dict(self) -> dict:
        return {k: v for k, v in {
            "name": self.name,
            "liquid_handler": self.liquid_handler,
            "liquid_handler_model": self.liquid_handler_model,
            "has_single_channel_module": self.has_single_channel_module,
            "plate_reader": self.plate_reader,
            "bli": self.bli,
            "plate_format": self.plate_format,
            "default_vendor": self.default_vendor,
            "operator": self.operator,
        }.items()} | self.extra

    def describe(self) -> str:
        if not self.configured:
            return ("lab: not configured -- using defaults. "
                    f"Write {path()} to set instruments.")
        lines = [f"lab: {self.name}  (from {self.source})",
                 f"  liquid handler: {self.liquid_handler}"
                 + (f" ({self.liquid_handler_model})" if self.liquid_handler_model else "")
                 + ("  + single-channel module" if self.has_single_channel_module else "")]
        for label, value in (("plate reader", self.plate_reader),
                             ("BLI", self.bli), ("vendor", self.default_vendor)):
            if value:
                lines.append(f"  {label}: {value}")
        return "\n".join(lines)


def path() -> Path:
    """Where the lab config lives. PLATEFORGE_LAB wins, else the data root."""
    override = os.environ.get("PLATEFORGE_LAB")
    return Path(override).expanduser() if override else paths.data_root() / FILENAME


def load(overrides: dict | None = None) -> Lab:
    """Read the lab config, applying any per-run overrides on top.

    An override is a deliberate exception and the caller is expected to record
    it on the run that used it -- `steps.record(..., observed={"lab": ...})`
    or the run manifest. A silent override is how a worklist ends up written
    for the wrong machine.
    """
    known = {k for k in DEFAULTS if not k.startswith("_")}
    data, source = dict(DEFAULTS), "defaults"
    where = path()
    if where.exists():
        try:
            loaded = json.loads(where.read_text())
            data.update(loaded)
            source = str(where)
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"{where} is not readable JSON: {exc}") from None
    if overrides:
        data.update(overrides)
        source = f"{source} + overrides"

    extra = {k: v for k, v in data.items()
             if k not in known and not k.startswith("_")}
    return Lab(source=source, extra=extra,
               **{k: data[k] for k in known if k in data})


def write(lab: Lab | dict, where: str | Path | None = None) -> Path:
    """Write a lab config. Creates the file the agent will later read."""
    target = Path(where) if where else path()
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = lab.as_dict() if isinstance(lab, Lab) else dict(lab)
    target.write_text(json.dumps({"_comment": DEFAULTS["_comment"]} | payload,
                                 indent=2, sort_keys=False) + "\n")
    return target


def instrument_for(job: str, lab: Lab | None = None,
                   override: str | None = None) -> str:
    """Which instrument does this job here. `override` wins and is returned.

    `job` is one of 'liquid_handler', 'plate_reader', 'bli'.
    """
    if override:
        return override
    lab = lab or load()
    value = getattr(lab, job, None)
    if not value:
        raise ValueError(
            f"this lab has no {job} configured. Write {path()} with a "
            f"{job!r} entry, or pass one explicitly for this run.")
    return value
