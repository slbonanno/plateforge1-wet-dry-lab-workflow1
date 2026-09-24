"""The steps themselves: transfection through ELISA input.

One class per step, each registered. They are separate from `steps.py` on
purpose — that file is the mechanism and this file is one lab's protocol. A
different lab's protocol is a different file, not a fork of the machinery.

Volumes and incubation times are **parameters with declared defaults**, never
hardcoded constants, because they are wet-lab decisions (Q19). Every default
here is a placeholder that says so in `assumed`; a campaign supplies real
numbers through a protocol JSON (rule 8).

Nothing here writes an instrument file. A step says what transfers it needs;
`emit` turns those into a worklist for a named machine, and only once a real
example of that machine's format exists (rule 6).
"""
from __future__ import annotations

from . import plates as platemod
from .steps import Plan, Step, Transfer, register


# --- helpers ---------------------------------------------------------------

def _defaults(**kw):
    """Mark a step's numeric defaults as placeholders, not decisions."""
    return kw


# --- the chain -------------------------------------------------------------

@register
class Transfect(Step):
    """DNA plate in, culture plate out. Same layout, new container."""
    name = "transfect"
    consumes = "dna"
    produces = "cells"
    robot = False
    assumed = _defaults(dna_ng_per_well=200.0, cells_per_well=1e5,
                        volume_ul=200.0, days=5)

    def plan(self, source: platemod.Plate, *, volume_ul: float = 200.0,
             days: int = 5, barcode: str | None = None, **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("cells", source.plate_format, barcode=barcode,
                               label="culture")
        self.carry_over(source, dest, volume_ul=None)
        for i in dest.wells.index:
            if dest.wells.loc[i, "role"] != "empty":
                dest.wells.loc[i, "volume_ul"] = volume_ul
        return self.new_plan(
            [source.plate_id], [dest], params={"volume_ul": volume_ul, "days": days} | params,
            instructions=[
                f"Transfect {len(source.occupied)} wells from {source.plate_id}.",
                f"Express {days} days before harvest.",
            ])


@register
class HarvestSupernatant(Step):
    """Draw supernatant off the culture into a fresh plate. 1:1."""
    name = "harvest_supernatant"
    consumes = "cells"
    produces = "supernatant"
    robot = True
    assumed = _defaults(volume_ul=500.0)

    def plan(self, source: platemod.Plate, *, volume_ul: float = 500.0,
             barcode: str | None = None, **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("supernatant", source.plate_format,
                               barcode=barcode, label="supernatant")
        transfers = self.carry_over(source, dest, volume_ul=volume_ul)
        return self.new_plan(
            [source.plate_id], [dest], transfers=transfers,
            params={"volume_ul": volume_ul} | params, instrument="liquid_handler",
            instructions=[f"Transfer {volume_ul} uL supernatant, 1:1, "
                          f"{source.plate_id} -> {dest.plate_id}.",
                          "Avoid the cell layer; do not touch the bottom."])


@register
class AddBeads(Step):
    """Add pre-washed ProteinA Dynabeads to every occupied well.

    A reagent addition, not a plate-to-plate transfer: one trough, 96
    destinations. Stays on the same plate -- the beads join the supernatant
    rather than moving it -- so the output is the same container with a new
    content kind.
    """
    name = "add_beads"
    consumes = "supernatant"
    produces = "beads"
    robot = True
    assumed = _defaults(bead_slurry_ul=25.0, incubate_minutes=60,
                        temperature_c=4)

    def plan(self, source: platemod.Plate, *, bead_slurry_ul: float = 25.0,
             incubate_minutes: int = 60, temperature_c: int = 4,
             barcode: str | None = None, **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("beads", source.plate_format,
                               barcode=barcode or source.barcode,
                               label="capture")
        self.carry_over(source, dest, volume_ul=None)
        for i in dest.wells.index:
            if dest.wells.loc[i, "role"] != "empty":
                existing = source.wells.loc[i, "volume_ul"] or 0.0
                dest.wells.loc[i, "volume_ul"] = existing + bead_slurry_ul
        return self.new_plan(
            [source.plate_id], [dest],
            params={"bead_slurry_ul": bead_slurry_ul,
                    "incubate_minutes": incubate_minutes,
                    "temperature_c": temperature_c} | params,
            instrument="liquid_handler",
            instructions=[
                f"Pre-wash ProteinA Dynabeads; dispense {bead_slurry_ul} uL "
                f"to each of {len(source.occupied)} wells from one trough.",
                f"Incubate {incubate_minutes} min at {temperature_c} C with mixing.",
            ])


@register
class Immunoprecipitate(Step):
    """Capture and wash. A human does this one; it records the same way."""
    name = "immunoprecipitate"
    consumes = "beads"
    produces = "beads"
    robot = False
    assumed = _defaults(washes=3, wash_volume_ul=200.0)

    def plan(self, source: platemod.Plate, *, washes: int = 3,
             wash_volume_ul: float = 200.0, barcode: str | None = None,
             **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("beads", source.plate_format,
                               barcode=barcode or source.barcode, label="washed")
        self.carry_over(source, dest, volume_ul=None)
        return self.new_plan(
            [source.plate_id], [dest],
            params={"washes": washes, "wash_volume_ul": wash_volume_ul} | params,
            instructions=[
                "Magnet; remove depleted supernatant.",
                f"Wash {washes}x with {wash_volume_ul} uL, magnet between each.",
                "Record any well that ran dry or lost beads as a deviation.",
            ])


@register
class ElutePlate(Step):
    """Low-pH elution and immediate neutralisation, on the bead plate."""
    name = "elute_neutralise"
    consumes = "beads"
    produces = "eluate"
    robot = True
    assumed = _defaults(elution_ul=50.0, neutralisation_ul=5.0,
                        incubate_minutes=5)

    def plan(self, source: platemod.Plate, *, elution_ul: float = 50.0,
             neutralisation_ul: float = 5.0, incubate_minutes: int = 5,
             barcode: str | None = None, **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("eluate", source.plate_format,
                               barcode=barcode or source.barcode, label="eluate")
        self.carry_over(source, dest, volume_ul=None)
        total = elution_ul + neutralisation_ul
        for i in dest.wells.index:
            if dest.wells.loc[i, "role"] != "empty":
                dest.wells.loc[i, "volume_ul"] = total
        return self.new_plan(
            [source.plate_id], [dest],
            params={"elution_ul": elution_ul,
                    "neutralisation_ul": neutralisation_ul,
                    "incubate_minutes": incubate_minutes} | params,
            instrument="liquid_handler",
            instructions=[
                f"Dispense {elution_ul} uL low-pH elution buffer to each well.",
                f"Incubate {incubate_minutes} min with mixing.",
                f"Add {neutralisation_ul} uL neutralisation buffer to each well "
                "BEFORE removing from the beads -- IgG does not enjoy the wait.",
            ])


@register
class MagnetTransfer(Step):
    """Magnet, then move the neutralised eluate off the beads. The IgG prep.

    This is the plate that gets a new real-life barcode and carries the same
    sample map as the plate that was ordered -- which `plates.layout_matches`
    can assert rather than anyone having to trust it.
    """
    name = "magnet_transfer"
    consumes = "eluate"
    produces = "igg_prep"
    robot = True
    assumed = _defaults(volume_ul=50.0)

    def plan(self, source: platemod.Plate, *, volume_ul: float = 50.0,
             barcode: str | None = None, **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("igg_prep", source.plate_format,
                               barcode=barcode, label="igg")
        transfers = self.carry_over(source, dest, volume_ul=volume_ul)
        return self.new_plan(
            [source.plate_id], [dest], transfers=transfers,
            params={"volume_ul": volume_ul} | params,
            instrument="liquid_handler",
            instructions=[
                "Magnet until the supernatant is clear.",
                f"Transfer {volume_ul} uL, 1:1, to {dest.plate_id}. "
                "Do not carry beads over.",
            ])


@register
class QuantifyBLI(Step):
    """Read IgG concentration. The plate is the input; the export comes back.

    Plans a read rather than a transfer: nothing moves, so there are no
    transfers, only a sample map the instrument is set up against. The export
    is ingested by `record(observed=...)` once it exists (Q7).
    """
    name = "quantify_bli"
    consumes = "igg_prep"
    produces = "igg_prep"
    robot = True
    assumed = _defaults(dilution=1.0, sensor="ProteinA")

    def plan(self, source: platemod.Plate, *, sensor: str = "ProteinA",
             dilution: float = 1.0, **params) -> Plan:
        self.check_input(source)
        return self.new_plan(
            [source.plate_id], [], instrument="bli",
            params={"sensor": sensor, "dilution": dilution} | params,
            instructions=[
                f"Read {source.plate_id} on {sensor} sensors, 96-well format.",
                "Include an IgG standard curve.",
                "Nothing is transferred; the plate is read in place.",
            ],
            notes="Concentrations land on the same plate via record(observed=).")


@register
class Normalise(Step):
    """Dilute every well to one concentration. The variable-volume step.

    This is the only step whose volumes differ per well, which is exactly the
    capability a liquid handler either has or does not. `emit` decides how to
    express it; the step only says what is required.
    """
    name = "normalise"
    consumes = "igg_prep"
    produces = "dilution"
    robot = True
    assumed = _defaults(target_conc=10.0, target_volume_ul=100.0,
                        min_transfer_ul=2.0)

    def plan(self, source: platemod.Plate, *, target_conc: float = 10.0,
             target_volume_ul: float = 100.0, min_transfer_ul: float = 2.0,
             barcode: str | None = None, **params) -> Plan:
        self.check_input(source)
        dest = platemod.create("dilution", source.plate_format,
                               barcode=barcode, label="normalised")
        self.carry_over(source, dest, volume_ul=None)

        transfers, short = [], []
        index = {w: i for i, w in enumerate(dest.wells["well"])}
        for _, row in source.wells.iterrows():
            if row["role"] == "empty":
                continue
            have = row["concentration"]
            i = index[row["well"]]
            if have is None or (isinstance(have, float) and have != have):
                short.append((row["well"], "no concentration measured"))
                continue
            need = target_conc * target_volume_ul / have if have else 0.0
            if need < min_transfer_ul:
                short.append((row["well"], f"{need:.1f} uL below the "
                                           f"{min_transfer_ul} uL minimum"))
                need = min_transfer_ul
            if need > target_volume_ul:
                short.append((row["well"], f"too dilute: needs {need:.0f} uL "
                                           f"into {target_volume_ul} uL"))
                need = target_volume_ul
            dest.wells.loc[i, "volume_ul"] = target_volume_ul
            dest.wells.loc[i, "concentration"] = min(
                target_conc, have * need / target_volume_ul)
            dest.wells.loc[i, "conc_units"] = row["conc_units"] or "ug/mL"
            if short and short[-1][0] == row["well"]:
                dest.wells.loc[i, "flags"] = short[-1][1]
            transfers.append(Transfer(source.plate_id, row["well"],
                                      dest.plate_id, row["well"],
                                      round(need, 2)))

        return self.new_plan(
            [source.plate_id], [dest], transfers=transfers,
            params={"target_conc": target_conc,
                    "target_volume_ul": target_volume_ul,
                    "min_transfer_ul": min_transfer_ul,
                    "wells_flagged": len(short)} | params,
            instrument="liquid_handler",
            instructions=[
                f"Dilute each well to {target_conc} ug/mL in "
                f"{target_volume_ul} uL. Volumes differ per well.",
            ] + ([f"{len(short)} well(s) could not hit the target; see flags."]
                 if short else []))


# The order the steps are meant to run in. A protocol JSON may name a subset
# or insert its own; this is the default chain, not a constraint.
DEFAULT_CHAIN = ["transfect", "harvest_supernatant", "add_beads",
                 "immunoprecipitate", "elute_neutralise", "magnet_transfer",
                 "quantify_bli", "normalise"]
