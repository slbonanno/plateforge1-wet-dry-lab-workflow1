"""Instrument worklists: turning planned transfers into a file a machine reads.

A step declares the transfers it needs (`assay.steps.Transfer`). This turns
those into a file for one named instrument. The split matters: the protocol
does not know what is in the lab, and the instrument does not know why the
liquid is moving.

## Rule 6 applies here more than anywhere

A worklist emitted against a guessed format is a wasted plate of reagent, and
possibly a wasted week. So every emitter declares how well its format is
actually known:

    verified    a real file from this instrument is in fixtures/, and the
                emitter's test reads it. Nothing is verified yet.
    documented  the format is specified in public documentation from the
                vendor or a mature open-source implementation, and at least
                two independent sources agree on the field order. Emitted on
                request, stamped, never silently.
    sketch      the shape is known but the details are not. Refuses to emit
                and says what would unblock it.

`emit(..., allow_unverified=True)` is required for anything below `verified`,
and every such file is written with a companion `.CAVEAT.txt`.

## What is deliberately missing

Nothing here drives an instrument. These are files a human loads, reviews and
runs. That is not a limitation to be removed later -- a review step between a
computed volume and a moving pipette is worth keeping.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..core import registry, wells as wellmod

EMITTERS = registry.Registry("instrument worklist")

VERIFIED, DOCUMENTED, SKETCH = "verified", "documented", "sketch"


@dataclass
class Worklist:
    """A file for one instrument, and everything that qualifies it."""
    instrument: str
    software: str
    text: str
    suffix: str
    confidence: str
    lines: int = 0
    caveats: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)
    limits: dict = field(default_factory=dict)
    issues: list[str] = field(default_factory=list)

    @property
    def verified(self) -> bool:
        return self.confidence == VERIFIED

    def write(self, directory: str | Path, stem: str = "worklist") -> Path:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{stem}_{self.instrument}{self.suffix}"
        path.write_text(self.text)
        if not self.verified:
            path.with_suffix(".CAVEAT.txt").write_text(
                f"{self.instrument} / {self.software}\n"
                f"format confidence: {self.confidence}\n\n"
                + "\n".join(f"- {c}" for c in self.caveats)
                + "\n\nsources:\n"
                + "\n".join(f"  {s}" for s in self.sources) + "\n")
        return path

    def as_dict(self) -> dict:
        return {"instrument": self.instrument, "software": self.software,
                "confidence": self.confidence, "lines": self.lines,
                "suffix": self.suffix, "caveats": self.caveats,
                "issues": self.issues, "limits": self.limits}


# --- helpers ----------------------------------------------------------------

def _frame(transfers) -> pd.DataFrame:
    """Accept Transfer objects or a frame; return a frame."""
    if isinstance(transfers, pd.DataFrame):
        return transfers
    return pd.DataFrame([t.as_dict() for t in transfers])


def tecan_position(well: str, plate_format: int = 96) -> int:
    """Tecan numbers wells down each column, 1-based. A01=1, B01=2, A02=9."""
    well = wellmod.normalize(well)
    n_rows, _ = wellmod.dims(plate_format)
    row = ord(well[0]) - ord("A")
    col = int(well[1:]) - 1
    return col * n_rows + row + 1


def check_volumes(frame: pd.DataFrame, minimum: float, maximum: float) -> list[str]:
    out = []
    if frame.empty:
        return ["no transfers to emit"]
    low = frame[frame["volume_ul"] < minimum]
    high = frame[frame["volume_ul"] > maximum]
    if len(low):
        out.append(f"{len(low)} transfer(s) below the {minimum} uL minimum "
                   f"(smallest {low['volume_ul'].min()})")
    if len(high):
        out.append(f"{len(high)} transfer(s) above the {maximum} uL maximum "
                   f"(largest {high['volume_ul'].max()}); they need splitting")
    return out


# --- ours -------------------------------------------------------------------

@EMITTERS.register("generic", software="none", confidence=VERIFIED,
                   suffix=".csv")
def generic(transfers, *, source_label: str = "Source",
            dest_label: str = "Destination", plate_format: int = 96,
            **_) -> Worklist:
    """Our own columns. Nothing here is a guess about anyone's format."""
    frame = _frame(transfers)
    text = frame.to_csv(index=False) if len(frame) else "source_well,dest_well,volume_ul\n"
    return Worklist(instrument="generic", software="none", text=text,
                    suffix=".csv", confidence=VERIFIED, lines=len(frame))


# --- documented -------------------------------------------------------------

TECAN_SOURCES = [
    "https://robotools.readthedocs.io/en/latest/notebooks/02_Worklist_Basics.html",
    "https://github.com/Edinburgh-Genome-Foundry/Dioscuri",
]


@EMITTERS.register("tecan_evo", software="Freedom EVOware / FluentControl",
                   confidence=DOCUMENTED, suffix=".gwl")
def tecan_evo(transfers, *, source_label: str = "Source",
              dest_label: str = "Destination", plate_format: int = 96,
              liquid_class: str = "", wash_between: bool = True,
              min_ul: float = 1.0, max_ul: float = 950.0, **_) -> Worklist:
    """Tecan Gemini WorkList.

    One aspirate, one dispense, one wash per transfer. Records are
    semicolon-separated with ten fields after the letter:

        A;rack label;rack id;rack type;position;tube id;volume;
          liquid class;tip type;tip mask;forced rack type

    Positions are 1-based and run down each column, so A01=1 and A02=9 on a
    96-well plate. Two independent implementations agree on this layout; no
    real file from an instrument has been checked against it.
    """
    frame = _frame(transfers)
    lines = []
    for _, row in frame.iterrows():
        source = tecan_position(row["source_well"], plate_format)
        dest = tecan_position(row["dest_well"], plate_format)
        volume = f"{float(row['volume_ul']):.2f}"
        lines.append(f"A;{source_label};;;{source};;{volume};{liquid_class};;;")
        lines.append(f"D;{dest_label};;;{dest};;{volume};{liquid_class};;;")
        if wash_between:
            lines.append("W1;")
    return Worklist(
        instrument="tecan_evo", software="Freedom EVOware / FluentControl",
        text="\n".join(lines) + ("\n" if lines else ""), suffix=".gwl",
        confidence=DOCUMENTED, lines=len(lines),
        limits={"min_ul": min_ul, "max_ul": max_ul},
        issues=check_volumes(frame, min_ul, max_ul),
        sources=TECAN_SOURCES,
        caveats=[
            "Field order is taken from public documentation and an "
            "open-source implementation, not from a file this instrument "
            "produced. Check the first line against a known-good worklist.",
            "Rack labels must match the labware names in the EVOware/Fluent "
            "deck layout; they are passed through verbatim.",
            "Liquid class is left empty, so the instrument's default applies.",
        ])


OPENTRONS_SOURCES = ["https://docs.opentrons.com/v2/"]


@EMITTERS.register("opentrons", software="Opentrons Python Protocol API v2",
                   confidence=DOCUMENTED, suffix=".py")
def opentrons(transfers, *, source_label: str = "source",
              dest_label: str = "destination", plate_format: int = 96,
              pipette: str = "p300_single_gen2",
              tiprack: str = "opentrons_96_tiprack_300ul",
              labware: str = "corning_96_wellplate_360ul_flat",
              api_version: str = "2.14", **_) -> Worklist:
    """An Opentrons protocol, as Python.

    Unusual among these: the "worklist" is a program, and the API is public
    and versioned, so the shape is not in doubt. What IS in doubt is the deck
    -- slot numbers and labware names are this emitter's guesses and must be
    matched to the actual setup.
    """
    frame = _frame(transfers)
    body = [
        '"""Generated by plateforge. Review the deck layout before running."""',
        "from opentrons import protocol_api",
        "",
        f'metadata = {{"protocolName": "plateforge transfer",',
        f'            "apiLevel": "{api_version}"}}',
        "",
        "",
        "def run(protocol: protocol_api.ProtocolContext):",
        f'    tips = protocol.load_labware("{tiprack}", 1)',
        f'    source = protocol.load_labware("{labware}", 2)',
        f'    destination = protocol.load_labware("{labware}", 3)',
        f'    pipette = protocol.load_instrument("{pipette}", "right",',
        "                                       tip_racks=[tips])",
        "",
    ]
    for _, row in frame.iterrows():
        body.append(f'    pipette.transfer({float(row["volume_ul"]):.2f}, '
                    f'source["{wellmod.normalize(row["source_well"])}"], '
                    f'destination["{wellmod.normalize(row["dest_well"])}"])')
    if frame.empty:
        body.append("    pass")
    return Worklist(
        instrument="opentrons", software="Opentrons Python Protocol API v2",
        text="\n".join(body) + "\n", suffix=".py", confidence=DOCUMENTED,
        lines=len(frame), sources=OPENTRONS_SOURCES,
        issues=check_volumes(frame, 1.0, 300.0),
        limits={"min_ul": 1.0, "max_ul": 300.0},
        caveats=[
            "Deck slots (1 tips, 2 source, 3 destination) and labware names "
            "are placeholders. They must match the real deck.",
            "Wells are addressed as A01; Opentrons also accepts A1. If the "
            "labware definition disagrees, fix it here rather than in the file.",
            "Simulate with `opentrons_simulate` before running.",
        ])


INTEGRA_SOURCES = [
    "https://www.integra-biosciences.com/sites/default/files/documents/128951_V09_OI_VIALAB.pdf",
]


@EMITTERS.register("integra_vialab", software="VIALAB (D-ONE worklist import)",
                   confidence=SKETCH, suffix=".csv")
def integra_vialab(transfers, *, plate_format: int = 96,
                   worklist_kind: str = "worklist", **_) -> Worklist:
    """INTEGRA VIALAB worklist import.

    VIALAB imports three kinds of CSV -- Hit Picking, Normalization Worklist
    and Worklist -- semicolon separated, at most 384 lines. All three are
    **D-ONE module** features; a VIAFLO multichannel head has no documented
    CSV route and its variable volumes are set in the GUI.

    The column layout is not known. INTEGRA ships Excel templates that
    generate compliant files and they are behind an account login, so this
    emitter refuses until one is in `fixtures/integra/`.
    """
    frame = _frame(transfers)
    raise NotImplementedError(
        "integra_vialab has no real template to build against (Q15). "
        "VIALAB's CSV import is a D-ONE feature and INTEGRA's templates are "
        "behind a login. Put one in fixtures/integra/ -- see the README there "
        f"for what to capture -- and this becomes a two-hour job. "
        f"({len(frame)} transfers are waiting for it.)")


# --- sketches ---------------------------------------------------------------

def _closed_format(name: str, software: str, why: str, unblock: str):
    def emitter(transfers, **_) -> Worklist:
        raise NotImplementedError(f"{name}: {why} {unblock}")
    emitter.__name__ = name
    emitter.__doc__ = f"{software}. Not implemented: {why}"
    return EMITTERS.register(name, software=software, confidence=SKETCH,
                             suffix=".txt")(emitter)


_closed_format(
    "hamilton_venus", "Hamilton VENUS",
    "VENUS methods are a binary/proprietary project format, not a worklist "
    "file, so there is nothing to generate directly.",
    "PyHamilton drives VENUS from Python instead; that is an integration, not "
    "an emitter, and needs the instrument in the loop to develop against.")

_closed_format(
    "beckman_biomek", "Beckman Biomek",
    "Biomek's transfer-file format is not publicly specified.",
    "A real exported transfer file in fixtures/beckman/ would settle it.")

_closed_format(
    "agilent_bravo", "Agilent Bravo / VWorks",
    "VWorks protocols are XML but the schema is not published.",
    "A real .pro file in fixtures/agilent/ would settle it.")


# --- the front door ---------------------------------------------------------

def emit(transfers, instrument: str = "generic",
         allow_unverified: bool = False, **kwargs) -> Worklist:
    """Build one instrument's worklist.

    Anything short of `verified` has to be asked for. The default is our own
    format, which cannot be wrong about anyone else's.
    """
    meta = EMITTERS.meta(instrument)
    confidence = meta.get("confidence", SKETCH)
    if confidence != VERIFIED and not allow_unverified:
        raise ValueError(
            f"{instrument!r} format confidence is {confidence!r}: no file from "
            "this instrument has been checked against the emitter. Pass "
            "allow_unverified=True to emit it anyway (it will be stamped), or "
            "use instrument='generic'.")
    return EMITTERS.get(instrument)(transfers, **kwargs)


def available() -> pd.DataFrame:
    """Every registered instrument and how well its format is known."""
    return pd.DataFrame([
        {"instrument": name,
         "software": EMITTERS.meta(name).get("software"),
         "confidence": EMITTERS.meta(name).get("confidence"),
         "suffix": EMITTERS.meta(name).get("suffix")}
        for name in EMITTERS
    ]).sort_values(["confidence", "instrument"]).reset_index(drop=True)
