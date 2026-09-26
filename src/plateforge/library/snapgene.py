"""Reading SnapGene `.dna` files, so a real plasmid map can back the vector.

Written against `fixtures/backbone/pcDNA3.1.dna` (rule 6).

Until now `library.vector` described the backbone in words -- `Element("CMV
promoter", source="backbone", notes="not ordered; sequence not tracked
here")`. That was honest and it was also a hole: without the promoter and the
5' UTR there is no way to say where a sequencing primer lands, and no way to
report a plasmid length that means anything.

A `.dna` file is a flat run of segments. One byte of type, four bytes of
big-endian length, then the payload:

    0   the DNA itself: one flags byte (bit 0 = circular), then the bases
    5   primers, as XML
    6   notes, as XML
    8   additional sequence properties, as XML
    9   the file header: b"SnapGene", then three 16-bit version fields
    10  features, as XML

Segment 3 is the enzyme recognition list, which is SnapGene's own UI state and
is skipped. Feature coordinates in the XML are 1-based and inclusive, which is
converted here exactly once so that nothing downstream has to remember it.
"""
from __future__ import annotations

import struct
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

HEADER, DNA, PRIMERS, NOTES, PROPERTIES, FEATURES = 9, 0, 5, 6, 8, 10
COOKIE = b"SnapGene"


class NotSnapGene(ValueError):
    """The file does not carry the SnapGene header segment."""


@dataclass
class Feature:
    """One annotated stretch. `start`/`end` are 0-based, end-exclusive."""
    name: str
    kind: str
    start: int
    end: int
    forward: bool = True
    notes: str = ""

    def __len__(self) -> int:
        return self.end - self.start

    def sequence(self, plasmid: "Plasmid") -> str:
        return plasmid.sequence[self.start:self.end]


@dataclass
class Plasmid:
    name: str
    sequence: str
    circular: bool = True
    features: list[Feature] = field(default_factory=list)
    description: str = ""
    created_by: str = ""

    def __len__(self) -> int:
        return len(self.sequence)

    def feature(self, name: str) -> Feature | None:
        """One feature by name, case-insensitively. None if absent."""
        wanted = name.casefold()
        return next((f for f in self.features if f.name.casefold() == wanted), None)

    def slice(self, start: int, end: int) -> str:
        """Bases from start to end, wrapping the origin if the plasmid is circular.

        A sequencing primer near the end of the numbering reads straight
        through base 1, and a slice that silently returns a short string there
        is how an alignment quietly loses its reference.
        """
        n = len(self.sequence)
        if not self.circular or 0 <= start <= end <= n:
            return self.sequence[max(start, 0):min(end, n)]
        return "".join(self.sequence[i % n] for i in range(start, end))

    def summary(self) -> dict:
        return {"name": self.name, "bp": len(self.sequence),
                "topology": "circular" if self.circular else "linear",
                "features": len(self.features),
                "created_by": self.created_by}


def segments(path: str | Path):
    """Every segment in the file, as (type, payload)."""
    data = Path(path).read_bytes()
    at = 0
    while at + 5 <= len(data):
        kind = data[at]
        (length,) = struct.unpack(">I", data[at + 1:at + 5])
        yield kind, data[at + 5:at + 5 + length]
        at += 5 + length


def _text(node, default: str = "") -> str:
    """Un-escape the HTML SnapGene wraps its free text in."""
    if node is None:
        return default
    raw = node.get("text", "") if node.tag == "V" else (node.text or "")
    cleaned = ET.fromstring(f"<x>{raw}</x>").itertext() if "<" in raw else [raw]
    return "".join(cleaned).strip() or default


def _features(payload: bytes) -> list[Feature]:
    out = []
    for element in ET.fromstring(payload.decode("utf-8")).findall("Feature"):
        spans = element.findall("Segment")
        if not spans:
            continue
        starts, ends = [], []
        for span in spans:
            first, _, last = span.get("range", "").partition("-")
            if first and last:
                starts.append(int(first) - 1)     # 1-based inclusive -> 0-based
                ends.append(int(last))
        if not starts:
            continue
        note = ""
        for qualifier in element.findall("Q"):
            if qualifier.get("name") == "note":
                note = _text(qualifier.find("V"))
                break
        out.append(Feature(
            name=element.get("name", "?"),
            kind=element.get("type", "misc_feature"),
            start=min(starts), end=max(ends),
            # SnapGene writes directionality 1 forward, 2 reverse, 3 both.
            forward=element.get("directionality", "1") != "2",
            notes=note))
    return sorted(out, key=lambda f: (f.start, f.end))


def _notes(payload: bytes) -> dict:
    root = ET.fromstring(payload.decode("utf-8"))
    return {"name": _text(root.find("CustomMapLabel")),
            "description": _text(root.find("Description")),
            "created_by": _text(root.find("CreatedBy"))}


def read(path: str | Path, name: str = "") -> Plasmid:
    """One `.dna` file as a plasmid with its features."""
    path = Path(path)
    sequence, circular, features, notes, seen_header = "", True, [], {}, False

    for kind, payload in segments(path):
        if kind == HEADER:
            if not payload.startswith(COOKIE):
                raise NotSnapGene(f"{path.name}: header segment is not SnapGene")
            seen_header = True
        elif kind == DNA and payload:
            circular = bool(payload[0] & 1)
            sequence = payload[1:].decode("ascii").upper()
        elif kind == FEATURES:
            features = _features(payload)
        elif kind == NOTES:
            notes = _notes(payload)

    if not seen_header:
        raise NotSnapGene(f"{path.name}: no SnapGene header segment; not a .dna file")
    if not sequence:
        raise NotSnapGene(f"{path.name}: no DNA segment")

    return Plasmid(name=name or notes.get("name") or path.stem,
                   sequence=sequence, circular=circular, features=features,
                   description=notes.get("description", ""),
                   created_by=notes.get("created_by", ""))


def check(plasmid: Plasmid, expect_bp: int | None = None,
          expect_features: tuple[str, ...] = ()) -> list[str]:
    """What is wrong with a plasmid we are about to build on.

    A backbone is the one input where a silent error is unrecoverable: a
    one-base slip shifts every ORF built on it, survives every test that only
    looks at the insert, and shows up months later as an expression run that
    produced nothing. So the file gets checked against what it is supposed to
    be before anything uses it, and the caller is told rather than guessed at.
    """
    problems = []
    if expect_bp is not None and len(plasmid) != expect_bp:
        problems.append(f"expected {expect_bp} bp, found {len(plasmid)}")
    stray = set(plasmid.sequence) - set("ACGT")
    if stray:
        problems.append(f"non-ACGT bases present: {''.join(sorted(stray))}")
    for wanted in expect_features:
        if plasmid.feature(wanted) is None:
            problems.append(f"no feature named {wanted!r}")
    return problems
