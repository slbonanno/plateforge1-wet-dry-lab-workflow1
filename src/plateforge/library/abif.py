"""Reading `.ab1` trace files (ABIF), the format every Sanger vendor returns.

Written against three real traces from Azenta/GENEWIZ in `fixtures/sanger/`
(rule 6), not against the spec alone. What those files turned out to contain
changed the design of everything downstream, so it is worth stating plainly:

**The trace knows where it came from.** Alongside the basecalls it carries the
well it was run from (`TUBE`), the plate barcode and name (`CTID`, `CTNM`),
the run start date and time (`RUND1`/`RUNT1`), the instrument model, the dye
set, and the vendor's own name (`User`). All of it written by the sequencer.

That matters because the alternative -- parsing a filename -- is guesswork,
and because file modification times are rewritten by downloading. A rerun and
its original can be ordered correctly from inside the files even when the
filenames are identical and the download order was backwards.

## Format

A 4-byte magic `ABIF`, a version, then a directory entry at byte 6 that
describes the directory itself: 28 bytes, big-endian, giving the number of
entries and where they start. Every entry is another 28 bytes:

    name (4 bytes)  number (int32)  element type (int16)  element size (int16)
    n elements (int32)  data size (int32)  data offset (int32)  handle (int32)

A tag is (name, number), not name alone -- `RUND` appears four times in a real
file and they mean different things. Data of 4 bytes or fewer is stored inline
in the offset field rather than at the offset, which is the detail that
silently returns garbage if you miss it.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from datetime import date, datetime, time
from pathlib import Path

MAGIC = b"ABIF"
ENTRY = struct.Struct(">4slhhllll")
ENTRY_SIZE = 28

# ABIF element types. The ones this module actually decodes; anything else
# comes back as raw bytes rather than as a wrong guess.
BYTE, CHAR, WORD, SHORT, LONG = 1, 2, 3, 4, 5
FLOAT, DOUBLE, DATE, TIME, PSTRING, CSTRING = 7, 8, 10, 11, 18, 19

# Tags worth naming. Everything else stays reachable through `raw`.
SEQUENCE, QUALITY = ("PBAS", 1), ("PCON", 1)


class NotABIF(ValueError):
    """The file does not start with the ABIF magic."""


@dataclass
class Trace:
    """One `.ab1` file, read.

    `sample`, `well` and `plate` are what the instrument recorded, not what a
    filename claims. `path` is kept so a report can point at the file.
    """
    path: Path
    sequence: str
    quality: list[int]
    sample: str = ""
    well: str = ""
    plate_id: str = ""
    plate_name: str = ""
    instrument: str = ""
    dye_set: str = ""
    vendor: str = ""
    run_started: datetime | None = None
    raw: dict = field(default_factory=dict, repr=False)

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def mean_quality(self) -> float:
        return sum(self.quality) / len(self.quality) if self.quality else 0.0

    def bases_at_least(self, phred: int = 20) -> int:
        return sum(1 for q in self.quality if q >= phred)

    def tag(self, name: str, number: int = 1):
        """Any ABIF tag, decoded if the type is one we know."""
        return self.raw.get((name, number))

    def summary(self) -> dict:
        return {
            "file": self.path.name, "sample": self.sample, "well": self.well,
            "plate_id": self.plate_id, "bases": len(self.sequence),
            "mean_quality": round(self.mean_quality, 1),
            "q20_bases": self.bases_at_least(20),
            "run_started": self.run_started.isoformat() if self.run_started else None,
            "instrument": self.instrument, "vendor": self.vendor,
        }


def _decode(kind: int, n_elements: int, payload: bytes):
    """One entry's payload, by ABIF element type."""
    if kind == CHAR:
        return payload[:n_elements]
    if kind == PSTRING:                       # Pascal string: length byte first
        return payload[1:1 + payload[0]].decode("latin-1") if payload else ""
    if kind == CSTRING:                       # C string: trailing NUL
        return payload[:max(n_elements - 1, 0)].decode("latin-1")
    if kind == BYTE:
        return list(payload[:n_elements])
    if kind in (SHORT, WORD):
        code = "h" if kind == SHORT else "H"
        return list(struct.unpack(f">{n_elements}{code}", payload[:2 * n_elements]))
    if kind == LONG:
        return list(struct.unpack(f">{n_elements}l", payload[:4 * n_elements]))
    if kind == FLOAT:
        return list(struct.unpack(f">{n_elements}f", payload[:4 * n_elements]))
    if kind == DATE:
        year, month, day = struct.unpack(">hBB", payload[:4])
        try:
            return date(year, month, day)
        except ValueError:
            return None
    if kind == TIME:
        hour, minute, second, _hundredths = struct.unpack(">4B", payload[:4])
        try:
            return time(hour, minute, second)
        except ValueError:
            return None
    return payload


def entries(path: str | Path) -> dict:
    """Every tag in the file, as {(name, number): value}.

    The whole directory, decoded. `read()` is the useful entry point; this is
    here because a vendor will eventually put something we care about in a tag
    nobody has looked at yet, and finding it should not need a new parser.
    """
    path = Path(path)
    data = path.read_bytes()
    if data[:4] != MAGIC:
        raise NotABIF(f"{path.name} does not start with {MAGIC!r}; not an .ab1 file")

    *_, count, _size, offset, _handle = ENTRY.unpack(data[6:6 + ENTRY_SIZE])
    out: dict = {}
    for index in range(count):
        at = offset + index * ENTRY_SIZE
        chunk = data[at:at + ENTRY_SIZE]
        if len(chunk) < ENTRY_SIZE:
            break                              # truncated directory; keep what parsed
        name, number, kind, _esize, n, dsize, doffset, _h = ENTRY.unpack(chunk)
        # Four bytes or fewer live in the offset field itself, not at it.
        payload = (data[at + 20:at + 24] if dsize <= 4
                   else data[doffset:doffset + dsize])
        out[(name.decode("latin-1"), number)] = _decode(kind, n, payload)
    return out


def _text(tags: dict, name: str, number: int = 1) -> str:
    value = tags.get((name, number))
    if isinstance(value, bytes):
        return value.decode("latin-1").strip()
    return str(value).strip() if value is not None else ""


def read(path: str | Path) -> Trace:
    """One trace, with the fields anything downstream actually needs."""
    path = Path(path)
    tags = entries(path)

    bases = tags.get(SEQUENCE, b"")
    sequence = bases.decode("latin-1") if isinstance(bases, bytes) else str(bases)
    # PCON is a char array whose *bytes* are Phred scores, not text. Decoding
    # it as a string is the obvious mistake and yields quality values in the
    # 30s-60s for every base, which looks plausible and is meaningless.
    scores = tags.get(QUALITY, b"")
    quality = list(scores) if isinstance(scores, (bytes, bytearray)) else list(scores or [])

    started = None
    day, clock = tags.get(("RUND", 1)), tags.get(("RUNT", 1))
    if isinstance(day, date) and isinstance(clock, time):
        started = datetime.combine(day, clock)
    elif isinstance(day, date):
        started = datetime.combine(day, time())

    return Trace(
        path=path,
        sequence=sequence,
        quality=quality,
        sample=_text(tags, "SMPL"),
        well=_text(tags, "TUBE"),
        plate_id=_text(tags, "CTID"),
        plate_name=_text(tags, "CTNM"),
        instrument=_text(tags, "MODL"),
        dye_set=_text(tags, "DySN"),
        vendor=_text(tags, "User"),
        run_started=started,
        raw=tags,
    )


def read_all(paths) -> list[Trace]:
    """Several traces. A file that is not ABIF is reported, not skipped."""
    traces, problems = [], []
    for path in paths:
        try:
            traces.append(read(path))
        except NotABIF as exc:
            problems.append(str(exc))
    if problems:
        raise NotABIF("; ".join(problems))
    return traces


def mott_trim(quality: list[int], cutoff: int = 20,
              minimum: int = 20) -> tuple[int, int]:
    """Where the good part of the read starts and ends (Mott's algorithm).

    The standard trim, and the reason to do it before aligning rather than
    after: untrimmed 5' noise aligns to *something*, and produces confident
    mismatches that read as real mutations. The three fixture traces start
    with 13-23 called Ns apiece, so this is not a hypothetical.

    Each base scores `p_error_cutoff - p_error(base)`; the returned window is
    the highest-scoring run. Returns `(0, 0)` when nothing clears the cutoff,
    which is a failed read and should be reported as one.
    """
    if not quality:
        return 0, 0
    limit = 10 ** (-cutoff / 10.0)
    best_score = best_start = best_end = 0
    running = 0.0
    start = 0
    for index, phred in enumerate(quality):
        running += limit - 10 ** (-phred / 10.0)
        if running < 0:
            running, start = 0.0, index + 1
        elif running > best_score:
            best_score, best_start, best_end = running, start, index + 1
    if best_end - best_start < minimum:
        return 0, 0
    return best_start, best_end


# --- writing a trace, for tests ---------------------------------------------

def _entry(name: str, number: int, kind: int, size: int, count: int,
           payload: bytes, offset: int) -> bytes:
    """One 28-byte directory entry. Short payloads live in the offset field."""
    inline = len(payload) <= 4
    tail = (payload.ljust(4, b"\x00") if inline
            else struct.pack(">l", offset))
    return (struct.pack(">4slhhll", name.encode("ascii"), number, kind, size,
                        count, len(payload)) + tail + b"\x00\x00\x00\x00")


def write(path: str | Path, sequence: str, quality: list[int] | None = None, *,
          sample: str = "", well: str = "", plate_id: str = "",
          plate_name: str = "", instrument: str = "3730",
          vendor: str = "plateforge-test",
          run_started: datetime | None = None) -> Path:
    """Write a minimal but real ABIF file.

    Not a mock: the bytes are laid out to the spec and read back by `read()`
    through exactly the same code path as a vendor's file. It exists because
    testing every verdict needs 96 traces with known defects, and because
    asserting on a fixture you generated with the parser you are testing
    proves nothing -- so the fixtures in `fixtures/sanger/` stay the ground
    truth for the *format*, and this covers the *logic*.
    """
    path = Path(path)
    sequence = sequence.upper()
    quality = list(quality if quality is not None else [40] * len(sequence))
    started = run_started or datetime(2026, 1, 1, 12, 0, 0)

    def pstring(text: str) -> bytes:
        raw = text.encode("latin-1")[:255]
        return bytes([len(raw)]) + raw

    def cstring(text: str) -> bytes:
        return text.encode("latin-1") + b"\x00"

    items = [
        ("PBAS", 1, CHAR, 1, len(sequence), sequence.encode("latin-1")),
        ("PCON", 1, CHAR, 1, len(quality), bytes(min(max(q, 0), 255) for q in quality)),
        ("SMPL", 1, PSTRING, 1, len(sample) + 1, pstring(sample)),
        ("TUBE", 1, PSTRING, 1, len(well) + 1, pstring(well)),
        ("CTID", 1, CSTRING, 1, len(plate_id) + 1, cstring(plate_id)),
        ("CTNM", 1, CSTRING, 1, len(plate_name) + 1, cstring(plate_name)),
        ("MODL", 1, CHAR, 1, len(instrument), instrument.encode("latin-1")),
        ("User", 1, PSTRING, 1, len(vendor) + 1, pstring(vendor)),
        ("RUND", 1, DATE, 2, 1, struct.pack(">hBB", started.year, started.month,
                                            started.day)),
        ("RUNT", 1, TIME, 1, 1, struct.pack(">4B", started.hour, started.minute,
                                            started.second, 0)),
    ]

    # Header, then the directory, then the payloads. The directory has to know
    # where the payloads land, so sizes are computed before anything is written.
    header = MAGIC + struct.pack(">h", 101)
    directory_at = len(header) + ENTRY_SIZE
    payload_at = directory_at + len(items) * ENTRY_SIZE

    entries_out, payloads_out, cursor = [], [], payload_at
    for name, number, kind, size, count, payload in items:
        entries_out.append(_entry(name, number, kind, size, count, payload, cursor))
        if len(payload) > 4:
            payloads_out.append(payload)
            cursor += len(payload)

    root = _entry("tdir", 1, LONG, ENTRY_SIZE, len(items),
                  b"\x00" * (len(items) * ENTRY_SIZE), directory_at)
    path.write_bytes(header + root + b"".join(entries_out) + b"".join(payloads_out))
    return path
