"""Reading what came off the plate reader.

A BioTek/Agilent Gen5 export is a workbook with one sheet per plate read --
everything scanned in that session, in order. Somewhere in each sheet is a
grid: a header row of column numbers, then rows labelled A..H.

## Why this does not parse the vendor's layout

Rule 6 says no emitter without a real example file, and readers are held to
the same standard. There is no real Gen5 export in `fixtures/` yet, and Gen5
exports vary a great deal -- matrix or column layout, with or without a
metadata preamble, one read per sheet or several stacked.

So this does not try to know Gen5's layout. It **finds the grid by its
shape**: a run of consecutive integers 1..12 (or 1..24) on one row, with
rows beneath it labelled A..H (or A..P), enclosing a rectangle of numbers.
That description is true of every microplate export anyone produces, which
makes it more robust than a parser written against one vendor's template and
then broken by a software update.

What it cannot do without a real file is read the **metadata** reliably --
which read is 450 nm, which sheet is which plate, what the timestamps mean.
It collects the text above each grid verbatim and hands it over unparsed, so
a human or the assignment step can use it without anyone pretending it was
understood.

## Then the grids have to become plates

A sheet is not a plate. `assay.assign` matches grids to the virtual plates
the pipeline already knows about -- by barcode, by name, or by asking the
classifier which one looks like the target.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..core import registry, wells as wellmod

READERS = registry.Registry("plate reader export")

# How well a reader's format is known. Same vocabulary as `emit.worklists`,
# plus one level those do not need.
VERIFIED = "verified"        # a real export from this instrument is in fixtures/
HEURISTIC = "heuristic"      # does not depend on the vendor layout at all
SKETCH = "sketch"            # neither

ROW_LABEL = re.compile(r"^\s*([A-P])\s*$", re.IGNORECASE)


@dataclass
class Grid:
    """One plate-shaped block of numbers found in a sheet."""
    sheet: str
    values: dict[str, float]          # well -> value, wells normalised
    plate_format: int
    header_row: int
    first_row: int
    preamble: list[str] = field(default_factory=list)
    label: str = ""

    @property
    def n_values(self) -> int:
        return sum(1 for v in self.values.values() if v == v)   # NaN-safe

    def frame(self, channel: str = "od450") -> pd.DataFrame:
        order = wellmod.all_wells(self.plate_format, "column")
        return pd.DataFrame({"well": order,
                             channel: [self.values.get(w) for w in order]})

    def describe(self) -> str:
        return (f"{self.sheet}" + (f" / {self.label}" if self.label else "")
                + f": {self.plate_format}-well, {self.n_values} values, "
                  f"row {self.header_row + 1}")


def _as_number(cell) -> float | None:
    if cell is None or (isinstance(cell, float) and cell != cell):
        return None
    if isinstance(cell, (int, float)):
        return float(cell)
    text = str(cell).strip()
    if not text:
        return None
    # Readers write saturation as text surprisingly often.
    if text.upper() in {"OVRFLW", "OVER", "SATURATED", "*****"}:
        return float("inf")
    try:
        return float(text)
    except ValueError:
        return None


def _header_columns(row: list) -> tuple[int, int] | None:
    """(first index, count) of a run of consecutive integers starting at 1."""
    best = None
    for start in range(len(row)):
        value = _as_number(row[start])
        if value != 1:
            continue
        count, i = 1, start + 1
        while i < len(row):
            nxt = _as_number(row[i])
            if nxt is None or nxt != count + 1:
                break
            count += 1
            i += 1
        if count >= 6 and (best is None or count > best[1]):
            best = (start, count)
    return best


def find_grids(rows: list[list], sheet: str = "") -> list[Grid]:
    """Every plate-shaped block in one sheet, by shape rather than by layout."""
    grids: list[Grid] = []
    used_rows: set[int] = set()

    for r, row in enumerate(rows):
        if r in used_rows:
            continue
        header = _header_columns(row)
        if header is None:
            continue
        first_col, n_cols = header

        # Row labels start in a column at or before the numbers begin.
        found: dict[str, float] = {}
        labels: list[str] = []
        body = r + 1
        while body < len(rows):
            candidate = rows[body]
            label = None
            for c in range(0, first_col + 1):
                if c < len(candidate):
                    match = ROW_LABEL.match(str(candidate[c] or ""))
                    if match:
                        label = match.group(1).upper()
                        break
            if label is None or label in labels:
                break
            labels.append(label)
            for j in range(n_cols):
                c = first_col + j
                value = _as_number(candidate[c]) if c < len(candidate) else None
                if value is not None:
                    found[f"{label}{j + 1:02d}"] = value
            body += 1

        if len(labels) < 4 or not found:
            continue

        n_rows = len(labels)
        plate_format = n_rows * n_cols
        if plate_format not in (6, 12, 24, 48, 96, 384, 1536):
            plate_format = 96 if plate_format <= 96 else 384

        # Empty cells arrive as NaN from pandas, and "nan nan nan" is not
        # metadata. Anything that survives here is text a human typed or the
        # instrument wrote.
        preamble = []
        for back in range(max(0, r - 6), r):
            parts = []
            for cell in rows[back]:
                if cell is None or (isinstance(cell, float) and cell != cell):
                    continue
                text = str(cell).strip()
                if text and text.lower() != "nan":
                    parts.append(text)
            if parts:
                preamble.append(" ".join(parts))

        grids.append(Grid(sheet=sheet, values=found, plate_format=plate_format,
                          header_row=r, first_row=r + 1, preamble=preamble,
                          label=preamble[-1][:80] if preamble else ""))
        used_rows.update(range(r, body))
    return grids


# --- the registered readers -------------------------------------------------

@READERS.register("biotek_gen5", software="BioTek / Agilent Gen5",
                  confidence=HEURISTIC, suffix=".xlsx")
def biotek_gen5(source: str | Path, **_) -> list[Grid]:
    """Every plate grid in every sheet of a Gen5 workbook, in sheet order.

    Sheet order is preserved because it is usually read order, which is the
    single most useful clue when matching grids to plates.
    """
    sheets = pd.read_excel(source, sheet_name=None, header=None, dtype=object)
    grids: list[Grid] = []
    for name, frame in sheets.items():
        grids.extend(find_grids(frame.values.tolist(), sheet=str(name)))
    return grids


@READERS.register("grid_csv", software="any CSV holding one plate grid",
                  confidence=HEURISTIC, suffix=".csv")
def grid_csv(source: str | Path, **_) -> list[Grid]:
    """One plate laid out as a grid in a CSV. The common hand-made case."""
    frame = pd.read_csv(source, header=None, dtype=object)
    return find_grids(frame.values.tolist(), sheet=Path(source).stem)


@READERS.register("long_csv", software="ours", confidence=VERIFIED,
                  suffix=".csv")
def long_csv(source: str | Path, well_column: str = "well",
             value_column: str = "od450", plate_column: str | None = None,
             **_) -> list[Grid]:
    """One row per well. Our own shape, so nothing here is a guess."""
    frame = pd.read_csv(source)
    groups = ([(str(k), g) for k, g in frame.groupby(plate_column)]
              if plate_column and plate_column in frame.columns
              else [(Path(source).stem, frame)])
    out = []
    for name, group in groups:
        values = {wellmod.normalize(str(r[well_column])): float(r[value_column])
                  for _, r in group.iterrows()}
        out.append(Grid(sheet=name, values=values,
                        plate_format=96 if len(values) <= 96 else 384,
                        header_row=-1, first_row=-1, label=name))
    return out


def read(source: str | Path, reader: str = "biotek_gen5",
         allow_unverified: bool = True, **kwargs) -> list[Grid]:
    """Parse an export into grids.

    Unlike `emit`, reading defaults to permitted: misreading a file is
    recoverable and the result is visible on screen, whereas writing a bad
    worklist moves liquid. The confidence is still reported.
    """
    meta = READERS.meta(reader)
    if meta.get("confidence") == SKETCH and not allow_unverified:
        raise ValueError(f"{reader!r} has no usable parser yet")
    return READERS.get(reader)(source, **kwargs)


def available() -> pd.DataFrame:
    return pd.DataFrame([
        {"reader": name,
         "software": READERS.meta(name).get("software"),
         "confidence": READERS.meta(name).get("confidence"),
         "suffix": READERS.meta(name).get("suffix")}
        for name in READERS
    ])


# --- something to test against, until a real export arrives -----------------

def write_gen5_like(path: str | Path, plates: list[tuple[str, pd.DataFrame]],
                    *, channel: str = "od450", wavelength: str = "450",
                    plate_format: int = 96, preamble: bool = True) -> Path:
    """Write a workbook shaped like a Gen5 export. SYNTHETIC.

    This is not a claim about Gen5's format -- it is a fixture so the grid
    finder can be tested, and so there is something concrete to compare a
    real export against when one arrives. It deliberately includes the
    awkward parts: a metadata preamble of varying height, the grid not
    starting at A1, and blank rows between blocks.
    """
    import openpyxl

    book = openpyxl.Workbook()
    book.remove(book.active)
    n_rows, n_cols = wellmod.dims(plate_format)
    # Excel refuses []:*?/\ in a sheet title, and a real plate name may well
    # contain one. The original goes into the preamble either way, which is
    # where the grid finder reads it from.
    forbidden = str.maketrans({c: "_" for c in "[]:*?/\\"})

    for i, (name, frame) in enumerate(plates):
        sheet = book.create_sheet(title=str(name).translate(forbidden)[:31])
        row = 1
        if preamble:
            for line in [("Software Version", "3.11.19"),
                         ("Experiment File Path:", f"C:\\Gen5\\{name}.xpt"),
                         ("Protocol File Path:", "C:\\Gen5\\IgG_ELISA.prt"),
                         ("Plate Number", name),
                         ("Date", "9/24/2026"),
                         ("Reader Type:", "Synergy H1"),
                         (f"Read {i + 1}:{wavelength}", None)]:
                sheet.cell(row=row, column=1, value=line[0])
                if line[1] is not None:
                    sheet.cell(row=row, column=2, value=line[1])
                row += 1
            row += 1                       # a blank line, as exports have

        offset = 2 if preamble else 1      # grids rarely start in column A
        for j in range(n_cols):
            sheet.cell(row=row, column=offset + 1 + j, value=j + 1)
        lookup = dict(zip(frame["well"], frame[channel]))
        for r in range(n_rows):
            label = wellmod.row_label(r + 1)
            sheet.cell(row=row + 1 + r, column=offset, value=label)
            for j in range(n_cols):
                value = lookup.get(f"{label}{j + 1:02d}")
                if value is not None:
                    sheet.cell(row=row + 1 + r, column=offset + 1 + j,
                               value=float(value))

    path = Path(path)
    book.save(path)
    return path
