"""Well coordinates. Canonical form is always zero-padded: A01, not A1.

Every boundary into the system calls normalize(). Nothing else parses wells.
"""
from __future__ import annotations

import re

ROWS = "ABCDEFGHIJKLMNOP"

FORMATS: dict[int, tuple[int, int]] = {
    6: (2, 3),
    12: (3, 4),
    24: (4, 6),
    48: (6, 8),
    96: (8, 12),
    384: (16, 24),
}

_RC = re.compile(r"^\s*([A-Za-z]{1,2})\s*0*(\d{1,2})\s*$")
_CR = re.compile(r"^\s*0*(\d{1,2})\s*([A-Za-z]{1,2})\s*$")


def dims(plate_format: int = 96) -> tuple[int, int]:
    if plate_format not in FORMATS:
        raise ValueError(f"unsupported plate format {plate_format}")
    return FORMATS[plate_format]


def normalize(raw: str, plate_format: int = 96) -> str:
    """'a1', 'A1', ' A 01 ', '1A' -> 'A01'. Raises if outside the format."""
    if raw is None:
        raise ValueError("well is None")
    m = _RC.match(str(raw)) or _CR.match(str(raw))
    if not m:
        raise ValueError(f"unparseable well {raw!r}")
    a, b = m.groups()
    row_s, col_s = (a, b) if a.isalpha() else (b, a)
    return from_rc(row_index(row_s), int(col_s), plate_format)


def row_index(row: str) -> int:
    """'A' -> 1, 'B' -> 2, 'AA' -> 27."""
    row = row.upper()
    n = 0
    for ch in row:
        if ch not in ROWS[:26] and not ch.isalpha():
            raise ValueError(f"bad row {row!r}")
        n = n * 26 + (ord(ch) - 64)
    return n


def row_label(index: int) -> str:
    """1 -> 'A'."""
    if index < 1:
        raise ValueError("row index is 1-based")
    return ROWS[index - 1]


def from_rc(row: int, col: int, plate_format: int = 96) -> str:
    n_rows, n_cols = dims(plate_format)
    if not (1 <= row <= n_rows) or not (1 <= col <= n_cols):
        raise ValueError(f"({row},{col}) outside {plate_format}-well plate")
    return f"{row_label(row)}{col:02d}"


def to_rc(well: str, plate_format: int = 96) -> tuple[int, int]:
    w = normalize(well, plate_format)
    return row_index(w[0]), int(w[1:])


def all_wells(plate_format: int = 96, order: str = "column") -> list[str]:
    """order='column' fills A01,B01,... (pipetting order); 'row' fills A01,A02,..."""
    n_rows, n_cols = dims(plate_format)
    if order == "column":
        return [from_rc(r, c, plate_format) for c in range(1, n_cols + 1) for r in range(1, n_rows + 1)]
    if order == "row":
        return [from_rc(r, c, plate_format) for r in range(1, n_rows + 1) for c in range(1, n_cols + 1)]
    raise ValueError("order must be 'column' or 'row'")


def is_edge(well: str, plate_format: int = 96) -> bool:
    r, c = to_rc(well, plate_format)
    n_rows, n_cols = dims(plate_format)
    return r in (1, n_rows) or c in (1, n_cols)
