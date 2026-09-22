"""Typed, human-readable object IDs.

One registry of prefixes. Nothing outside this module invents an ID.
An ID is minted once and stored; it is never recomputed downstream.
"""
from __future__ import annotations

import itertools
import re
import uuid
from datetime import datetime, timezone

_PREFIXES: dict[str, str] = {}

_ID_RE = re.compile(r"^(?P<prefix>[A-Z][A-Z0-9]{1,5})-(?P<body>[A-Za-z0-9_.:+-]+)$")


def register_prefix(prefix: str, description: str) -> str:
    """Declare an ID type. Adding a new object type is one call to this."""
    if not re.fullmatch(r"[A-Z][A-Z0-9]{1,5}", prefix):
        raise ValueError(f"bad prefix {prefix!r}: 2-6 chars, uppercase, starts with a letter")
    existing = _PREFIXES.get(prefix)
    if existing is not None and existing != description:
        raise ValueError(f"prefix {prefix!r} already registered as {existing!r}")
    _PREFIXES[prefix] = description
    return prefix


def known_prefixes() -> dict[str, str]:
    return dict(_PREFIXES)


def mint(prefix: str, body: str | None = None) -> str:
    """Create an ID. Body defaults to a short random token."""
    if prefix not in _PREFIXES:
        raise KeyError(f"unregistered prefix {prefix!r}; call register_prefix first")
    body = body if body is not None else uuid.uuid4().hex[:10]
    obj_id = f"{prefix}-{body}"
    if not _ID_RE.match(obj_id):
        raise ValueError(f"illegal id body {body!r}")
    return obj_id


_SEQ = itertools.count()


def mint_stamped(prefix: str, label: str) -> str:
    """Mint an ID whose body carries a label and a UTC timestamp.

    The timestamp is only second-resolution, so two mints in the same second
    would collide on their own. The counter makes that impossible within a
    process; the random token covers concurrent processes.
    """
    return mint(prefix, f"{label}-{stamp()}-{uuid.uuid4().hex[:6]}{next(_SEQ):x}")


def parse(obj_id: str) -> tuple[str, str]:
    """Split an ID. Raises on unknown prefixes so typos fail loudly."""
    m = _ID_RE.match(obj_id)
    if not m:
        raise ValueError(f"malformed id {obj_id!r}")
    prefix = m.group("prefix")
    if prefix not in _PREFIXES:
        raise KeyError(f"unknown id prefix {prefix!r} in {obj_id!r}")
    return prefix, m.group("body")


def prefix_of(obj_id: str) -> str:
    return parse(obj_id)[0]


def is_a(obj_id: str, prefix: str) -> bool:
    try:
        return prefix_of(obj_id) == prefix
    except (ValueError, KeyError):
        return False


def stamp() -> str:
    """UTC timestamp usable inside an ID body: 20260912T174501Z"""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# --- core object types -------------------------------------------------
SEQ = register_prefix("SEQ", "single antibody sequence record")
LIB = register_prefix("LIB", "a named, versioned set of sequences")
FMT = register_prefix("FMT", "a reformatted construct derived from a sequence")
CLN = register_prefix("CLN", "a clone tracked through the wet-lab workflow")
PLT = register_prefix("PLT", "a plate, real or virtual")
EXP = register_prefix("EXP", "an experiment definition")
RUN = register_prefix("RUN", "one execution of an analysis")
PCK = register_prefix("PCK", "a pick session selecting hits")
XFR = register_prefix("XFR", "a transfer list for a liquid handler")
RGT = register_prefix("RGT", "a reagent catalog entry")
LOT = register_prefix("LOT", "a physical lot of a reagent")
CAV = register_prefix("CAV", "a caveat raised by the reagent rules layer")
DOC = register_prefix("DOC", "a generated file or report")
