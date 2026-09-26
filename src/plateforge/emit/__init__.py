"""Files for external systems: instrument worklists, order forms, reports.

Every emitter declares how well it actually knows its format (rule 6). The
default is our own, which cannot be wrong about anyone else's -- except for
`sequencing`, where the vendors' real templates are in `fixtures/sequencing/`
and the emitters are checked against them.
"""
from __future__ import annotations

from . import sequencing, worklists          # noqa: F401  (registers emitters)

__all__ = ["sequencing", "worklists"]
