"""Files for external systems: instrument worklists, order forms, reports.

Every emitter declares how well it actually knows its format (rule 6). The
default is our own, which cannot be wrong about anyone else's.
"""
from __future__ import annotations

from . import worklists                       # noqa: F401  (registers emitters)

__all__ = ["worklists"]
