"""Capture which code produced a result.

A six-month-old artifact is only reproducible if you know what was running.
"""
from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]


@lru_cache(maxsize=1)
def code_version() -> str:
    """Short git SHA, with '+dirty' if the tree has uncommitted changes."""
    try:
        sha = subprocess.run(
            ["git", "-C", str(_REPO), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        if sha.returncode != 0:
            return "unknown"
        rev = sha.stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(_REPO), "status", "--porcelain"],
            capture_output=True, text=True, timeout=5,
        )
        return rev + ("+dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return "unknown"
