"""A run bundle: everything one execution produced, in one directory.

A figure on its own is not a result. Six months later the only question that
matters about a set of 96 sequences is *how these 96 and not others*, and that
answer has to live next to the files rather than in a terminal scrollback that
is long gone.

So every selection writes a bundle:

    $PLATEFORGE_DATA/outputs/<obj_id>/
      manifest.json         what was asked for, what happened, what came out
      selected.csv          the sequences, in rank order
      pool_summary.csv      what they were drawn from
      alignment_<gene>.png  one per V gene
      alignment_<gene>.fasta

The manifest is the point. It records the source and its query parameters, the
filter funnel that took the source down to the pool, the sampler's settings and
seed, the diversity report, and the code version -- so the run can be argued
with, and re-run.

Bundles are append-only in spirit: a re-run mints a new obj_id and gets its own
directory. Nothing here overwrites an earlier run's answer.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import paths, version


def _jsonable(obj: Any) -> Any:
    """Make numpy/pandas scalars and frames survive json.dump."""
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    for attr in ("to_dict", "item", "isoformat"):
        fn = getattr(obj, attr, None)
        if callable(fn):
            try:
                return _jsonable(fn("records") if attr == "to_dict"
                                 and hasattr(obj, "columns") else fn())
            except (TypeError, ValueError):
                continue
    return str(obj)


class Bundle:
    """Collects a run's outputs and writes the manifest when closed."""

    def __init__(self, obj_id: str, kind: str, label: str = ""):
        self.obj_id = obj_id
        self.kind = kind
        self.label = label
        self.dir = paths.output_dir(obj_id)
        self.started = datetime.now(timezone.utc)
        self.sections: dict[str, Any] = {}
        self.files: list[dict] = []

    def section(self, name: str, payload: Any) -> "Bundle":
        """Record one block of the manifest. Later calls merge into earlier ones."""
        existing = self.sections.get(name)
        if isinstance(existing, dict) and isinstance(payload, dict):
            existing.update(_jsonable(payload))
        else:
            self.sections[name] = _jsonable(payload)
        return self

    def note(self, text: str) -> "Bundle":
        self.sections.setdefault("notes", []).append(text)
        return self

    def add(self, path: str | Path, role: str = "") -> Path:
        """Register a file that is already in (or copied into) the bundle."""
        path = Path(path)
        self.files.append({
            "name": path.name,
            "role": role,
            "bytes": path.stat().st_size if path.exists() else None,
        })
        return path

    def write_text(self, name: str, text: str, role: str = "") -> Path:
        path = self.dir / name
        path.write_text(text)
        return self.add(path, role)

    def write_table(self, name: str, df, role: str = "") -> Path:
        path = self.dir / name
        df.to_csv(path, index=False)
        return self.add(path, role)

    def close(self) -> Path:
        manifest = {
            "obj_id": self.obj_id,
            "kind": self.kind,
            "label": self.label,
            "started_utc": self.started.isoformat(),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "code_version": version.code_version(),
            "sections": self.sections,
            "files": sorted(self.files, key=lambda f: f["name"]),
        }
        path = self.dir / "manifest.json"
        path.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")
        return path


def read_manifest(obj_id: str) -> dict:
    return json.loads((paths.output_dir(obj_id) / "manifest.json").read_text())


def describe(obj_id: str) -> str:
    """The run in a paragraph, for when a directory turns up unexplained."""
    m = read_manifest(obj_id)
    s = m.get("sections", {})
    src = s.get("source", {})
    sel = s.get("selection", {})
    lines = [
        f"{m['obj_id']}  ({m['kind']}, {m.get('label', '')})",
        f"  produced {m['finished_utc']} by code {m['code_version']}",
        f"  source: {src.get('kind', '?')} {src.get('ref', '')}"
        f" -> scanned {src.get('scanned', '?')}, kept {src.get('kept', '?')}",
        f"  selection: {sel.get('n_requested', '?')} requested, "
        f"{sel.get('n_selected', '?')} selected, seed {sel.get('seed', '?')}",
        f"  files: {len(m.get('files', []))}",
    ]
    return "\n".join(lines)
