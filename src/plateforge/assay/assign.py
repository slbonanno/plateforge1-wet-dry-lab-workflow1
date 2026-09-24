"""Matching what came off the reader to the plates the pipeline knows about.

A sheet is not a plate. The reader hands back grids in scan order with
whatever text happened to sit above them; the pipeline holds `PLT` artifacts
with barcodes, clone maps and lineage. Joining the two is its own step, and
it is the step where a quiet mistake is most expensive -- swap target and
control and every hit call inverts.

So the ordering of evidence is deliberate:

    1. the user said so                 explicit, and always wins
    2. a barcode appears in the sheet   near-certain
    3. a plate name matches             good
    4. the classifier recognises it     a suggestion, never a decision

Anything below "the user said so" is returned as a *proposal* with a
confidence and the reason for it. `Assignment.needs_confirmation` is true
whenever the pipeline guessed, and the caller is expected to ask. An agent
holding this should show its reasoning and wait, not assign silently.

The classifier only ever proposes which of a **pair** is the target. Asked
about one plate in isolation it will answer, and on a panel that yielded
nothing that answer is a coin flip -- see decision 0022. The honest use is
"these two plates, which is which", where the comparison carries most of the
information.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from ..core import ids, paths
from . import elisa

EXPLICIT, BARCODE, NAME, MODEL, UNRESOLVED = (
    "explicit", "barcode", "name", "model", "unresolved")

# Ordered best first. A later source never overrides an earlier one.
PRECEDENCE = [EXPLICIT, BARCODE, NAME, MODEL, UNRESOLVED]


@dataclass
class Assignment:
    """One grid, and what we think it is."""
    grid_index: int
    sheet: str
    plate_id: str | None = None
    plate_kind: str | None = None          # target | control | unknown
    source: str = UNRESOLVED
    confidence: float = 0.0
    reason: str = ""
    alternatives: list[str] = field(default_factory=list)

    @property
    def needs_confirmation(self) -> bool:
        return self.source in (MODEL, UNRESOLVED)

    def as_dict(self) -> dict:
        return {"grid_index": self.grid_index, "sheet": self.sheet,
                "plate_id": self.plate_id, "plate_kind": self.plate_kind,
                "source": self.source, "confidence": round(self.confidence, 4),
                "reason": self.reason, "alternatives": self.alternatives}


def _text_of(grid) -> str:
    return " ".join([grid.sheet, grid.label] + list(grid.preamble)).lower()


def by_barcode(grids, barcodes: dict[str, str]) -> dict[int, Assignment]:
    """Find a known barcode written anywhere in the sheet or its preamble."""
    out = {}
    for i, grid in enumerate(grids):
        haystack = _text_of(grid)
        hits = [(bc, pid) for bc, pid in barcodes.items()
                if bc and bc.lower() in haystack]
        if len(hits) == 1:
            bc, pid = hits[0]
            out[i] = Assignment(i, grid.sheet, plate_id=pid, source=BARCODE,
                                confidence=0.99,
                                reason=f"barcode {bc!r} appears in the sheet")
        elif len(hits) > 1:
            out[i] = Assignment(i, grid.sheet, source=UNRESOLVED,
                                confidence=0.0,
                                reason=f"{len(hits)} barcodes appear in one sheet",
                                alternatives=[pid for _, pid in hits])
    return out


def by_name(grids, names: dict[str, str]) -> dict[int, Assignment]:
    """Match a plate name or label appearing in the sheet."""
    out = {}
    for i, grid in enumerate(grids):
        haystack = _text_of(grid)
        hits = [(n, pid) for n, pid in names.items()
                if n and len(n) > 3 and n.lower() in haystack]
        if len(hits) == 1:
            n, pid = hits[0]
            out[i] = Assignment(i, grid.sheet, plate_id=pid, source=NAME,
                                confidence=0.8,
                                reason=f"name {n!r} appears in the sheet")
    return out


def kind_from_text(grids) -> dict[int, Assignment]:
    """The sheet often just says. Cheap, and worth trying before a model."""
    target_words = ("target", "antigen", "experimental", "test", "expt")
    control_words = ("control", "ctrl", "hla", "irrelevant", "negative", "neg")
    out = {}
    for i, grid in enumerate(grids):
        haystack = _text_of(grid)
        is_control = any(w in haystack for w in control_words)
        is_target = any(w in haystack for w in target_words)
        if is_control and not is_target:
            out[i] = Assignment(i, grid.sheet, plate_kind="control",
                                source=NAME, confidence=0.75,
                                reason="the sheet says control")
        elif is_target and not is_control:
            out[i] = Assignment(i, grid.sheet, plate_kind="target",
                                source=NAME, confidence=0.75,
                                reason="the sheet says target")
    return out


# --- the model --------------------------------------------------------------

FEATURE_ORDER_KEY = "features"


@dataclass
class PlateModel:
    """A fitted target-vs-control classifier, and what it was fitted on."""
    model: object
    features: list[str]
    accuracy: float
    trained_on: dict
    model_id: str = ""

    def predict(self, feature_rows: pd.DataFrame) -> pd.Series:
        X = feature_rows.reindex(columns=self.features).fillna(-1.0)
        return pd.Series(self.model.predict_proba(X)[:, 1],
                         index=feature_rows.index, name="p_target")


def train(pairs: list, *, random_state: int = 0) -> PlateModel:
    """Fit the target-vs-control classifier on simulated pairs.

    Grouped by pair so a plate never appears in the fold that scored it --
    the two plates of one experiment share their clones and their scenario,
    and splitting them across folds inflates accuracy by a lot.
    """
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import GroupKFold, cross_val_predict

    table = elisa.feature_table(pairs)
    features = [c for c in table.columns
                if c not in ("pair_id", "scenario", "plate_kind", "is_target")]
    X = table[features].fillna(-1.0)
    y = table["is_target"].astype(int)

    predicted = cross_val_predict(
        GradientBoostingClassifier(random_state=random_state), X, y,
        groups=table["pair_id"], cv=GroupKFold(5))
    accuracy = float((predicted == y).mean())

    model = GradientBoostingClassifier(random_state=random_state).fit(X, y)
    return PlateModel(model=model, features=features, accuracy=accuracy,
                      trained_on={
                          "SIMULATED": True,
                          "pairs": len(pairs),
                          "plates": int(len(table)),
                          "scenarios": table.drop_duplicates("pair_id")["scenario"]
                                            .value_counts().to_dict(),
                          "cv": "GroupKFold(5), grouped by pair",
                          "fitted_utc": datetime.now(timezone.utc).isoformat(
                              timespec="seconds"),
                      })


def save(model: PlateModel, label: str = "plate-kind") -> Path:
    """Persist a fitted model as an artifact directory."""
    import joblib

    model.model_id = model.model_id or ids.mint_stamped("DOC", label)
    directory = paths.output_dir(model.model_id)
    joblib.dump(model.model, directory / "model.joblib")
    (directory / "model.json").write_text(json.dumps({
        "model_id": model.model_id,
        "kind": "plate_target_vs_control",
        # Not rounded: a recorded metric that changes when you reload it is
        # a small dishonesty, and rounding belongs at the point of display.
        "accuracy": model.accuracy,
        FEATURE_ORDER_KEY: model.features,
        "trained_on": model.trained_on,
        "warning": "fitted on SIMULATED plates; re-fit on real data before "
                   "trusting it on real data",
    }, indent=2) + "\n")
    return directory


def load(model_id: str) -> PlateModel:
    import joblib

    directory = paths.output_dir(model_id)
    meta = json.loads((directory / "model.json").read_text())
    return PlateModel(model=joblib.load(directory / "model.joblib"),
                      features=meta[FEATURE_ORDER_KEY],
                      accuracy=meta["accuracy"],
                      trained_on=meta["trained_on"], model_id=model_id)


def latest(kind: str = "plate_target_vs_control") -> str | None:
    """The most recently saved model of this kind, if any."""
    root = paths.data_root() / "outputs"
    best = None
    for candidate in sorted(root.glob("DOC-*")):
        meta = candidate / "model.json"
        if not meta.exists():
            continue
        try:
            if json.loads(meta.read_text()).get("kind") == kind:
                best = candidate.name
        except (OSError, json.JSONDecodeError):
            continue
    return best


# --- putting it together ----------------------------------------------------

def propose(grids, *, barcodes: dict | None = None, names: dict | None = None,
            explicit: dict[int, str] | None = None,
            model: PlateModel | None = None,
            channel: str = elisa.CHANNEL,
            roles: dict[int, pd.Series] | None = None,
            n_samples: int | dict[int, int] | None = None) -> list[Assignment]:
    """Best available guess for every grid, with its reason and confidence.

    Evidence is applied best-first and never overwritten by something weaker,
    so an explicit choice survives a confident-looking model.
    """
    found: dict[int, Assignment] = {}

    def offer(candidates: dict[int, Assignment]) -> None:
        for i, assignment in candidates.items():
            existing = found.get(i)
            if existing is None:
                found[i] = assignment
            elif (PRECEDENCE.index(assignment.source)
                  < PRECEDENCE.index(existing.source)):
                found[i] = assignment
            elif existing.plate_kind is None and assignment.plate_kind:
                existing.plate_kind = assignment.plate_kind

    if explicit:
        offer({i: Assignment(i, grids[i].sheet, plate_id=value,
                             plate_kind=value if value in ("target", "control")
                             else None,
                             source=EXPLICIT, confidence=1.0,
                             reason="given by the user")
               for i, value in explicit.items()})
    if barcodes:
        offer(by_barcode(grids, barcodes))
    if names:
        offer(by_name(grids, names))
    offer(kind_from_text(grids))

    undecided = [i for i in range(len(grids))
                 if found.get(i) is None or found[i].plate_kind is None]
    if model is not None and undecided:
        rows, notes = [], {}
        for i in undecided:
            frame = grids[i].frame(channel).dropna(subset=[channel])
            role = (roles or {}).get(i)
            note = ""
            if role is None:
                # No plate map yet -- this step is what produces one. The
                # count of samples is the next best thing and whoever ran the
                # plate knows it; without either, the features are computed
                # over every well and that is said out loud, because on a
                # half-filled plate they are then wrong rather than noisy.
                count = (n_samples.get(i) if isinstance(n_samples, dict)
                         else n_samples)
                role = elisa.infer_roles(frame[channel], count)
                note = ("; occupancy from n_samples" if count
                        else "; PLATE MAP UNKNOWN, features over all wells")
            notes[i] = note
            rows.append(elisa.plate_features(frame[channel], role))
        scores = model.predict(pd.DataFrame(rows, index=undecided))
        for i, p in scores.items():
            kind = "target" if p >= 0.5 else "control"
            confidence = float(max(p, 1 - p))
            reason = (f"model: p(target)={p:.2f}; trained on simulated plates, "
                      f"held-out accuracy {model.accuracy:.0%}"
                      + notes.get(i, ""))
            if i in found and found[i].plate_id:
                found[i].plate_kind = kind
                found[i].reason += f"; {reason}"
                found[i].confidence = min(found[i].confidence, confidence)
            else:
                found[i] = Assignment(i, grids[i].sheet, plate_kind=kind,
                                      source=MODEL, confidence=confidence,
                                      reason=reason)

    return [found.get(i) or Assignment(i, grids[i].sheet, source=UNRESOLVED,
                                       reason="nothing matched")
            for i in range(len(grids))]


def pair_up(assignments: list[Assignment], model: PlateModel | None = None,
            grids=None, channel: str = elisa.CHANNEL) -> list[tuple[int, int]]:
    """Pair consecutive grids as target/control, in scan order.

    Reads are usually run target-then-control, or all targets then all
    controls; scan order is the strongest structural clue available and this
    only commits to the pairing, not to which member is which.
    """
    pairs = []
    for i in range(0, len(assignments) - 1, 2):
        pairs.append((i, i + 1))
    return pairs


def report(assignments: list[Assignment]) -> pd.DataFrame:
    """What was decided, how, and what still needs a human."""
    return pd.DataFrame([a.as_dict() for a in assignments])


def unresolved(assignments: list[Assignment]) -> list[Assignment]:
    return [a for a in assignments if a.needs_confirmation]
