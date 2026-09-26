"""Deciding which clones bound the target.

A hit is a well whose signal on the target antigen beats its **own** signal
on the control antigen. That framing does most of the work: it cancels
anything that affects a clone equally on both plates -- how much IgG was
prepped, how sticky the clone is, how the row was pipetted.

What it does not cancel is anything that affects the two *plates*
differently, which is why a caller cannot be a single threshold.

## A caller is allowed to say no

The rule that matters most here is not which threshold to use. It is that on
some plates **no threshold is correct**, and the honest output is a refusal.

An over-developed plate is the clearest case. TMB develops toward the reader
ceiling, so once a well saturates its reading stops carrying information
about how much antibody is bound. A strong binder and a moderate one both
read 4.0, the ratio against control collapses toward 1, and a caller applying
its usual threshold quietly returns a short, confident, wrong list. The
answer is "re-read this plate sooner", not a list of hits.

So every caller returns a plate-level `verdict` alongside its per-well calls,
and `uncallable` is a normal outcome rather than an error. Anything that is
not `callable` should reach a human before it reaches a cherry-pick.

## Adding a caller

One function plus a decorator (rule 5). It receives the paired frame and the
plate context, and returns a `CallSet`. Nothing existing moves.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..core import registry

CALLERS = registry.Registry("hit caller")

TARGET, CONTROL = "od450_target", "od450_control"

CALLABLE, DEGRADED, UNCALLABLE = "callable", "degraded", "uncallable"

# Defaults, and they are defaults rather than findings. A ratio of 3 is the
# number people reach for; the evaluation in decision 0024 is what says when
# it holds.
DEFAULTS = {
    "min_ratio": 3.0,
    "min_signal": 0.2,          # absolute OD floor for a hit
    "z_threshold": 3.5,
    "saturation_od": 3.6,       # above this a reading has lost its headroom
    "reader_resolution": 0.001,  # plate readers report OD to three decimals
    # Dynamic range -- top decile of the target plate over its background --
    # is what actually predicts whether calling works. Measured over 555
    # simulated plates that had at least one binder:
    #
    #   range  3-6x   F1 0.28     range 12-30x  F1 0.57
    #   range  6-12x  F1 0.34     range  >30x   F1 0.81
    #
    "uncallable_dynamic_range": 10.0,
    "degraded_dynamic_range": 25.0,
}


@dataclass
class CallSet:
    """Per-well calls, plus what the caller thought of the plate as a whole."""
    caller: str
    calls: pd.DataFrame               # well, clone_id, is_hit, score, flag
    verdict: str = CALLABLE
    reasons: list[str] = field(default_factory=list)
    params: dict = field(default_factory=dict)

    @property
    def hits(self) -> pd.DataFrame:
        return self.calls[self.calls["is_hit"]]

    @property
    def trustworthy(self) -> bool:
        return self.verdict == CALLABLE

    def as_dict(self) -> dict:
        return {"caller": self.caller, "verdict": self.verdict,
                "n_hits": int(self.calls["is_hit"].sum()),
                "n_wells": int(len(self.calls)),
                "reasons": self.reasons, "params": self.params}


# --- what the plate itself says about whether it can be read ---------------

def plate_health(frame: pd.DataFrame, blanks: pd.Series | None = None,
                 **params) -> tuple[str, list[str]]:
    """Whether this pair is in a state where any threshold means something.

    Checked before any calling, because these are properties of the plate
    rather than of a rule, and every caller inherits the same answer.
    """
    p = DEFAULTS | params
    reasons: list[str] = []
    target = frame[TARGET].astype(float)
    verdict = CALLABLE

    background = (float(blanks.median()) if blanks is not None and len(blanks)
                  else float(np.median(np.sort(target)[:max(3, len(target) // 10)])))
    top = float(target.nlargest(max(1, len(target) // 10)).mean())
    dynamic = top / background if background > 0 else float("inf")

    if dynamic < p["uncallable_dynamic_range"]:
        verdict = UNCALLABLE
        reasons.append(
            f"top-decile signal is only {dynamic:.1f}x background. Below "
            f"{p['uncallable_dynamic_range']:.0f}x there is no range to call "
            "hits in, and measured F1 on such plates is about 0.3 whatever "
            "threshold is used -- whether or not anything bound.")
    elif dynamic < p["degraded_dynamic_range"]:
        verdict = DEGRADED
        reasons.append(
            f"top-decile signal is {dynamic:.1f}x background; calls here are "
            "usable but measurably worse (F1 about 0.57 against 0.81 on "
            "plates with more range).")

    # Saturation is FLAGGED, never refused on. Measured across 555 plates,
    # more saturation went with BETTER calls -- F1 0.54 at under 2% saturated
    # against 0.91 above 20% -- because a well only saturates when something
    # really bound. What is lost is the ability to RANK among saturated
    # wells, not the ability to call them.
    saturated = float((target >= p["saturation_od"]).mean())
    if saturated > 0:
        reasons.append(
            f"{saturated:.0%} of target wells are at the reader ceiling: they "
            "are called, but cannot be ranked against each other.")
    return verdict, reasons


def _skeleton(frame: pd.DataFrame, score: pd.Series, is_hit: pd.Series,
              **params) -> pd.DataFrame:
    p = DEFAULTS | params
    calls = pd.DataFrame({
        "well": frame["well"].values,
        "clone_id": frame["clone_id"].values,
        TARGET: frame[TARGET].values,
        CONTROL: frame[CONTROL].values,
        "score": np.asarray(score, dtype=float),
        "is_hit": np.asarray(is_hit, dtype=bool),
    })
    calls["flag"] = np.where(
        calls[TARGET] >= p["saturation_od"], "saturated",
        np.where(calls[TARGET] < p["min_signal"], "below floor", ""))
    return calls


# --- the callers ------------------------------------------------------------

@CALLERS.register("fixed_ratio", needs_blanks=False)
def fixed_ratio(frame: pd.DataFrame, blanks: pd.Series | None = None,
                **params) -> CallSet:
    """target / control >= k, with an absolute floor. The classic.

    The floor matters more than the ratio. Two wells both near background can
    have a ratio of 5 by noise alone, and without a floor a clean plate
    returns a dozen hits that are nothing.
    """
    p = DEFAULTS | params
    ratio = frame[TARGET] / frame[CONTROL].clip(lower=1e-6)
    is_hit = (ratio >= p["min_ratio"]) & (frame[TARGET] >= p["min_signal"])
    verdict, reasons = plate_health(frame, blanks, **params)
    return CallSet("fixed_ratio", _skeleton(frame, ratio, is_hit, **params),
                   verdict, reasons,
                   {"min_ratio": p["min_ratio"], "min_signal": p["min_signal"]})


@CALLERS.register("ratio_over_blank", needs_blanks=True)
def ratio_over_blank(frame: pd.DataFrame, blanks: pd.Series | None = None,
                     **params) -> CallSet:
    """Like `fixed_ratio`, but the floor is measured rather than assumed.

    The empty wells of a partly-filled plate are this plate's background, on
    this day. A fixed 0.2 OD floor is wrong on a plate whose blank is 0.02
    and wrong again on one whose blank is 0.25; a multiple of the measured
    blank is right on both.
    """
    p = DEFAULTS | params
    background = (float(blanks.median()) if blanks is not None and len(blanks)
                  else float(np.median(np.sort(frame[TARGET])[:max(3, len(frame) // 10)])))
    floor = max(background * p.get("blank_multiple", 4.0), 0.05)
    ratio = frame[TARGET] / frame[CONTROL].clip(lower=1e-6)
    is_hit = (ratio >= p["min_ratio"]) & (frame[TARGET] >= floor)
    verdict, reasons = plate_health(frame, blanks, **params)
    return CallSet("ratio_over_blank", _skeleton(frame, ratio, is_hit, **params),
                   verdict, reasons,
                   {"min_ratio": p["min_ratio"], "measured_floor": round(floor, 4),
                    "background": round(background, 4)})


@CALLERS.register("robust_z", needs_blanks=False)
def robust_z(frame: pd.DataFrame, blanks: pd.Series | None = None,
             **params) -> CallSet:
    """Robust z of (target - control) against the plate's own spread.

    Median and MAD rather than mean and SD, because the hits are exactly the
    outliers and letting them set the scale is how a plate with many hits
    calls none of them. Scales itself to the plate, so it survives a shift
    that a fixed ratio does not -- and fails differently, by finding
    "outliers" on a plate where nothing bound.
    """
    p = DEFAULTS | params
    delta = (frame[TARGET] - frame[CONTROL]).astype(float)
    centre = float(delta.median())
    spread = float((delta - centre).abs().median()) * 1.4826
    if spread <= 1e-9:
        # MAD is exactly zero when over half the wells read identically --
        # a plate where almost nothing expressed. Falling back to the plate
        # SD here is the one thing that must not happen: the SD is inflated
        # by the very wells being tested, so a plate of 80 blanks and 16
        # strong binders gets a spread wide enough to call none of them.
        # Use the spread of the middle half instead, which the outliers
        # cannot reach, and floor it at the reader's own resolution.
        middle = delta[delta.between(delta.quantile(0.25), delta.quantile(0.75))]
        spread = max(float(middle.std(ddof=0)) if len(middle) > 1 else 0.0,
                     p["reader_resolution"])
    z = (delta - centre) / spread
    is_hit = (z >= p["z_threshold"]) & (frame[TARGET] >= p["min_signal"])
    verdict, reasons = plate_health(frame, blanks, **params)
    return CallSet("robust_z", _skeleton(frame, z, is_hit, **params),
                   verdict, reasons,
                   {"z_threshold": p["z_threshold"], "median_delta": round(centre, 4),
                    "mad_scaled": round(spread, 4)})


@CALLERS.register("ratio_and_z", needs_blanks=True)
def ratio_and_z(frame: pd.DataFrame, blanks: pd.Series | None = None,
                **params) -> CallSet:
    """Both, and a well has to satisfy each.

    They fail in different directions -- ratio over-calls near background, z
    over-calls on a plate with nothing in it -- so requiring both costs a
    little recall and removes most of each one's characteristic mistake.
    """
    p = DEFAULTS | params
    by_ratio = ratio_over_blank(frame, blanks, **params)
    by_z = robust_z(frame, blanks, **params)
    is_hit = by_ratio.calls["is_hit"].values & by_z.calls["is_hit"].values
    verdict, reasons = plate_health(frame, blanks, **params)
    calls = _skeleton(frame, by_z.calls["score"], is_hit, **params)
    calls["ratio"] = by_ratio.calls["score"].values
    return CallSet("ratio_and_z", calls, verdict, reasons,
                   by_ratio.params | by_z.params)


def call(frame: pd.DataFrame, caller: str = "ratio_and_z",
         blanks: pd.Series | None = None, **params) -> CallSet:
    """Run one caller over a paired frame."""
    return CALLERS.get(caller)(frame, blanks, **params)


def available() -> pd.DataFrame:
    return pd.DataFrame([{"caller": name,
                          "uses_blank_wells": CALLERS.meta(name).get("needs_blanks")}
                         for name in CALLERS])


# --- measuring a caller instead of arguing about it ------------------------

def score_calls(calls: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """Precision, recall and F1 against known binders.

    Only possible on simulated data, which is the point of having it: a
    calling rule should be chosen by measurement, and real plates never come
    with an answer key.
    """
    merged = calls.merge(truth[["well", "is_binder"]], on="well", how="left")
    predicted = merged["is_hit"].fillna(False).astype(bool)
    actual = merged["is_binder"].fillna(False).astype(bool)

    tp = int((predicted & actual).sum())
    fp = int((predicted & ~actual).sum())
    fn = int((~predicted & actual).sum())
    tn = int((~predicted & ~actual).sum())
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if tp and precision == precision and recall == recall else 0.0)
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "precision": precision, "recall": recall, "f1": f1,
            "n_binders": int(actual.sum())}


def evaluate(pairs, callers: list[str] | None = None,
             respect_verdict: bool = True, **params) -> pd.DataFrame:
    """Run every caller over every pair and score it against the truth.

    `respect_verdict` is the interesting switch. With it on, a caller that
    declared a plate uncallable returns no hits for it -- which costs recall
    and is the behaviour you want, because the alternative is a confident
    wrong list. Turning it off measures what the thresholds would have done
    if nobody was allowed to refuse.
    """
    callers = callers or list(CALLERS)
    rows = []
    for pair in pairs:
        truth = pair.truth()
        blanks = None
        if "role" in pair.target.columns:
            empty = pair.target[pair.target["role"] != "sample"]
            blanks = empty["od450"] if len(empty) else None
        for name in callers:
            result = call(truth, name, blanks, **params)
            calls = result.calls
            if respect_verdict and result.verdict == UNCALLABLE:
                calls = calls.assign(is_hit=False)
            rows.append({
                "pair_id": pair.pair_id, "scenario": pair.scenario,
                "caller": name, "verdict": result.verdict,
                "fill": float((pair.target["role"] == "sample").mean())
                if "role" in pair.target.columns else 1.0,
            } | score_calls(calls, truth))
    return pd.DataFrame(rows)


def summarise(scored: pd.DataFrame, by: str = "caller") -> pd.DataFrame:
    """Mean precision/recall/F1, and how often each caller refused."""
    grouped = scored.groupby(by, observed=True)
    out = grouped[["precision", "recall", "f1"]].mean().round(3)
    out["refused"] = grouped["verdict"].apply(
        lambda v: round(float((v == UNCALLABLE).mean()), 3))
    out["plates"] = grouped.size()
    return out.sort_values("f1", ascending=False)
