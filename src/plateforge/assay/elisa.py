"""Simulated ELISA plates. SYNTHETIC DATA — never present this as real.

Same rule as `library.synth`: this exists so the analysis path can be built
and tested before a real plate reader export arrives, and every artifact it
produces is stamped `simulated`. A synthetic plate presented as real is the
most damaging mistake this repo can make.

## What is being modelled

A pair of plates read at 450 nm, one coated with the **target antigen** and
one with an **irrelevant control antigen** (HLA-His), the same 96 clones in
the same wells on both. That pairing is the whole point: a clone is a hit
because its signal on target exceeds its signal on control, not because it is
high in absolute terms.

The chain of causes, in order:

    clone  ->  is it expressed at all?          (failed wells)
           ->  does it bind the target?         (specific affinity)
           ->  does it stick to anything?       (non-specific background)
    plate  ->  how well did the antigen coat?
           ->  how long was TMB left in?
           ->  how well was it washed?
           ->  evaporation at the edges, reader drift, pipetting scatter

OD is not linear in bound antibody. TMB develops toward a ceiling, so a plate
left too long compresses everything near the top of the reader's range and
destroys the very contrast the assay is for. That is modelled explicitly:

    OD = blank + (ceiling - blank) * (1 - exp(-k * bound * minutes))

which is why an over-developed plate looks saturated rather than simply
brighter, and why the control plate's *background* rises at the same time.

## Why this is hard to classify, which is the point

Telling an experimental plate from a control plate by its signal pattern is
easy when both are well-behaved: target is right-skewed with a clear hit
population, control is a tight low blob. It stops being easy when TMB runs
long — the experimental plate saturates into a tight HIGH blob, and a control
plate with poor washing becomes broad and high. The confusable cases are
generated deliberately, in known proportions, because a classifier trained
only on clean plates will fail on exactly the plates someone needs help with.

Every number below is a placeholder with a plausible magnitude, declared in
`ASSUMED`. None of them has been fitted to a real plate.
"""
from __future__ import annotations

import math
import random
from dataclasses import asdict, dataclass, field

import pandas as pd

from ..core import wells as wellmod

CHANNEL = "od450"

# Placeholder magnitudes. Plausible, not measured -- the same discipline as
# `assay.protocol.assumed`. Fit these the first time a real export lands.
ASSUMED = {
    "reader_ceiling_od": 4.0,          # most readers clip around here
    "blank_od": 0.05,                  # substrate on a clean, well-washed plate
    "nominal_tmb_minutes": 15.0,
    "hit_rate": 0.18,                  # fraction of clones that bind the target
    "expression_failure_rate": 0.05,   # wells with no usable IgG
    "pipetting_cv": 0.07,
    "edge_effect": 0.12,               # evaporation lifts the outer ring
    "well_failure_rate": 0.02,         # a well that got no detection antibody
}


# --- what a clone is --------------------------------------------------------

@dataclass
class CloneBehaviour:
    """The latent truth about one clone. This is what a classifier never sees."""
    clone_id: str
    expressed: bool
    affinity: float            # specific binding to the TARGET antigen
    control_affinity: float    # specific binding to the CONTROL antigen
    stickiness: float          # non-specific binding, hits both plates

    @property
    def is_binder(self) -> bool:
        return self.expressed and self.affinity > 0.0

    @property
    def is_control_binder(self) -> bool:
        return self.expressed and self.control_affinity > 0.0


# How many of a plate's clones actually bind the target. A single fixed rate
# is the thing that makes a simulated dataset too easy: it guarantees every
# target plate has a hit population and every control plate does not, so a
# classifier learns "does this plate have a right tail" and learns nothing
# about the assay. Real panels vary enormously, and a panel with NO hits is
# common -- it is most of what screening a naive repertoire feels like.
HIT_RATE_MIX = [
    (0.00, 0.00, 0.14),        # nothing bound. Happens, and often.
    (0.01, 0.04, 0.26),        # a handful
    (0.05, 0.15, 0.34),        # a normal campaign
    (0.15, 0.35, 0.20),        # a good day
    (0.35, 0.60, 0.06),        # an enriched or pre-panned library
]

# Some clones bind the control antigen for real. Anti-tag binders are the
# usual culprit -- both antigens carry a His tag -- and a few clones simply
# prefer HLA. Without them a control plate can never have a right tail, and
# the classification task collapses to "which plate has outliers".
CONTROL_BINDER_RATE = 0.02
ANTI_TAG_RATE = 0.012          # binds the tag, so it lights up BOTH plates


def draw_hit_rate(rng: random.Random) -> float:
    """How productive this particular panel was."""
    lo, hi, _ = rng.choices(HIT_RATE_MIX, weights=[w for *_, w in HIT_RATE_MIX])[0]
    return rng.uniform(lo, hi)


def draw_clones(clone_ids: list[str], rng: random.Random, *,
                hit_rate: float | None = None,
                expression_failure_rate: float = ASSUMED["expression_failure_rate"],
                control_binder_rate: float = CONTROL_BINDER_RATE,
                anti_tag_rate: float = ANTI_TAG_RATE,
                ) -> dict[str, CloneBehaviour]:
    """Decide, once per clone, what it actually is.

    A clone's behaviour belongs to the clone, not the plate: the same clone
    read on two plates on two days is the same molecule. Drawing it here is
    what makes a pair of plates *paired*.
    """
    if hit_rate is None:
        hit_rate = draw_hit_rate(rng)
    out = {}
    for clone_id in clone_ids:
        expressed = rng.random() > expression_failure_rate
        binder = expressed and rng.random() < hit_rate
        control_binder = expressed and rng.random() < control_binder_rate
        anti_tag = expressed and rng.random() < anti_tag_rate

        # Affinities span orders of magnitude; a few are very good.
        affinity = math.exp(rng.gauss(0.2, 1.1)) if binder else 0.0
        control_affinity = math.exp(rng.gauss(-0.3, 1.0)) if control_binder else 0.0
        if anti_tag:
            # The tag is on both antigens, so this clone is a hit on both
            # plates and a hit on neither question anyone is asking.
            shared = math.exp(rng.gauss(0.0, 0.9))
            affinity += shared
            control_affinity += shared
        # Almost everything is a little sticky; a few are a lot.
        stickiness = math.exp(rng.gauss(-3.2, 0.95))
        out[clone_id] = CloneBehaviour(clone_id, expressed, affinity,
                                       control_affinity, stickiness)
    return out


# --- what a plate is --------------------------------------------------------

@dataclass
class PlateConditions:
    """How this particular plate was run. Applies to the whole plate."""
    well_failure_rate: float = ASSUMED["well_failure_rate"]
    tmb_minutes: float = ASSUMED["nominal_tmb_minutes"]
    coating: float = 1.0               # antigen density, relative to nominal
    wash_quality: float = 1.0          # 1.0 clean; lower leaves background
    blank_od: float = ASSUMED["blank_od"]
    ceiling_od: float = ASSUMED["reader_ceiling_od"]
    edge_effect: float = ASSUMED["edge_effect"]
    row_drift: float = 0.0             # reader or dispense gradient down rows
    col_drift: float = 0.0
    cv: float = ASSUMED["pipetting_cv"]
    label: str = "nominal"

    def as_dict(self) -> dict:
        return asdict(self)


# Named plate behaviours. Each is a thing that actually happens, and the
# reason it is here is that it changes the SHAPE of the plate's distribution,
# which is what a classifier has to survive.
SCENARIOS: dict[str, dict] = {
    "nominal": {
        "_comment": "everything went right",
        "target": {}, "control": {},
    },
    "tmb_overdeveloped": {
        "_comment": "TMB left in far too long: target saturates into a tight "
                    "high blob, control background climbs. The confusable case.",
        "target": {"tmb_minutes": 45.0},
        "control": {"tmb_minutes": 45.0, "blank_od": 0.18},
    },
    "tmb_underdeveloped": {
        "_comment": "read too early; real hits are present but compressed low",
        "target": {"tmb_minutes": 5.0}, "control": {"tmb_minutes": 5.0},
    },
    "poor_wash_control": {
        "_comment": "control plate washed badly: broad, high, and easily "
                    "mistaken for a target plate",
        "target": {},
        "control": {"wash_quality": 0.45, "blank_od": 0.22},
    },
    "poor_wash_both": {
        "_comment": "whole run washed badly",
        "target": {"wash_quality": 0.5, "blank_od": 0.2},
        "control": {"wash_quality": 0.45, "blank_od": 0.22},
    },
    "weak_coating": {
        "_comment": "target antigen coated poorly; hits are muted and the two "
                    "plates look alike",
        "target": {"coating": 0.25}, "control": {},
    },
    "edge_evaporation": {
        "_comment": "plate sat uncovered; the outer ring reads high on both",
        "target": {"edge_effect": 0.45}, "control": {"edge_effect": 0.45},
    },
    "reader_gradient": {
        "_comment": "dispense or read gradient across the plate",
        "target": {"row_drift": 0.05, "col_drift": 0.02},
        "control": {"row_drift": 0.05, "col_drift": 0.02},
    },
    "sloppy_pipetting": {
        "_comment": "high well-to-well scatter",
        "target": {"cv": 0.28}, "control": {"cv": 0.28},
    },
}

# How often each scenario turns up. Deliberately not uniform: most plates are
# fine, and the confusable ones are a minority -- which is also what makes a
# naive classifier look good on paper and fail in the hands.
SCENARIO_WEIGHTS = {
    "nominal": 0.40, "tmb_overdeveloped": 0.12, "tmb_underdeveloped": 0.07,
    "poor_wash_control": 0.10, "poor_wash_both": 0.06, "weak_coating": 0.08,
    "edge_evaporation": 0.07, "reader_gradient": 0.05, "sloppy_pipetting": 0.05,
}


def conditions_for(scenario: str, plate_kind: str,
                   rng: random.Random | None = None) -> PlateConditions:
    """Build one plate's conditions from a named scenario, with jitter."""
    if scenario not in SCENARIOS:
        raise KeyError(f"unknown scenario {scenario!r}; "
                       f"known: {sorted(SCENARIOS)}")
    base = PlateConditions(label=scenario)
    for key, value in SCENARIOS[scenario][plate_kind].items():
        setattr(base, key, value)
    if rng is not None:
        # No two plates are identical even when nothing went wrong.
        base.tmb_minutes *= math.exp(rng.gauss(0, 0.08))
        base.coating *= math.exp(rng.gauss(0, 0.12))
        base.blank_od *= math.exp(rng.gauss(0, 0.15))
        base.cv *= math.exp(rng.gauss(0, 0.10))
    return base


# --- the optics -------------------------------------------------------------

DEVELOPMENT_RATE = 0.055       # per unit bound per minute; sets the time scale
CROSS_REACTIVITY = 0.85        # how much of a clone's stickiness the control
                               # antigen also picks up


def optical_density(bound: float, conditions: PlateConditions) -> float:
    """TMB develops toward the reader ceiling, never past it.

    This is the single most important non-linearity in the whole simulation.
    A plate left too long does not get proportionally brighter -- it runs out
    of headroom, and every well above a certain amount of bound antibody ends
    up at the same reading. That is why an over-developed plate loses its hits
    rather than exaggerating them.
    """
    developed = 1.0 - math.exp(-DEVELOPMENT_RATE * max(bound, 0.0)
                               * conditions.tmb_minutes)
    return conditions.blank_od + (conditions.ceiling_od - conditions.blank_od) * developed


def blank_od(conditions: PlateConditions, well: str, rng: random.Random,
             plate_format: int = 96) -> float:
    """An empty well: buffer and substrate, no antibody.

    Not zero, and not noiseless. These wells are the plate's own measured
    background, which is more useful than any assumed blank -- see
    `plate_features`.
    """
    value = conditions.blank_od * math.exp(rng.gauss(0, 0.12))
    if conditions.edge_effect and wellmod.is_edge(well, plate_format):
        value += conditions.edge_effect * rng.uniform(0.5, 1.0) * 0.4
    return round(max(value + rng.gauss(0, 0.008), 0.0), 4)


def well_od(behaviour: CloneBehaviour, well: str, on_target: bool,
            conditions: PlateConditions, rng: random.Random,
            plate_format: int = 96) -> float:
    """One well's reading.

    A small fraction of wells simply fail -- a bubble, a missed dispense, no
    detection antibody -- and read at background whatever is in them. This
    matters far more than its rate suggests, because the failures are
    INDEPENDENT between the two plates. A failed CONTROL well makes an
    ordinary clone's ratio explode, which is where real false positives come
    from; a failed TARGET well hides a genuine binder. Without it, a ratio
    threshold is never wrong and every caller looks perfect.
    """
    if rng.random() < conditions.well_failure_rate:
        return blank_od(conditions, well, rng, plate_format)
    if not behaviour.expressed:
        bound = 0.0
    else:
        specific = behaviour.affinity if on_target else behaviour.control_affinity
        residual = behaviour.stickiness * (1.0 if on_target else CROSS_REACTIVITY)
        # Washing removes non-specific binding, not specific binding.
        residual /= max(conditions.wash_quality, 0.05)
        bound = conditions.coating * (specific + residual)

    value = optical_density(bound * math.exp(rng.gauss(0, conditions.cv)),
                            conditions)

    row = ord(wellmod.normalize(well)[0]) - ord("A")
    col = int(wellmod.normalize(well)[1:]) - 1
    value += conditions.row_drift * row + conditions.col_drift * col
    if conditions.edge_effect and wellmod.is_edge(well, plate_format):
        value += conditions.edge_effect * rng.uniform(0.5, 1.0)

    value += rng.gauss(0, 0.012)                       # read noise
    return round(min(max(value, 0.0), conditions.ceiling_od), 4)


# --- a pair -----------------------------------------------------------------

@dataclass
class PlatePair:
    """One experiment: the same clones on target and on control."""
    pair_id: str
    scenario: str
    target: pd.DataFrame          # well, clone_id, od450
    control: pd.DataFrame
    behaviours: dict[str, CloneBehaviour]
    target_conditions: PlateConditions
    control_conditions: PlateConditions

    def truth(self) -> pd.DataFrame:
        """What the classifier is not allowed to see."""
        target = self.target[self.target["role"] == "sample"]
        control = self.control[self.control["role"] == "sample"]
        merged = target.merge(control, on=["well", "clone_id", "role"],
                              suffixes=("_target", "_control"))
        merged["ratio"] = (merged["od450_target"]
                           / merged["od450_control"].clip(lower=1e-6))
        merged["delta"] = merged["od450_target"] - merged["od450_control"]
        merged["is_binder"] = [self.behaviours[c].is_binder
                               for c in merged["clone_id"]]
        merged["is_control_binder"] = [self.behaviours[c].is_control_binder
                                       for c in merged["clone_id"]]
        merged["expressed"] = [self.behaviours[c].expressed
                               for c in merged["clone_id"]]
        merged["pair_id"] = self.pair_id
        merged["scenario"] = self.scenario
        return merged

    def long(self) -> pd.DataFrame:
        """Both plates stacked, with the label a classifier must predict."""
        frames = []
        for kind, frame in (("target", self.target), ("control", self.control)):
            part = frame.copy()
            part["plate_kind"] = kind
            part["pair_id"] = self.pair_id
            part["scenario"] = self.scenario
            frames.append(part)
        return pd.concat(frames, ignore_index=True)


def simulate_pair(clone_ids: list[str], *, seed: int = 0,
                  scenario: str | None = None, pair_id: str | None = None,
                  plate_format: int = 96,
                  wells_order: list[str] | None = None,
                  **clone_kwargs) -> PlatePair:
    """One target/control pair over the same plate map.

    Fewer clones than wells is normal, and the plate still comes back with
    every well in it -- because that is what a plate reader hands you. The
    unused wells are blanks, marked `empty`, and they matter twice over: they
    drag every distribution statistic toward the blank if you forget they are
    there, and they are a free measurement of this plate's own background if
    you remember.
    """
    rng = random.Random(seed)
    if scenario is None:
        names, weights = zip(*SCENARIO_WEIGHTS.items())
        scenario = rng.choices(names, weights=weights, k=1)[0]

    all_positions = wells_order or wellmod.all_wells(plate_format, "column")
    used = all_positions[:len(clone_ids)]
    clone_of = dict(zip(used, clone_ids))
    behaviours = draw_clones(clone_ids, rng, **clone_kwargs)

    frames, conditions = {}, {}
    for kind, on_target in (("target", True), ("control", False)):
        cond = conditions_for(scenario, kind, rng)
        conditions[kind] = cond
        rows = []
        for well in all_positions:
            clone = clone_of.get(well)
            if clone is None:
                rows.append({"well": well, "clone_id": None, "role": "empty",
                             CHANNEL: blank_od(cond, well, rng, plate_format)})
            else:
                rows.append({"well": well, "clone_id": clone, "role": "sample",
                             CHANNEL: well_od(behaviours[clone], well,
                                              on_target, cond, rng,
                                              plate_format)})
        frames[kind] = pd.DataFrame(rows)

    return PlatePair(pair_id=pair_id or f"PAIR-{seed:06d}", scenario=scenario,
                     target=frames["target"], control=frames["control"],
                     behaviours=behaviours,
                     target_conditions=conditions["target"],
                     control_conditions=conditions["control"])


# How full a plate is. Most runs fill it, but a campaign tail, a repeat of a
# few clones, or a pilot leaves wells empty -- and a model that has only seen
# full plates mis-reads those badly (see `plate_features`).
FILL_MIX = [(1.00, 0.55), (0.90, 0.10), (0.75, 0.12), (0.50, 0.13),
            (0.25, 0.07), (0.12, 0.03)]


def simulate_many(clone_ids: list[str], n_pairs: int = 100, *,
                  seed: int = 0, scenario: str | None = None,
                  vary_fill: bool = True, plate_format: int = 96,
                  **kwargs) -> list[PlatePair]:
    """A campaign's worth of plates, each with its own scenario and clones.

    Clone behaviour is redrawn per pair on purpose: a training set built from
    one set of clones teaches a classifier those clones, not the assay. Fill
    fraction varies for the same reason -- a training set of full plates
    teaches it that a pile of blanks means something.
    """
    out = []
    for i in range(n_pairs):
        pair_seed = seed * 100_000 + i
        subset = clone_ids
        if vary_fill:
            chooser = random.Random(pair_seed ^ 0x5EED)
            fraction = chooser.choices([f for f, _ in FILL_MIX],
                                       weights=[w for _, w in FILL_MIX])[0]
            n = max(4, int(round(plate_format * fraction)))
            subset = clone_ids[:min(n, len(clone_ids))]
        out.append(simulate_pair(subset, seed=pair_seed, scenario=scenario,
                                 pair_id=f"PAIR-{seed:04d}-{i:05d}",
                                 plate_format=plate_format, **kwargs))
    return out


# --- features a classifier could use ---------------------------------------

def plate_features(values: pd.Series, roles: pd.Series | None = None) -> dict:
    """Shape of one plate's signal distribution, with no labels in it.

    Deliberately scale-aware AND scale-free: absolute level alone separates
    clean plates and fails on over-developed ones, while shape alone fails on
    the weak-coating case. A classifier needs both to survive the mixture.

    **`roles` is not optional in practice.** Computed over all 96 wells of a
    half-filled plate, every statistic here is wrong -- not noisy, wrong. The
    blanks take over the median, so on measured data with identical biology
    and 40 samples instead of 96:

        median            0.203 -> 0.061
        p99_over_median   6.46  -> 15.04
        top_over_median   4.40  -> 10.67
        frac_near_blank   0.29  -> 0.72

    A model trained on full plates then reads a half-filled CONTROL plate as
    a target plate, because a pile of blanks below a normal signal looks
    exactly like a hit population above a low background. So the shape
    statistics are computed over sample wells only.

    The empty wells are then an asset rather than a nuisance: they are this
    plate's own background, measured on this plate on this day, so
    `signal_over_blank` is a real ratio rather than one against an assumed
    constant. `fill_fraction` is reported because a nearly empty plate has
    genuinely less evidence in it and a model should be able to say so.
    """
    values = pd.Series(values).astype(float).reset_index(drop=True)
    if roles is None:
        samples, blanks = values, pd.Series(dtype=float)
    else:
        roles = pd.Series(roles).reset_index(drop=True)
        samples = values[roles == "sample"]
        blanks = values[roles != "sample"]

    v = samples.astype(float)
    n = len(v)
    if not n:
        return {}
    mean, sd = float(v.mean()), float(v.std(ddof=0))
    q = {p: float(v.quantile(p)) for p in (0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)}
    median = q[0.5]
    mad = float((v - median).abs().median())
    top = v.nlargest(max(1, n // 10)).mean()
    return {
        "mean": mean, "sd": sd, "cv": sd / mean if mean else 0.0,
        "min": float(v.min()), "max": float(v.max()),
        "median": median, "mad": mad,
        "iqr": q[0.75] - q[0.25],
        "p90_over_median": q[0.9] / median if median else 0.0,
        "p99_over_median": q[0.99] / median if median else 0.0,
        "top_decile_mean": float(top),
        "top_over_median": float(top) / median if median else 0.0,
        "skew": float(v.skew()) if n > 2 else 0.0,
        "kurtosis": float(v.kurtosis()) if n > 3 else 0.0,
        "frac_over_1": float((v > 1.0).mean()),
        "frac_over_2": float((v > 2.0).mean()),
        "frac_near_ceiling": float((v > 3.6).mean()),
        "frac_near_blank": float((v < 0.15).mean()),
        "range_over_median": (float(v.max()) - float(v.min())) / median if median else 0.0,
        # what the plate itself says its background is, and how much of the
        # plate was actually used
        "n_samples": float(n),
        "fill_fraction": float(n / len(values)) if len(values) else 1.0,
        "blank_median": float(blanks.median()) if len(blanks) else float("nan"),
        "blank_sd": float(blanks.std(ddof=0)) if len(blanks) > 1 else float("nan"),
        "signal_over_blank": (median / float(blanks.median())
                              if len(blanks) and blanks.median() else float("nan")),
        "top_over_blank": (float(top) / float(blanks.median())
                           if len(blanks) and blanks.median() else float("nan")),
    }


def infer_roles(values, n_samples: int | None = None) -> pd.Series:
    """Which wells held samples, when the plate map is not to hand.

    **You cannot tell an empty well from a failed sample well by OD alone.**
    A well with no antibody in it and a well whose clone did not express read
    the same, because they are the same thing optically. That is a fact about
    the assay, not a gap in the algorithm, and every threshold rule tried here
    ran into it: splitting on the largest multiplicative gap, on a multiple of
    the blank, and on the tightness of the low population all landed between
    74% and 90% per-well agreement, and all of them made FULL plates worse by
    calling an ordinary low tail empty.

    So there is no silent guessing. Two honest inputs, in order:

    * the plate map, from `plates.load(...).wells` -- exact, and the pipeline
      already has it;
    * `n_samples`, which whoever ran the plate knows. The n highest wells are
      taken as the samples, which measured 87-100% per-well agreement across
      fills, and exactly 100% on a full plate.

    With neither, every well is called a sample and the caller is expected to
    record that the map was unknown. That is not a good answer, but it is a
    visible one.
    """
    v = pd.Series(values).astype(float).reset_index(drop=True)
    n = len(v)
    if n_samples is None or n_samples >= n:
        return pd.Series(["sample"] * n)
    if n_samples <= 0:
        return pd.Series(["empty"] * n)
    occupied = set(v.nlargest(n_samples).index)
    return pd.Series(["sample" if i in occupied else "empty" for i in range(n)])


def feature_table(pairs: list[PlatePair]) -> pd.DataFrame:
    """One row per plate: its features, and whether it was the target plate."""
    rows = []
    for pair in pairs:
        for kind, frame in (("target", pair.target), ("control", pair.control)):
            rows.append({"pair_id": pair.pair_id, "scenario": pair.scenario,
                         "plate_kind": kind,
                         "is_target": kind == "target"}
                        | plate_features(frame[CHANNEL],
                                         frame.get("role")))
    return pd.DataFrame(rows)


# --- looking at it ----------------------------------------------------------

def heatmap(pairs: list["PlatePair"], out, n: int = 4, plate_format: int = 96):
    """Example pairs as plate heatmaps, target above control.

    The point of looking is to catch a simulation that is obviously wrong --
    a control plate brighter than its target everywhere, a hit population
    that is really a gradient. Numbers hide that; a plate picture does not.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from ..core import wells as _w

    n = min(n, len(pairs))
    n_rows, n_cols = _w.dims(plate_format)
    fig, axes = plt.subplots(2, n, figsize=(3.6 * n, 5.4), dpi=150,
                             squeeze=False)
    ceiling = max(p.target_conditions.ceiling_od for p in pairs[:n])

    for j, pair in enumerate(pairs[:n]):
        for i, (kind, frame) in enumerate((("target", pair.target),
                                           ("control", pair.control))):
            grid = np.full((n_rows, n_cols), np.nan)
            for _, row in frame.iterrows():
                well = _w.normalize(row["well"])
                grid[ord(well[0]) - ord("A"), int(well[1:]) - 1] = row[CHANNEL]
            ax = axes[i][j]
            im = ax.imshow(grid, vmin=0, vmax=ceiling, cmap="magma")
            ax.set_xticks([]); ax.set_yticks([])
            hits = sum(b.is_binder for b in pair.behaviours.values())
            ax.set_title(f"{kind}" + (f"  ({hits} binders)" if i == 0 else ""),
                         fontsize=9, color="#0b0b0b")
            if i == 0:
                ax.text(0.5, 1.28, pair.scenario, transform=ax.transAxes,
                        ha="center", fontsize=9, color="#52514e")
    fig.colorbar(im, ax=axes, shrink=0.7, label="OD 450 nm")
    fig.suptitle("SIMULATED ELISA — not real data", fontsize=11,
                 color="#52514e", y=1.02)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    from pathlib import Path as _P
    return _P(out)
