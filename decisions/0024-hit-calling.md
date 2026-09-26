# 0024 — Hit calling, and letting a caller refuse

Date: 2026-09-24
Status: accepted
Extends: 0022, 0023

Every number in this record comes out of `scripts/call_hits.py` and can be
re-run on the machine that owns the data:

```
python scripts/call_hits.py --pairs 600 --seed 17
```

The plates are simulated. That is the only reason the numbers exist at all —
a real plate never comes with an answer key, so a calling rule can only be
calibrated where the truth is known and then carried across. Every threshold
below should be re-measured against real plates once we have them (Q7).

## A hit is a well that beats its own control well

Not a well that is high. The pair is the unit (0022), so anything that
affects a clone equally on both plates — how much IgG was prepped, how sticky
the clone is, how that row was pipetted — cancels before any threshold is
applied. What does not cancel is anything that hit the two *plates*
differently, and that is the whole reason a caller cannot be one number.

## Four callers, because they fail differently

One registry entry each (rule 5); adding a fifth moves nothing.

| caller | rule | characteristic failure |
|---|---|---|
| `fixed_ratio` | target/control ≥ 3, target ≥ 0.2 OD | the 0.2 floor is wrong on a plate whose blank is 0.02, and wrong again at 0.25 |
| `ratio_over_blank` | same ratio, floor = 4× the plate's own blanks | needs empty wells, so it is unavailable on a full plate |
| `robust_z` | robust z of (target − control), MAD not SD | finds "outliers" on a plate where nothing bound |
| `ratio_and_z` | both must agree | costs recall |

Measured over 600 pairs, refusals respected:

| caller | precision | recall | F1 |
|---|---|---|---|
| `robust_z` | 0.745 | 0.742 | **0.632** |
| `ratio_over_blank` | 0.869 | 0.646 | 0.626 |
| `fixed_ratio` | 0.843 | 0.646 | 0.621 |
| `ratio_and_z` | 0.875 | 0.626 | 0.612 |

They are within noise of each other on F1 and *not* within noise on the
trade-off, which is the actual finding: `robust_z` finds about 15% more
binders and is wrong about twice as often. The default is `ratio_and_z`
because a false positive here costs a rescreen of a clone that does not bind,
and the next step after a call is expensive. A cherry-pick that is short is
easier to live with than one that is dirty.

## The simulator had to get worse first

The first evaluation returned precision of exactly 1.000, which is not a
result, it is a bug report about the data. Nothing in the simulator could
produce a *spuriously* high ratio: every well's reading was a faithful
function of what the clone did, so a well only cleared 3× when something
really bound.

Real plates have wells that simply fail — a bubble, a missed dispense, no
detection antibody — and read at background whatever is in them.
`well_failure_rate = 0.02` (0022's `ASSUMED`) adds those, and the property
that matters is that **the failures are independent between the two plates**.
A failed *control* well makes an ordinary clone's ratio explode, which is
where real false positives come from. A failed *target* well hides a real
binder. With that in, precision landed at 0.81–0.88 depending on the caller,
and the callers started to separate from each other.

A simulation that cannot generate the mistake you are trying to measure will
tell you your rule never makes it.

## The MAD-is-zero trap

`robust_z` fell back to the plate SD when the MAD came out at exactly zero,
which happens when over half the wells read identically — a plate where
almost nothing expressed. That fallback is the one thing that must not
happen: the SD is inflated by the very wells being tested, so a plate of 80
blanks and 16 strong binders gets a spread wide enough to call **none** of
them. It now falls back to the spread of the middle half, which the outliers
cannot reach, floored at the reader's own resolution (0.001 OD). A test
pins the 80/16 case.

## A caller is allowed to say no, and the evidence says which plates

`callable` / `degraded` / `uncallable`, returned alongside the per-well calls.
`uncallable` is a normal outcome, not an error. Anything that is not
`callable` reaches a human before it reaches a cherry-pick.

**The first design was backwards.** The obvious refusal is saturation: TMB
develops toward the reader ceiling, so an over-developed well stops carrying
information and the ratio collapses toward 1. Measured, saturation goes with
*better* calls, on plates that had at least one binder:

| saturated target wells | F1 | plates |
|---|---|---|
| <2% | 0.604 | 258 |
| 2–10% | 0.817 | 230 |
| 10–20% | 0.905 | 47 |
| >20% | 0.915 | 4 |

(`ratio_over_blank`, plates with at least one binder. The top row is four
plates and carries no weight on its own; the trend across the first three,
which are 535 plates between them, is the finding.)

The reasoning was right about the mechanism and wrong about the consequence.
A well only reaches the ceiling when something really bound, so saturation is
mostly a marker of a plate with strong binders on it. What is genuinely lost
is the ability to **rank** saturated wells against each other, which matters
when picking the top 8 of 20 and not at all when asking which wells bound.
So saturation is now **flagged per well and reported in the reasons, never
refused on**.

What does predict whether calling works is dynamic range — the top decile of
the target plate over its background:

| top-decile / background | F1 | precision | recall | plates |
|---|---|---|---|---|
| <6× | 0.450 | 1.000 | 0.417 | 4 |
| 6–12× | 0.463 | 0.764 | 0.423 | 52 |
| 12–30× | 0.611 | 0.836 | 0.569 | 186 |
| >30× | 0.843 | 0.951 | 0.778 | 297 |

Hence `uncallable_dynamic_range = 10`, `degraded_dynamic_range = 25`. The
resulting tiers, scored with nobody allowed to refuse, so that the refusal
can be judged on what it is throwing away:

| verdict | F1 | true positives | false positives | plates |
|---|---|---|---|---|
| `callable` | 0.812 | 3932 | 134 | 350 |
| `degraded` | 0.478 | 525 | 76 | 191 |
| `uncallable` | 0.242 | 35 | 7 | 59 |

(One caller, `ratio_over_blank`, not pooled over four — pooled, every plate is
counted once per caller and the totals read as plate counts while being row
counts.)

On the refused plates `ratio_over_blank` would have found 35 real binders and
invented 7; `robust_z`, over the same 59 plates, would have found 50 and
invented 64. Refusing costs real hits. It is still right: those plates should
be re-read, and a re-read recovers the 35 without the 7 — and without
`robust_z`'s 64, which is the number that would actually have been paid,
since a refused plate is exactly where the loosest caller looks best.

`degraded` is deliberately not a refusal. Its calls are measurably worse and
still worth having, and the verdict travels with them so nobody has to guess.

## Why not a trained classifier

It would fit. The thresholds are interpretable, they carry a reason string a
lab member can act on ("re-read this plate sooner"), and they do not need
retraining when the antigen changes. A model here would be calibrated on
simulated plates and quietly wrong on real ones, which is the failure this
repo is built to avoid. Revisit once real reads exist.

## What real data would change

The thresholds, probably all of them. Specifically:

- **A real Gen5 export** (`fixtures/gen5/`) fixes the background estimate.
  We currently take the lowest decile when blank wells are not marked, which
  is an assumption about how plates get laid out.
- **A plate with a known positive control** would replace `min_signal = 0.2`,
  which is the weakest number here — it is the one people reach for, not one
  we measured.
- **A plate read twice, before and after over-development**, would test the
  saturation finding directly instead of inferring it from the simulation's
  own optics.

## Figures

`assay.figures` draws what a table describes badly: `elisa_pair` (target,
control, and the difference on a diverging scale), `hit_calls` (which wells
were found, missed and invented, plus the threshold as a line on a log-log
scatter), `clone_journey`, `pipeline_map`. Colour comes from `core.style`,
which is the single palette for this repo and is CVD-validated there.

Two rules those figures follow that are easy to get wrong, and were:

**Bar length is linear, always.** `pipeline_map` spans six orders of
magnitude and originally log-scaled the bar widths. That is the same
arithmetic and not the same picture — a bar reads as proportional to its
value, so 96 looked like two thirds of 616,809. It is now a log **axis**,
with ticks, so the reader sees the fall-off instead of being walked past it.

**An empty well is not a weak well.** Painted on the same ramp as the
samples, a 12-sample plate looked like a 96-sample plate that mostly failed.
Non-sample wells are NaN and the colormap paints them as "nothing here", with
the count in the legend.
