# 0022 — Simulated ELISA data, and why it is deliberately hard

Date: 2026-09-24
Status: accepted

## Context

The analysis path — hit calling, QC, the eventual classifier — has to be
built before a real plate reader export exists. `library.synth` set the
precedent: synthetic OAS units let the whole ingest path be tested with no
network and no downloaded data. `assay.elisa` is the same idea one module
along.

The design is a **pair**: the same 96 clones on a target-antigen plate and on
a control-antigen plate (HLA-His), position for position. A clone is a hit
because it beats *its own* control well, not because it is high in absolute
terms. The pair, not the plate, is the unit of the experiment.

## The optics matter more than the statistics

OD is not linear in bound antibody. TMB develops toward the reader's ceiling:

    OD = blank + (ceiling - blank) * (1 - exp(-k * bound * minutes))

Everything interesting follows from that one line. A plate left too long does
not get proportionally brighter — it runs out of headroom, every well above
some amount of bound antibody lands at the same reading, and the plate loses
its hits rather than exaggerating them. Simulating this as `signal * factor`
would produce data on which every method works and no method transfers.

## The first version was too easy, which is a bug

A fixed hit rate gave every target plate a hit population and every control
plate none. A gradient-boosted classifier hit **99.6%**, and two single
features (`p99_over_median`, `top_over_median`) were perfect on their own.
That is not a good dataset; it is a dataset that teaches "which plate has a
right tail".

Two changes, both things that are true of real screens:

**Panels vary enormously, and many yield nothing.** The hit rate is drawn per
pair from a mixture that includes zero. Screening a naive repertoire and
getting no binders is common, and it is exactly the plate someone needs help
reading.

**Control plates have real binders too.** A few clones bind the control
antigen, and a few bind the *tag* both antigens carry, lighting up both
plates. Without them a control plate can never have a right tail and the task
collapses.

Accuracy went from 99.6% to **78%**, and the difficulty landed where it
should:

| binders on the plate | accuracy |
|---|---|
| 0 | 0.69 |
| 1–2 | 0.65 |
| 3–10 | 0.74 |
| 11–30 | 0.91 |
| 30+ | 0.89 |

| hardest scenarios | | easiest | |
|---|---|---|---|
| TMB underdeveloped | 0.59 | TMB overdeveloped | 0.94 |
| poor wash, both plates | 0.63 | poor wash, control only | 0.88 |

**Near-chance on a panel with no binders is the correct answer, not a
failure.** A target plate with nothing bound to it *is* a control plate. Any
model scoring well there is reading an artefact — a coating difference, a
scenario fingerprint — and will not transfer. This is worth stating loudly
because it is the number someone will try to optimise.

## Scenarios

Nine, weighted so most plates are fine and the confusable ones are a
minority: nominal, TMB over- and under-developed, poor wash (control only,
and both), weak coating, edge evaporation, reader gradient, sloppy pipetting.
Each changes the *shape* of the distribution, which is what a classifier has
to survive.

## Stamping

Every artifact says `SIMULATED` — in the manifest, in `summary.md`, and
across the top of the figure. CLAUDE.md already calls presenting synthetic
data as real the most damaging mistake this repo can make; that applies to
plates as much as to figures.

## Not yet decided

The parameters in `ASSUMED` are plausible magnitudes, not fits. The first
real export should be used to check the blank, the ceiling, the development
rate and the hit-rate mixture, and this record updated with what was wrong.
