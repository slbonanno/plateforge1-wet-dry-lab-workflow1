# 0023 — Partial plates, reading exports, and matching grids to plates

Date: 2026-09-24
Status: accepted
Extends: 0022

## A 96-well plate rarely holds 96 samples, and that broke the features

`plate_features` computed its statistics over every well it was given. On a
half-filled plate the empty wells take over the median, and every ratio built
on that median inflates. Measured on the same plate, identical biology, 40
samples instead of 96:

| | 96 samples | 40 samples |
|---|---|---|
| median | 0.203 | 0.061 |
| p90/median | 2.60 | 6.43 |
| p99/median | 6.46 | **15.04** |
| top-decile/median | 4.40 | 10.67 |
| fraction near blank | 0.29 | 0.72 |

Not noisy — wrong. And wrong in the direction that matters: a pile of blanks
under a normal signal looks exactly like a hit population over a low
background, so a model trained on full plates reads a half-filled **control**
plate as a target plate.

Three changes:

**Plates emit every well.** `simulate_pair` returns all 96 rows with a `role`
of `sample` or `empty`, because that is what a plate reader hands you. Empty
wells read at background with real scatter, not zero.

**Features take roles** and compute shape statistics over sample wells only.

**Empty wells become an asset.** They are this plate's own background,
measured on this plate on this day, so `signal_over_blank` is a real ratio
rather than one against an assumed constant. `fill_fraction` is reported
because a nearly empty plate genuinely holds less evidence and a model should
be able to say so. On a full plate the blank features are NaN by
construction; they are filled with -1 for training, so a tree can split on
"this plate had no blanks" rather than being handed a value it cannot use.

Fill now varies across a simulated campaign, because a training set of full
plates teaches a model that a pile of blanks means something.

Accuracy moved from 78% to **75%**, and degrades with fill — 0.78 at full,
0.64 at 12% — which is correct.

## You cannot tell an empty well from a failed sample well by OD alone

This was worth establishing rather than assuming. A well with no antibody and
a well whose clone did not express are the same thing optically. Three rules
were tried for inferring occupancy from the numbers: split on the largest
multiplicative gap; threshold at a multiple of the blank; require the low
population to be tight and separated. All landed between **74% and 90%**
per-well agreement, and all made *full* plates worse by calling an ordinary
low tail empty — relative feature error on a full plate went from 0.00 with
no inference to 0.25 with an unconditional threshold.

So `infer_roles` does not guess silently. It takes the plate map when there
is one, or `n_samples` — which whoever ran the plate knows, and which
measured **87–100%** per-well agreement, exactly 100% on a full plate. With
neither, every well is a sample and the caller is told the map was unknown,
in the assignment's own reason string.

## Reading exports: find the grid by shape, not by vendor layout

`assay.readers` parses a Gen5 workbook by looking for a run of consecutive
integers on a row with `A`..`H` beneath it, enclosing a rectangle of numbers.
That description is true of every microplate export anyone produces, which
makes it more robust than a parser written against one vendor's template and
broken by the next software update. Confidence is `heuristic`, a level
`emit` does not need.

What it cannot do without a real file is read the **metadata** — which read
is 450 nm, which sheet is which plate. It collects the text above each grid
verbatim and hands it over unparsed, so nobody pretends it was understood.

`write_gen5_like` produces a workbook of that shape for testing, including
the awkward parts: a metadata preamble, a grid that does not start in column
A, blank rows between blocks. It is a fixture, not a claim about Gen5.

## A sheet is not a plate

`assay.assign` matches grids to `PLT` artifacts, best evidence first, and a
weaker source never overrides a stronger one:

| | |
|---|---|
| the user said so | always wins |
| a barcode appears in the sheet | 0.99 |
| a plate name matches | 0.8 |
| the sheet says "control" or "target" | 0.75 |
| the classifier recognises it | a suggestion |

**A model guess is never a decision.** `needs_confirmation` is true for
anything the pipeline guessed, and the reason string carries the probability,
the model's held-out accuracy and whether the plate map was known. Two
barcodes in one sheet is `unresolved` with both alternatives listed, not a
coin flip.

Swapping target and control inverts every hit call, so this is the step where
a quiet mistake is most expensive. An agent holding this should show its
reasoning and wait.

## The model is an artifact

`assign.train/save/load` persists a fitted classifier under a `DOC` id with
its feature order, held-out accuracy and what it was trained on — stamped
`SIMULATED`, with `"re-fit on real data before trusting it on real data"` in
the file. Accuracy is stored unrounded: a recorded metric that changes when
you reload it is a small dishonesty.

## Reproduction

`scripts/run_all.py` regenerates everything on the machine running it and
writes a bundle with the code version, timings and headline numbers. It
earned its place immediately by catching a NaN crash in the ELISA baseline
that the test suite did not.
