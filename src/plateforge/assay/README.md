# assay

Experiment definition through hit picking.

**Inputs:** `LIB`/`FMT` ids, experiment JSON, reader export files
**Outputs:** `EXP`, `PLT`, `RUN`, `PCK` artifacts

## What exists

The plate lineage layer (decision 0019). **The plate is not the sample**: a
`CLN` clone is the through-line, plates are containers it passes through, and
every well records which clone it descends from and which well it came from.

| File | Owns |
|---|---|
| `plates.py` | plates and wells; `trace`, `where_is`, `states_of`, `layout_matches` |
| `steps.py` | the step registry, and `plan` / `commit` / `record` |
| `protocol.py` | the working chain: transfect → harvest → beads → IP → elute → magnet → BLI → normalise |

`plan` and `record` are separate because they happen hours apart, and because
a human-performed step then has the same shape as a robot one — the protocol
does not fork on who did the work.

```python
from plateforge.assay import plates, steps, protocol

plan = steps.get("harvest_supernatant").plan(culture, volume_ul=500.0)
steps.commit(plan)                     # plates, transfers, a 'planned' row
steps.record(plan, operator="shivan", lots=["LOT-42"],
             deviations="column 7 short")

plates.trace(igg_plate.plate_id, "D07")     # back to the ordering plate
```

Adding a step is a new class plus `@register` (rule 5). Removing one is not
calling it. Nothing downstream depends on how many steps preceded it.

No instrument files are written yet — a step says what transfers it needs, and
`emit` turns those into a worklist once a real example of the format exists.

## Simulated ELISA

`elisa.py` generates paired target/control plates (decision 0022). **Synthetic
— never present as measurements.**

```bash
python scripts/simulate_elisa.py --pairs 600 --register 2
```

Writes a run bundle: `reads.parquet`, `truth.parquet` (with the latent
per-clone truth), `plate_features.parquet` for a classifier, a heatmap of four
example pairs, and a baseline.

The dataset is deliberately hard. A first version with a fixed hit rate was
99.6% separable and taught nothing; drawing the hit rate per panel (including
zero) and giving the control antigen real binders brought it to 78%, with
near-chance accuracy on panels that yielded nothing — which is the right
answer, since a target plate with nothing bound to it *is* a control plate.

## Reading real exports

```python
from plateforge.assay import readers, assign

grids = readers.read("scan.xlsx", "biotek_gen5")     # one per sheet
model = assign.load(assign.latest())
proposals = assign.propose(grids, barcodes={...}, model=model, n_samples=40)
assign.report(proposals)                              # what, how, how sure
assign.unresolved(proposals)                          # what needs you
```

Grids are found by **shape** — consecutive column numbers with A..H beneath —
not by any vendor's layout, so the parser survives a software update. A model
guess is never a decision: `needs_confirmation` is true for anything guessed,
and swapping target for control inverts every hit call.

**A 96-well plate rarely holds 96 samples**, and features computed over empty
wells are wrong, not noisy (decision 0023). Pass the plate map, or the sample
count. Without either the pipeline says the map was unknown rather than
guessing.

## Reproducing everything

```bash
python scripts/run_all.py            # offline stages
python scripts/run_all.py --full     # including the real OAS download
```

Writes a bundle with the code version, timings and headline numbers, so the
record of a run lives on the machine that ran it.

Still planned:

- experiment definition from JSON: plate format, replicates, control layout,
  cloning notes
- plate layout and clone registration
- reader adapters (registry) parsing exports into the `assay` store
- control resolution (registry): mirrored plate, on-plate wells, shared reference
- QC and hit calling (registry): fold-over-control, absolute cutoff, robust z
- plots for review
- pick sessions: prompt for N, rank, enforce diversity, emit a `PCK`
