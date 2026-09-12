# assay

Experiment definition through hit picking.

**Inputs:** `LIB`/`FMT` ids, experiment JSON, reader export files
**Outputs:** `EXP`, `PLT`, `RUN`, `PCK` artifacts

Planned scope:

- experiment definition from JSON: plate format, replicates, control layout,
  cloning notes
- plate layout and clone registration
- reader adapters (registry) parsing exports into the `assay` store
- control resolution (registry): mirrored plate, on-plate wells, shared reference
- QC and hit calling (registry): fold-over-control, absolute cutoff, robust z
- plots for review
- pick sessions: prompt for N, rank, enforce diversity, emit a `PCK`
