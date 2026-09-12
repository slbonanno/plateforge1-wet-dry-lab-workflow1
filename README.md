# plateforge

A toolkit for running antibody discovery workflows end to end — from sequence
library through plate-based screening to the files that drive a liquid handler.

> **Status: early.** The core layer is built and tested. The domain modules are
> scaffolded with their contracts declared and their implementations pending.
> The architecture is the deliverable right now, not the biology.

## Why

Antibody discovery generates a specific kind of mess. A screening campaign
produces plates of clones, each read on an instrument that exports its own
dialect of CSV, analysed against controls that live somewhere else, filtered
down to hits, rearrayed onto new plates by a robot that wants yet another file
format, and finally submitted for sequencing on a form designed by a vendor.
Provenance leaks at every step, and by the time a clone reaches sequencing it is
often hard to reconstruct which well on which plate it came from, or what
reagents were used to call it a hit.

Most of that work is deterministic. It is also error-prone, tedious, and done by
hand at most benches. `plateforge` treats the whole path as one traceable
pipeline, with the clone — not the plate — as the entity that persists.

## Design

Five modules, deliberately independent:

| Module | Responsibility |
|---|---|
| `core` | IDs, well coordinates, named stores, the artifact ledger, extension registries |
| `library` | Sequence acquisition, classification, and in-silico construct reformatting |
| `assay` | Experiment definition, plate layout, reader ingest, QC, hit calling, picking |
| `reagents` | Reagent catalog, lot tracking, and rules that raise assay caveats |
| `emit` | Worklists, sequencing forms, and reports for external systems |

Four ideas hold it together:

**Modules never import each other.** A module produces something, registers it
in a ledger with a typed ID, and returns that ID. The next module resolves the
ID. Provenance becomes a queryable graph rather than a call stack, and every
module runs and tests independently.

**The clone is the durable entity.** IDs are minted once and carried forward. A
clone stays traceable from growth plate through ELISA plate through rearray to
sequencing tube, across as many physical plate changes as the campaign takes.

**Extension happens through registries, not edits.** A new instrument format, a
new hit-calling strategy, a new construct scaffold, or a new reagent rule is a
new file plus a one-line registration. Existing code does not move.

**Formats are verified, never guessed.** Anything that parses or writes a file
for an external system requires a real example of that format, committed as a
fixture and documented, before the code is written. Guessed formats are how a
tool like this quietly produces wrong results.

## Install

```bash
git clone https://github.com/slbonanno/plateforge1-wet-dry-lab-workflow1.git
cd plateforge
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

Data lives outside the repository:

```bash
export PLATEFORGE_DATA=~/plateforge-data
```

## Layout

```
src/plateforge/
  core/       ids, wells, paths, stores, bulk, artifact ledger, registries
  library/    sequence acquisition and reformatting
  assay/      experiments, plates, ingest, analysis, picking
  reagents/   catalog, lots, caveat rules
  emit/       worklists, forms, reports
fixtures/     real example files from external systems
docs/formats/ notes on each external format
decisions/    dated architecture decision records
config/       JSON experiment definitions and schemas
tests/
```

`CLAUDE.md` holds the working rules that keep the modules independent.
`decisions/` explains why the structure is what it is.

## Roadmap

- [x] Core: typed IDs, well normalisation, paths, SQLite stores, parquet bulk tables, artifact ledger
- [ ] `library`: OAS interface, sequence pool, diversity-aware sampling
- [ ] `library`: in-silico reformatting (scFv, VHH, VH-only, CDRH3 grafting)
- [ ] `assay`: experiment definition from JSON, plate layout
- [ ] `assay`: reader adapters, control resolution, hit calling, QC
- [ ] `assay`: review plots and interactive pick sessions
- [ ] `reagents`: catalog, lot tracking, caveat rules
- [ ] `emit`: liquid handler worklists, sequencing forms, run summaries

## Conventions worth knowing

- Wells are always zero-padded: `A01`, never `A1`
- Unanticipated fields go in a JSON `meta` column, queried with `json_extract`,
  and get promoted to real columns once they recur
- Experiment configuration is JSON validated against a schema, so an agent can
  write the same file a human does

## License

MIT
