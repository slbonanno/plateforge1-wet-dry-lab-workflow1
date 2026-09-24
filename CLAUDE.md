# plateforge — working agreement

Read this before changing anything. These rules exist so the repo still makes
sense months from now when someone (probably me) drops back into one module
with no memory of the rest.

## What this is

A toolkit for antibody discovery workflows, built in independent modules:

| Module | Owns |
|---|---|
| `core` | IDs, well coordinates, paths, stores, bulk tables, the artifact ledger, registries |
| `library` | Sequence acquisition (OAS), classification, in-silico reformatting |
| `assay` | Experiment definition, plate layout, reader ingest, analysis, hit calling, picking |
| `reagents` | Reagent catalog, lots, usage, and the rules layer that raises assay caveats |
| `emit` | Files for external systems: liquid handler worklists, sequencing forms, reports |

## Non-negotiable rules

1. **Modules never import each other.** `library`, `assay`, `reagents`, and
   `emit` may all import `core`. They may not import one another. If module B
   needs something from module A, it takes an `obj_id` and resolves it through
   `core.artifacts`. If you feel the urge to add a cross-module import, that is
   a signal the artifact contract needs a new type, not a new import.

2. **IDs are minted once, in `core.ids`, and never recomputed.** Every object
   type has a registered prefix. Adding a type is one `register_prefix()` call.
   Never build an ID with an f-string outside `core.ids`.

3. **Wells are always `A01`, never `A1`.** Every entry point calls
   `core.wells.normalize()`. No other module parses a well string.

4. **No module hardcodes a path.** Databases come from
   `core.stores.connect(name)`; every other path comes from `core.paths`.
   Adding a store is one `register_store()` call plus a `.sql` file in
   `core/schema/`. Large tables go to parquet via `core.bulk`.

5. **Extension happens through registries, not edits.** New reader, construct
   format, hit caller, control resolver, emitter, or reagent rule = a new file
   plus one `@REGISTRY.register("name")` decorator. You should almost never
   need to modify an existing handler to add a new one.

6. **Never write a file emitter without a real example file in `fixtures/`.**
   Format guesses are how this project fails. See `docs/formats/README.md`.
   The same applies to readers: `docs/formats/oas.md` records what was checked
   against a real shard and what is still assumed. When an external source is
   unreachable, write down what was tried and what it returned — a future
   session should not have to rediscover that OPIG serves 403.

7. **Unanticipated fields go in `meta_json`, not a new column and not prose.**
   Query them with `json_extract`. When a key recurs enough to matter, promote
   it to a real column and backfill. `notes` is for human sentences only.

8. **Config is JSON, validated against a schema.** Experiment definitions must
   be machine-writable, because eventually an agent writes them. JSON has no
   comments — use a `"_comment"` key.

## Conventions

- Python 3.11+, `src/` layout, install with `pip install -e ".[dev]"`
- Data lives outside the repo: set `PLATEFORGE_DATA`, defaults to `./data/`
- `pytest` from the repo root; every module needs tests that run without
  network access or instrument files
- Type hints on public functions; `from __future__ import annotations` at top
- Deterministic functions do not touch the ledger; the caller registers results

9. **Germline calls are parsed only in `library.germlines`.** OAS `v_call` can
   be a comma-separated tie; treating it as a gene name directly produces silent
   garbage categories.

10. **Artifacts are never modified.** A re-run produces a new artifact and
   `artifacts.supersede(old, new)` records the replacement. `artifacts.current()`
   returns what has not been replaced.

## Decision log

Anything structural gets a dated note in `decisions/`. Read those before
proposing a redesign — the reasoning is usually there.

**`decisions/0004-open-questions.md` is a living register of what is
deliberately undecided.** Open questions belong there, not in code comments and
not in someone's head. Read it when returning to the project after a gap. Do not
resolve one silently — write the decision record.

## README.md is author-maintained — do not edit it

`README.md` belongs to the repository owner and is written by hand. Never edit
it, and never include it in a patch. When a change alters what the project can
do or how it is run, say so in the commit message and in the relevant module
README under `src/plateforge/*/README.md`, which are fair game.

## Current state

`core`, `library`, `assay` and `emit` are built and tested (311 tests).
`reagents` is still a stub.

`library` covers OAS ingest, the sequence pool, diversity-aware sampling,
germline-anchored alignment (`library.msa`, decision 0014), run bundles
(`library.selection` + `core.runs`, decision 0015) and the ordering path
(`library.vector` / `cloning` / `vendors` / `ordering`, decisions 0016-0018,
0020). It
is tested entirely against synthetic OAS-format units (`library.synth`), so the
suite runs with no network and no downloaded data. A real OAS fixture is still
outstanding (Q9).

Generated figures go under `figures/` (gitignored), split into `synthetic/`
and `real/`. Committed README figures live in `docs/figures/`. Never mix them:
a synthetic figure presented as real is the most damaging mistake this repo
can make.

Germline references come from the IMGT tables bundled with `anarci`
(`pip install -e ".[imgt]"`); IMGT's own download is unreachable from here and
no longer on the critical path (decision 0017).

`docs/pipeline.md` is the two-command tour of what exists end to end.

`assay` is the plate lineage layer (decision 0019): **the clone is tracked,
plates are containers, steps are registry entries.** The working chain runs
transfect → harvest → beads → IP → elute → magnet → BLI → normalise, and a
well traces back to the plate it was ordered on.

`emit` writes instrument worklists, one registered emitter per instrument,
each declaring how well its format is known — `verified`, `documented` or
`sketch` (decision 0021). Nothing is `verified` but our own format, and a test
asserts that. Tecan `.gwl` and Opentrons protocols are `documented`; INTEGRA,
Hamilton, Beckman and Agilent refuse and say what would unblock them.
`fixtures/integra/` and `fixtures/octet/` list exactly what to capture.

`core.lab` holds which instruments this lab has, so the agent's "which robot
do you have?" is asked once and stored.

`assay.elisa` simulates paired target/control ELISA plates (decision 0022),
stamped SIMULATED everywhere. It is deliberately hard to classify: panels vary
from zero binders upward and the control antigen has real binders of its own.

`assay.readers` parses plate exports by finding the grid by shape rather than
by vendor layout, and `assay.assign` matches grids to `PLT` artifacts with
explicit precedence — a model guess always needs confirmation (decision 0023).
A 96-well plate rarely holds 96 samples, and features computed over empty
wells are wrong rather than noisy, so the plate map or the sample count has to
be supplied.

`scripts/run_all.py` reproduces every claim on the machine running it.

The next piece of work is hit calling, then `reagents`, and closing Q15/Q7
with real files — including a real Gen5 export.

Known-undecided areas are listed in `decisions/0004-open-questions.md`. Schemas
in `core/schema/` are expected to gain columns; that is normal. What is expensive
to change, and therefore stable, is the meaning of an ID and the artifact
contract itself.
