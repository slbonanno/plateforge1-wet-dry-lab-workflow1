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

9. **Artifacts are never modified.** A re-run produces a new artifact and
   `artifacts.supersede(old, new)` records the replacement. `artifacts.current()`
   returns what has not been replaced.

## Decision log

Anything structural gets a dated note in `decisions/`. Read those before
proposing a redesign — the reasoning is usually there.

**`decisions/0004-open-questions.md` is a living register of what is
deliberately undecided.** Open questions belong there, not in code comments and
not in someone's head. Read it when returning to the project after a gap. Do not
resolve one silently — write the decision record.

## Current state

`core` is built and tested (13 tests). Every other module is a stub. The
immediate next piece of work is `library`: interfacing with OAS, downloading and
classifying sequences, and building a diversity-aware sampling layer.

Known-undecided areas are listed in `decisions/0004-open-questions.md`. Schemas
in `core/schema/` are expected to gain columns; that is normal. What is expensive
to change, and therefore stable, is the meaning of an ID and the artifact
contract itself.
