# 0002 — Several named SQLite stores, not one database

Date: 2026-09-12
Status: accepted

## Context

The sequence library is a slow-changing, reusable resource sampled across many
unrelated experiments. Assay data is per-experiment and churns. Reagent data is
reference data with its own lifecycle. Forcing these into one file couples
things with genuinely different update patterns and backup needs.

## Decision

Named stores registered in `core.stores`, each with its own `.sql` schema and
its own file under `PLATEFORGE_DATA`. The ledger is itself a store and spans
all of them. Cross-store queries use SQLite `ATTACH`, or join in pandas on IDs.

## Consequences

- Adding a store later is one `register_store()` call plus a schema file.
- Foreign keys cannot span stores; referential integrity across stores is the
  ledger's job, not SQLite's.
- The library can be copied, versioned, or shared without dragging assay data.
