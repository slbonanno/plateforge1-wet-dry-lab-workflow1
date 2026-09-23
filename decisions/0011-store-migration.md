# 0011 — Additive schema migration, and keeping the two halves of a store in step

Date: 2026-09-22
Status: accepted
Closes: Q8

## Context

Two failures on the same real run.

`CREATE TABLE IF NOT EXISTS` does not alter an existing table, so adding
`n_ambiguous` to the schema broke every database that already held data:
`table sequences has no column named n_ambiguous`. CLAUDE.md already said
schemas are expected to gain columns, but nothing implemented it.

Then, deleting the SQLite file to work around that left the parquet half
untouched. Ingest guards against re-ingest by checking seq_ids in SQLite, saw
an empty index, and appended a second copy of everything. The only visible
symptom was a per-germline count that was exactly double.

## Decision

`stores.migrate()` compares the columns a schema declares against those a table
has and issues `ALTER TABLE ADD COLUMN` for the difference. It runs on every
`connect()`, after tables are created but before indexes, since an index on a
not-yet-added column would fail.

**Additive only.** Renames, type changes and drops are not attempted: those
need a considered migration, and silently rewriting a table is how data
disappears. Constraints SQLite cannot add this way (PRIMARY KEY, UNIQUE) are
stripped, and NOT NULL without a default is dropped rather than failing on a
populated table.

`bulk.write(key=...)` deduplicates on append, so the parquet half cannot double
even when the index that normally guards it is missing. `pool.consistency()`
reports drift between the two halves and `pool.repair()` drops orphaned parquet
rows — never index rows, which are the record of what was ingested.

## Consequences

- Adding a column is a schema edit again, as intended.
- A store survives losing its index; re-ingest restores it without duplication.
- Drift is detectable rather than silent. The script prints it.
- Non-additive migrations remain unsolved and will need a real versioned
  migration story when one is first required.
