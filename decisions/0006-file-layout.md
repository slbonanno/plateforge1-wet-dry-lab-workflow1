# 0006 — Data file layout and bulk storage split

Date: 2026-09-12
Status: accepted

## Decision

All data lives under `PLATEFORGE_DATA`, never in the repository:

```
stores/            named SQLite databases
bulk/              large tables as parquet
raw/<source>/      downloaded external files, read-only after landing
outputs/<obj_id>/  everything produced for one artifact
```

Paths are constructed only in `core.paths`. Nothing else builds a data path.

SQLite holds what is queried and joined; parquet holds what is scanned in bulk.
The sequence pool is the table expected to outgrow SQLite, so the split exists
now rather than as a later migration. The pattern is a SQLite index table with
IDs and filter columns, alongside a parquet file with the full rows.

## Consequences

- Moving to an external drive is one environment variable, no code change.
- `raw/` is never edited, so any parsing bug can be re-run against the original.
- One directory per artifact means a file on disk always traces back to a ledger
  entry, and orphaned output is visible as a directory with no matching row.


## Amendment, 2026-09-23 — the data tree moved into the checkout

`PLATEFORGE_DATA` used to default to `~/plateforge-data`, set by each script.
Two problems: a project was two directories in two places, and the scripts
each restated the default, so they could drift.

The default is now `<repo>/plateforge-data/`, computed in `core.paths` from
the location of the package — not from the working directory. A cwd-relative
default means the same command writes to a different store depending on where
you ran it, and ending up with two half-populated pools is a bad afternoon.
The directory is gitignored; a pool of a few hundred thousand sequences has no
business in git.

`PLATEFORGE_DATA` still wins when set, so an external drive is one export
away. Existing data moves with `mv ~/plateforge-data <repo>/plateforge-data`.
