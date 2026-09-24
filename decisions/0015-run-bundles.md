# 0015 — Every selection writes a run bundle

Date: 2026-09-23
Status: accepted

## Context

A set of 96 sequences is not a result. The result is *these 96 and not
others*, and that answer lived in terminal scrollback: the study, the filter
funnel, the sampler's seed, the diversity report. Re-running the script with a
different `--limit` produced a different 96 with nothing to distinguish the two
afterwards.

The alignments had the same problem in a smaller way — they were written to
`figures/`, which is gitignored and shared across runs, so the second run
overwrote the first run's evidence.

## Decision

Every selection writes a **bundle** under `outputs/<LIB id>/`:

```
manifest.json           source, filter funnel, selection, diversity, alignments
selected.csv            the sequences, in rank order
pool_summary.csv        what they were drawn from
filter_funnel.csv       how the source narrowed to the pool
alignment_<gene>.png    one per V gene
alignment_<gene>.fasta  the same alignment, openable elsewhere
alignment_<gene>.aln.txt
input_summary.png
summary.md              a README block describing this run, generated
```

`core.runs.Bundle` writes it; `library.selection.run()` fills it. The bundle
directory is keyed on the `LIB` id, so a re-run mints a new id and gets its own
directory. Nothing overwrites an earlier run's answer — consistent with rule 10
(artifacts are never modified).

The manifest records `code_version`, so a six-month-old bundle names the commit
that produced it.

`summary.md` is generated from the closed manifest, never written by hand. A
hand-written description of a query goes stale the first time the query changes
and nobody notices; this one cannot describe a run other than the one it sits
in.

## Why the manifest and not just the ledger

`core.artifacts` already records lineage, and it stays the index. But the
ledger is a database that needs the repo and the right `PLATEFORGE_DATA` to
read. The manifest is a file sitting next to the figures it explains, readable
by anyone who is handed the directory — including the person who finds it on a
shared drive in a year. The ledger is for queries; the manifest is for the
directory to be self-describing. They are allowed to overlap.

## Consequences

- `scripts/fetch_oas.py` and `scripts/order.py` print the bundle path.
- `core.runs.describe(obj_id)` summarises an unexplained directory.
- `oas.Filter.as_dict()` exists so the filter can go in the manifest verbatim.
- Alignments per V gene are produced by `figures.alignments_per_gene`, which is
  the per-family deliverable, rather than one figure for the busiest gene.
