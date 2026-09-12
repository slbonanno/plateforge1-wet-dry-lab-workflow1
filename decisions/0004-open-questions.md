# 0004 — Open questions register

Date: 2026-09-12
Status: living document

Not everything can be decided before starting. The risk is not that questions
are open — it is that they are open and *forgotten*, and get answered by
accident in three different places.

Anything undecided goes here with a note on what is blocked by it. Review this
file when returning to the project after a gap. When a question is answered, it
gets its own numbered decision record and is struck from this list.

## Open

| # | Question | Blocks | Notes |
|---|---|---|---|
| Q1 | What makes a sampled sequence set "diverse enough"? | `library` sampling | Needs V-family spread, CDRH3 length range, pairwise identity ceiling |
| Q2 | Plate schema details: well roles, replicate representation, control annotation | `assay` plate layout | Deliberately not hardened; additive columns expected |
| Q3 | When do reagent caveats fire — experiment definition, analysis, or both? | `reagents` rules layer | Probably both, attaching to different artifact types |
| Q4 | Entry point: notebook, CLI, or both | `assay` pick session | Deferred; logic stays in functions so a CLI is a thin wrapper later |
| Q5 | Experiment JSON schema fields | `assay` experiment definition | Permissive schema now, tightened as fields are discovered |
| Q6 | Paired-chain support | `library` scFv path | Unpaired heavy only for now, fixed VL |
| Q7 | Reader format specifics for each instrument | `assay` ingest | Blocked on real export files; see docs/formats |

## Closed

See numbered decision records 0001-0003, 0005-0006.
