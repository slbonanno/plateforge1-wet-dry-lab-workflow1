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
| Q2 | Plate schema details: well roles, replicate representation, control annotation | `assay` plate layout | Deliberately not hardened; additive columns expected |
| Q3 | When do reagent caveats fire — experiment definition, analysis, or both? | `reagents` rules layer | Probably both, attaching to different artifact types |
| Q4 | Entry point: notebook, CLI, or both | `assay` pick session | Deferred; logic stays in functions so a CLI is a thin wrapper later |
| Q5 | Experiment JSON schema fields | `assay` experiment definition | Permissive schema now, tightened as fields are discovered |
| Q6 | Paired-chain support | `library` scFv path | Unpaired heavy only; light chain fixed to IGKV1-39 (see 0007) |
| Q7 | Reader format specifics for each instrument | `assay` ingest | Blocked on real export files; see docs/formats |
| Q8 | Schema migrations for stores that already hold data | all modules | `CREATE TABLE IF NOT EXISTS` does not alter existing tables; fine while stores are disposable, needs a plan before real data accumulates |
| Q9 | Real OAS fixture | `library` ingest | Still none. Blocked on OPIG 403; the HuggingFace mirror is now the route. Capture a shard subset as a fixture on the first successful pull |
| Q12 | Column-aligning ragged real sequences | `library` alignment figure | `sequence_alignment_aa` is not fixed-width: IGHV3-23 shows 54 distinct lengths over 16,497 sequences. `ANARCI_numbering` carries per-residue IMGT positions and is the proper fix |
| Q11 | Paired-chain data | `library` scFv path | The mirror splits heavy and light into separate folders, so they are not paired. True paired units are a different layout and need their own adapter |

## Closed

See numbered decision records 0001-0003, 0005-0009.

- **Q10** (which ANARCI flags are real liabilities) — closed by 0009.

- **Q1** (diversity criteria) — closed by 0008.
- **Q1a** (which germlines) — closed by 0007.
