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
| Q2 | Replicate representation, and control annotation beyond a well role | `assay` plate layout | Wells now carry `role` (sample/control/blank/empty) and per-well lineage (0019). Replicates and control *resolution* are still open |
| Q3 | When do reagent caveats fire — experiment definition, analysis, or both? | `reagents` rules layer | Probably both, attaching to different artifact types |
| Q4 | Entry point: notebook, CLI, or both | `assay` pick session | Deferred; logic stays in functions so a CLI is a thin wrapper later |
| Q5 | Experiment JSON schema fields | `assay` experiment definition | Permissive schema now, tightened as fields are discovered |
| Q6 | Paired-chain support | `library` scFv path | Unpaired heavy only; light chain fixed to IGKV1-39 (see 0007) |
| Q20 | A real Gen5 export | `assay.readers` metadata | The grid finder works by shape and needs no fixture, but which sheet is which plate, which read is 450 nm and what the timestamps mean are unparsed until one real workbook is in `fixtures/gen5/` |
| Q7 | Reader and worklist formats per instrument | `emit` output, `assay` ingest | Structure is built (0021): every emitter declares `verified` / `documented` / `sketch`. Tecan and Opentrons are documented; INTEGRA, Hamilton, Beckman and Agilent refuse. Nothing is verified — one real file per instrument closes it, see `fixtures/*/README.md` |
| Q19 | Protocol parameters: bead and buffer volumes, incubation times, wash counts | `assay.protocol` defaults | Each step now declares its placeholders in an `assumed` dict and takes them as arguments, so nothing is hardcoded. Still needs real numbers, and a campaign protocol JSON to carry them (rule 8) |
| Q9 | Real OAS fixture | `library` ingest | Still none, and now overdue: real ingest works, so a few thousand rows should be committed as a fixture |
| Q15 | Vendor order form layouts | ordering | IDT, GenScript, Genewiz each differ, and all three keep their upload templates behind a login. `vendors.idt_eblocks` is a believed-two-column guess that refuses to emit without `allow_unverified=True` and stamps a caveat file when it does. State and capture instructions: `docs/formats/vendors.md` |
| Q16 | Vector and adapter sequences for routes other than the working one | construct design | Closed for `pcdna-scfv-vk-v1` + Golden Gate by 0016. Gibson arms are emitted from the declared vector context, which is only as good as that declaration; a real plasmid map would settle it |
| Q14 | Short CDRH3s inflate pairwise identity | `library` diversity spec | Two 5-mers differing at one position are 0.8 identical. The spec's threshold is length-blind; a length floor is the current workaround |
| Q13 | Non-additive schema migration | all modules | Additive column adds are handled (0011). Renames, type changes and drops are not, and will need a versioned migration when first required |
| Q18 | Should the picks resemble panning output rather than naive repertoire? | `library` selection | The current 96 are drawn from an unsorted naive repertoire: near-germline, unclonal, CDRH3s that no selection has enriched. That is the right input for a plumbing test and the wrong one for a mock campaign |
| Q11 | Paired-chain data | `library` scFv path | The mirror splits heavy and light into separate folders, so they are not paired. True paired units are a different layout and need their own adapter |

## Closed

See numbered decision records 0001-0003, 0005-0023.

- **Q17** (is IMGT reachable) — closed by 0017: it is not, from here, and it
  no longer matters. The tables ship inside the `anarci` package.
- **Q8** (schema migration, additive) — closed by 0011.
- **Q12** (aligning ragged sequences) — closed by 0010, then re-answered by
  0014: IMGT numbering classifies, it does not align. Ragged reads now align
  to the germline with gaps and nothing is excluded.

- **Q10** (which ANARCI flags are real liabilities) — closed by 0009.

- **Q1** (diversity criteria) — closed by 0008.
- **Q1a** (which germlines) — closed by 0007.
