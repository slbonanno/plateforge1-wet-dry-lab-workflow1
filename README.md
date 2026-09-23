# plateforge

This project is a predecessor to an interactive experience (agent) that helps plan experiments, generate virtual representations, track workflows and sample IDs across different (connected) experiments, reformat and process data when necessary to generate files to march along the workflow, and generate low-lift and standard figures for basic experimental workflows in the wet lab.

As an example, this project is designed for running antibody discovery workflows.

> **Status:** Data layers v0 are built, architecture is sketched out.  development is focused on a functional workflow from start to finish for standardized input.  Expanding tolerance of different inputs/formats at different steps of the workflow will be most of the focus of dev for this project.

## Basic workflow:

#### 1. obtain antibody sequences
1. antibodies are discovered or sequences are generated - an input list for an experiment should at minimum have a full aa sequence.  For this project, VH's will be sourced from OAS.
2. curate sequence list - trim if necessary

#### 2. process for ordering (cloning)
1. sequences are built into a table and given a master sequence ID - immutable and forever unique
2. [expand] generate DNA sequences to be synthesized. assume scFv with fixed VL
outs:
2A) order forms for IDT, GenScript, and Genewiz
2B) table with metadata for sequences: cloneID, clone_aa_seq, clone_DNA_seq (Hu codon opt), Frag_to_order (Gibson/HiFi/InFusion adapters), Hu_germline (query IMGT), 
2C) alignment for the clones to be synthesized
3. generate plate map for DNA sequences as expected from vendor

#### 3. generate experiment for cloning in wet lab
1. using plate map from 2), calculate master mixes required for cloning.
*Golden Gate/TypeIIs is most robust for quick-cloning:
pre-digested vector, DNA fragments (digest delivered frags in-plate), T4 ligase.
Built-in selection against wrong products: vector recircularizes without insert - recognition site regenerated and re-cut in same rxn, equilibrium is pushed toward correct ligation*
2. [expand] generate Sanger sequencing order forms (plate format) for Genewiz, ELIM, other providers
3. 


## Data architecture v0

Five modules, deliberately independent:

| Module | Responsibility |
|---|---|
| `core` | IDs, well coordinates, named stores, the artifact ledger, extension registries |
| `library` | Sequence acquisition, classification, and in-silico construct reformatting |
| `assay` | Experiment definition, plate layout, reader ingest, QC, hit calling, picking |
| `reagents` | Reagent catalog, lot tracking, and rules that raise assay caveats |
| `emit` | Worklists, sequencing forms, and reports for external systems |

Four ideas hold it together:

**Modules never import each other.** A module produces something, registers it in a ledger with a typed ID, and returns that ID. The next module resolves the ID. Provenance becomes a queryable graph rather than a call stack, and every module runs and tests independently.

**The clone is the durable entity.** IDs are minted once and carried forward. A clone stays traceable from growth plate through ELISA plate through rearray to sequencing tube, across as many physical plate changes as the campaign takes.

**Extension happens through registries, not edits.** A new instrument format, a new hit-calling strategy, a new construct scaffold, or a new reagent rule is a new file plus a one-line registration. Existing code does not move.

**Formats are verified, never guessed.** Anything that parses or writes a file for an external system requires a real example of that format, committed as a fixture and documented, before the code is written. Guessed formats are how a tool like this quietly produces wrong results.

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
- [x] `library`: OAS interface, sequence pool, diversity-aware sampling
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
