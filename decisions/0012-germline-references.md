# 0012 — Germline references come from a source registry, and say where they came from

Date: 2026-09-22
Status: accepted
Partially closes: the `Hu_germline (query IMGT)` requirement

## Context

Two things need a germline sequence by gene name: the fixed light chain an scFv
is built against, and the `Hu_germline` annotation on every ordered clone.

No pip package bundles germline data (`receptor_utils` and `airr` were checked
and do not). IMGT/GENE-DB and OGRDB are plain HTTP, and both returned 403 from
one machine here — the same failure mode as OPIG's OAS downloads. Availability
of a remote reference cannot be assumed.

## Decision

A registry of sources, tried in order, each result carrying its provenance:

| Source | Network | Authoritative | Covers |
|---|---|---|---|
| `fasta` | no | yes | whatever the IMGT file holds |
| `url` | yes | yes | all of IMGT/GENE-DB |
| `pool` | no | **no** | genes in the pool, spans the reads covered |

The `pool` backend takes OAS's own `germline_alignment_aa` for a gene, maps
those onto IMGT columns via `library.imgt`, and takes a per-column consensus.
That is IgBlast's IMGT-derived germline rather than a guess, and it needs
nothing external — but it is inference from observed data, and it stops where
read coverage stopped.

`Germline.is_derived` and `Germline.support` make the difference visible, and
`coverage_report()` states it per gene. A germline assembled from 12 reads and
one downloaded from IMGT must not be used interchangeably without knowing
which is which.

## Consequences

- The project works with no IMGT access, and improves when IMGT is available.
- `scripts/germlines.py --download` caches IMGT once; when the host refuses, it
  says so and points at the manual download rather than failing silently.
- Pool-derived germlines are typically 5' truncated, because the amplicon
  starts inside FR1 (see 0009). Anything needing a complete FR1 must use an
  authoritative source, and should check `is_derived` before trusting a result.
- Adding OGRDB or a local IgBLAST database later is one more registered source.
