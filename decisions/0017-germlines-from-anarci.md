# 0017 — Germline references come from the tables bundled with anarci

Date: 2026-09-23
Status: accepted
Closes: Q17

## Context

Decision 0012 built a source registry for germline references because there was
"no offline germline database in any pip package" and IMGT/GENE-DB was a plain
HTTP download that may or may not be reachable. That turned out to be wrong in
a useful way: `anarci` is on PyPI, and it ships the IMGT germline tables as
Python data.

Measured from this environment, IMGT itself is unreachable — GENE-DB returns
`403 Forbidden` through the proxy, and its web interface is disallowed by
robots.txt. So the download route was never going to work here, and Q17 ("is
IMGT reachable from the working machine?") has been open since 2026-09-22
precisely because nobody could answer it.

## Decision

`germline_db.from_anarci` reads `anarci.germlines.all_germlines` and is first
in the default resolve order: `("anarci", "fasta", "pool")`.

It is authoritative, complete for human V and J genes, and needs no network.
`pip install -e ".[imgt]"` pulls it in; without it, everything falls back to
the previous behaviour and says so.

## What this changes downstream

The alignment reference stops being "OAS germline consensus, read span" and
becomes the actual IMGT reference. Two consequences worth stating:

1. An IMGT V gene stops at the conserved `...YYC` plus two residues. On its own
   it gives CDR3 and FR4 no germline at all, so every sequence's C-terminus
   would render as one long insertion. `msa.reference_for` therefore appends
   the **modal J gene these reads call**, and the figure names both alleles.
2. The D/N region between V and J stays an insertion, because it genuinely is
   not germline-encoded. That is the correct picture of a CDR3, not a gap in
   the reference.

## Consequences

- Figures now read `IMGT IGHV3-23*01 + IGHJ4*01` rather than a consensus label.
- `germline_db.known_genes()` lists everything the tables carry (55 human IGHV).
- The light chain in `library.vector` is built from these tables rather than
  typed in from memory.
- `scripts/germlines.py --download` is no longer on the critical path. It stays,
  because a GENE-DB FASTA on disk is still a valid and useful source.
