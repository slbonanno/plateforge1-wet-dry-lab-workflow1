# 0014 — Align to the germline with gaps; numbering belongs to the reference

Date: 2026-09-23
Status: accepted
Supersedes: 0010 (for the figure; see "What 0010 keeps" below)

## Context

0010 replaced string-position layout with ANARCI IMGT numbering, which was the
right fix for the wrong problem. Numbering *classifies* residues — it says
which structural position each one occupies. That is what you want when
comparing across germlines. It is not what an alignment is.

It also carried three costs:

- A sequence without `ANARCI_numbering` could not appear at all. Rows written
  before the column existed rendered as blank axes (see 0011), and any source
  that does not ship ANARCI output is simply unusable.
- Nothing in the picture was a gap. Positions a sequence did not cover and
  positions it genuinely lacked looked the same, because both were "no residue
  at this IMGT column".
- The germline row was assembled per column from whatever OAS reported, which
  is a reasonable germline but not *a sequence* — it could not be exported,
  handed to another tool, or checked.

And the report it drove (`column_aligned`) measured whether
`sequence_alignment_aa` is fixed width, which is never true for real reads and
never was the thing the figure needed. That label was wrong from v7 onward.

## Decision

Align each V gene's sequences **to that gene's germline**, pairwise, with gaps,
and merge the pairwise results into one column space anchored on the germline.
`library.msa` does this: Needleman–Wunsch with affine gaps and BLOSUM62, end
gaps free on the query so a 5′ truncated read keeps register instead of paying
for deletions it never had.

The column space is built from the germline outward:

- Every germline residue owns exactly one column, in order.
- Where a sequence carries residues the germline does not, an insertion column
  opens and every other row shows a gap there.
- **Only the germline row is numbered**, by its own residue index. An
  insertion column is not a germline position, so it has no number. A
  variant's residues are never numbered, because numbering them would assert a
  position the germline does not define.

The germline is a real sequence, resolved by `msa.reference_for`: an
IMGT/GENE-DB FASTA when one is on disk, otherwise a consensus of the germline
calls OAS ships with these very reads, labelled as such on the figure so no
picture can imply a reference it did not have.

Every alignment is written as PNG **and** FASTA **and** plain text, so it can
be opened in Geneious or Jalview and the picture can be checked.

## Why not Clustal or MAFFT

Both are fine and both are supported (`backend="mafft"`, `backend="clustalo"`,
used if the binary is on PATH). Neither is the default, for three reasons:

1. They are an install. This repo runs from a clean `pip install -e .` and
   should keep doing so.
2. A tree-guided MSA optimises a different objective. For a set of sequences
   that are all one germline's children, a reference-anchored alignment and a
   progressive MSA agree; where they disagree, the reference-anchored one is
   the one whose columns mean "germline position N", which is the whole point.
3. A general aligner puts the germline wherever the guide tree does. Here it
   is in column space by construction, which is what makes the numbering safe.

## Colour means one thing

A fill means *this differs from the germline*. So the germline row is never
filled — it is grey letters on white, the baseline the eye reads everything
else against. Filling it, as the first version did, made the reference the
loudest thing in the picture and left the reader decoding a legend to find the
mutations.

An inserted residue is coloured even though there is nothing above it to differ
from: the germline has no residue at that column, which is itself the
divergence. A gap where the germline does have a residue is tinted, not
coloured — there is no residue to colour.

The summary panel uses the same rule, so one reading applies to every alignment
this repo produces.

## What 0010 keeps

IMGT numbering is still correct and still used — for *classification*, not
layout: `library.imgt` maps residues to IMGT positions, and
`germline_db.from_pool` uses it to build a per-position germline consensus
across reads of differing span. What changed is that the alignment figure no
longer depends on it, so a pool with no ANARCI output still aligns.

## Consequences

- No sequence is excluded from an alignment for being short or unnumbered.
- `figures.alignment` takes `backend=` and `fasta_path=` and writes sidecar
  files. `figures.draw_alignment` and `_alignment_from_numbering` are gone;
  `figures.build_msa` / `figures.draw_msa` replace them.
- `hf.alignment_report`'s `column_aligned` no longer gates anything.
- The summary panel shows 50 sequences by default, not 10.
