# 0018 — Gap penalties follow IMGT regions

Date: 2026-09-23
Status: accepted
Extends: 0014

## Context

0014 aligns each sequence to its germline with one gap-open penalty for the
whole molecule. That is wrong in a specific, predictable way.

An antibody V region is not homogeneous. FR1–FR3 are germline-encoded and
near-invariant; a gap in one is almost always the aligner making a mistake
rather than biology. CDR3 is the product of V(D)J recombination, exonuclease
chew-back and N-addition — it has **no germline to be indel-free against**. The
IMGT reference makes this literal: a V gene stops at position 106, a J gene
starts at 115, and positions 107–114 have no germline residue at all.

With one penalty, opening the junction costs framework prices (−11), so the
aligner buys mismatches instead — packing junction residues into germline
columns and scoring them as framework substitutions.

## Decision

`msa.region_gap_profile` produces a (open, extend) pair per reference position
from the reference's IMGT numbering, which the `anarci` tables supply directly
(their strings are IMGT-gapped to 128, so index + 1 is the number).

| Region | IMGT | open | extend |
|---|---|---|---|
| FR1, FR2, FR3, FR4 | 1–26, 39–55, 66–104, 118–128 | −11 | −1 |
| CDR1, CDR2 | 27–38, 56–65 | −6 | −0.8 |
| CDR3 | 105–117 | −1.5 | −0.2 |

A seam takes the more permissive of the two regions it sits between, so the
V/J junction is relaxed from its first column rather than one residue late.

`align_pair` gained a `gap_profile` argument; `build(region_aware=True)` is the
default and is silently inert when the reference carries no numbering.

## What it changes, measured

On 360 sequences across IGHV3-23, IGHV1-69 and IGHV3-53:

| | mean identity to germline | total alignment width |
|---|---|---|
| uniform | 0.9512 | 480 |
| region-aware | **0.9557** | 501 |

Identity rises and the alignment gets *wider*. Both are the intended effect:
junction residues get their own insertion columns instead of being forced into
germline columns and counted as framework substitutions. The extra columns are
CDR3 columns, which is where the diversity actually is.

On unmutated flanks the two agree — the difference appears exactly where real
data lives, in somatically mutated sequences with long junctions.

## What it does not change

Framework penalties are untouched, so a genuine framework deletion is scored
exactly as before. This loosens the junction; it does not loosen the alignment.

## Alternative considered

A profile HMM with per-position emission and transition probabilities (HMMER,
or an antibody-specific model) is the principled version of this and would
also handle CDR1/CDR2 length classes. It needs training data and a dependency.
The region table is three lines of definition and captures most of the effect;
revisit if the alignments are ever load-bearing for anything but figures.


## Addendum, same day — the junction was still being torn apart

The region profile fixed the *cost* of opening the junction but not what
happened to it afterwards. Two things were still wrong, and together they made
CDR3 render as ragged expansion rather than alignment.

**Germline CDR3 residues were splitting the block.** A V gene contributes the
A-R of C-A-R and a J gene its F-D-Y; ANARCI's tables number them 105-106 and
115-117. Left in the reference they are anchors, and the aligner happily
matched arbitrary junction residues against them — chopping one CDR3 into
three or four separate insertion blocks, each laid out independently. The
reference now drops every residue numbered 105–117, leaving exactly one seam
between the conserved cysteine (104) and tryptophan (118). Nothing real is
lost: those residues *are* CDR3, and each sequence shows its own.

**Insertion blocks were left-justified.** `before[i].ljust(width)` padded on
the right and never aligned segments to each other, so a 10-mer and a 25-mer
junction shared their first column and nothing else. Insertions inside a CDR
are now placed from both ends with the gap in the middle (`msa.place`), which
is what IMGT itself does — its CDR3 insertion order runs 111, 111A, 111B …
112B, 112A, 112, accumulating outward and meeting in the middle. A CDR3 has
two anchors, not one.

**Deletions were relaxed along with insertions.** They should not be: "this
read is missing framework" is a strong claim at any position. The profile now
governs insertions only.

Measured on the same pools, IGHV3-23 at 50 sequences:

| | alignment width | CDR3 columns |
|---|---|---|
| before | 157 | 41 |
| after | 126 | 19 |

19 columns for junctions of 8–19 residues is the floor — you cannot show a
19-mer in fewer. About a third of that block is gaps, which is simply what it
costs to put an 8-mer and a 19-mer in shared columns, and is now visibly
length variation rather than a rendering artefact. `VYYCAR` lines up on the
left and `DYWGQGT` on the right, in every row.


## Addendum 2 — the two outliers

Even with one contiguous block, a real run showed two sequences out of ninety-
six with their whole junction one column left, and one right, of the stack.

The cause is an unoccupied anchor column. When a read's conserved cysteine
(IMGT 104) is missing or mutated, the pairwise aligner has a free column
adjacent to the junction and puts the junction's first residue in it — that
sequence's junction then starts one column left of everyone else's and is one
residue short, so it gains an interior gap too. The tryptophan (118) does the
same thing in mirror image.

**Scoring cannot fix this**, and the obvious attempt makes it worse. Forbidding
any non-anchor residue at the anchor columns — matching by residue identity —
lets a junction that happens to contain its own W capture the tryptophan
column from across the block, tearing the junction in two. Measured: it turned
a clean 12-residue junction into `ARGL--...--REYWKIDY` split across the
anchor.

So it is fixed structurally, after alignment and before layout.
`msa.consolidate_junction` enforces the definition rather than inferring it:
an anchor column holds its own residue or nothing, and everything between the
two anchors is junction. Anything sitting in an anchor column that is not that
anchor is moved into the junction block, in order.

That makes the invariant unconditional — every junction starts in the first
CDR3 column and ends in the last — regardless of what the pairwise aligner
did, and regardless of whether the read has its anchors at all.
