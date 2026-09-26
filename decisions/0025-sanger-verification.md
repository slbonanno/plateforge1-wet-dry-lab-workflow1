# 0025 — Ordering Sanger, reading `.ab1` back, and deciding what is in each well

Date: 2026-09-25
Status: accepted
Extends: 0016, 0020, 0021

Real files, in `fixtures/` (rule 6):

- three Azenta/GENEWIZ traces — `fixtures/sanger/194-{M13F,M13R,3xP3}.ab1`
- the pcDNA3.1(+) map — `fixtures/backbone/pcDNA3.1.dna`
- three vendors' own order-form templates — `fixtures/sequencing/`

Every format claim below is checked against those by the test suite. This is
the first part of the repo where nothing is a guess at someone else's layout.

## The trace knows where it came from

The design assumption going in was that a rerun would have to be identified by
parsing filenames, and that ordering attempts would need file modification
times — which downloading rewrites, and often backwards, since the rerun tends
to be fetched first.

Opening a real `.ab1` made that moot. Alongside the basecalls, an ABIF file
carries the well it was run from (`TUBE`), the plate barcode and name (`CTID`,
`CTNM`), the run start date and time (`RUND1`/`RUNT1`), the instrument model,
the dye set, and the vendor's own name (`User` — literally `Azenta/GENEWIZ`).
All written by the sequencer, none of it touched by downloading.

So the well comes from the instrument, not from a filename, and attempts are
ordered by the run timestamp. The filename is a fallback and a cross-check.

## Two real reads of one molecule check the whole chain

`194-M13F` and `194-M13R` are the same insert read from both ends. Their
overlap aligns at **100% identity over 253 bp, zero mismatches, zero gaps**.

Nothing short of a correct ABIF parse, a correct quality decode, a correct
Mott trim, a correct reverse complement and a correct aligner produces that.
It is one test and it covers the lot, on real data, and it is worth more than
any number of assertions against files we generated ourselves.

## Screening on the insert, not the contig

Every well's expected sequence shares the same vector flanks. Screening a read
against all 96 expected *contigs* therefore ranks them all equally and finds
no swaps at all — the shared leader and constant region dominate.

Measured: an M13F read from an unrelated construct still shares **302 exact
12-mers with pcDNA3.1**, purely through common backbone. So the screen runs on
the **insert** alone, which is the only part where the 96 references differ.

## Verdicts, because pass/fail is not what gets asked

| verdict | meaning | usable |
|---|---|---|
| `exact` | the insert is what we ordered | yes |
| `silent` | bases differ, protein does not | yes |
| `missense` | the protein differs | only if you meant it |
| `indel` | in-frame insertion or deletion | alive, not what was ordered |
| `frameshift` | indel not a multiple of three | dead |
| `mixed` | two colonies; ambiguity codes through a good trace | re-streak |
| `no_insert` | vector all the way through | empty backbone |
| `wrong_clone` | matches a different well's insert | **a plate swap** |
| `low_quality` | never cleared the quality floor, or too little coverage | re-run |
| `unreadable` | no usable alignment | re-run |

`wrong_clone` is the one built for deliberately. It costs one k-mer screen
against the other 95 inserts, and it is the failure that silently poisons
everything downstream, because every later step faithfully carries the wrong
antibody under the right name. When it fires, the report names the clone the
well actually holds and the identity against it, so a swap is distinguishable
from a bad read at a glance.

Differences outside the insert are dropped. A mismatch in the constant region
is the vector's problem, and reporting it once per well is how a real backbone
defect gets lost in 96 copies of the same noise.

## Two alignment bugs worth recording

**Affine traceback has to be state-aware.** A single pointer per cell says
which state won there, but following it during traceback lets a gap close and
immediately reopen. The score is unaffected — the DP is correct — but a single
6 bp deletion came back as a 2 bp and a 4 bp deletion with three matched bases
between them. That is not a cosmetic difference: it turns one in-frame
deletion into something that reads like a frameshift. The traceback now
carries the current state and records, per cell, whether each gap was opened
or extended.

**A cheap gap open splits real indels.** Even with a correct traceback, at
`open = -7` the pieces of a split gap pick up enough coincidental matches
between them to pay for the extra opens. Raised to `-14`, with mismatch at
`-4`. Checked across deletions of 1, 3, 6 and 12 bp and insertions of 1 and
9 bp: each is now reported as exactly one event at the right position.

**Local at both ends, not global.** The 5′ end of a real trace is noise — the
three fixture traces start with 13 to 23 called Ns apiece — and scoring that
as gaps invents a page of indels. The alignment may begin and end anywhere.

## Reruns: solved where the data allows, flagged where it does not

Grouping is by the **sample name the submitter typed**, with recognised rerun
markers stripped (`_RR`, `-rerun`, `_redo`, `run2`, `v2`, `(2)`); ordering is
by the run timestamp. Not by plate and well — a rerun is usually cherry-picked
into a fresh plate, so it arrives with a different barcode *and* a different
well, and the sample name is the only thing that follows it. A test pins that
case.

Three deliberate limits:

- **The marker list is conservative and always will be.** `PLATE_2` is a
  plate, not a second attempt. Over-grouping silently drops a distinct sample
  as superseded; under-grouping merely asks a human. A bare `_2` suffix is
  therefore *not* treated as a rerun marker, even though it sometimes is one.
- **Nothing is deleted.** Superseded traces are returned so the caller can
  record the replacement (rule 10), not discarded.
- **The later run is not assumed to be the better one.** `prefer="later"` is
  the default because a rerun is usually done for a reason the trace cannot
  see, but when the replacement has fewer Q20 bases the report says so in
  capitals. `prefer="best"` picks by Q20 count instead.

Anything that cannot be ordered — a missing timestamp, two runs sharing one —
goes to `unresolved` for a human. **This remains the weakest part of the
module** and it is weak for a reason: without a submission manifest tying
sample names to wells, no rule can be certain, and inventing one would be
worse than asking. A vendor manifest would close it properly (Q23).

## Order forms: three vendors, three plate conventions

| | well spelling | fill order | shape |
|---|---|---|---|
| GENEWIZ / Azenta *(default)* | `A01` | both, as two columns | one table, 500 rows |
| ELIM Biopharm | `A1` | row-major | one table, 96 rows, dropdowns |
| UC Berkeley | `A1` | column-major | printed form, two side-by-side blocks |

That table is the entire argument for rule 3. A column-major plate pasted into
a row-major form is a 96-well transposition that every later step preserves
and nothing can detect.

Azenta's form is the one that fails quietly, because it offers `Well (H)` and
`Well (V)` as separate columns — two orderings of the same 96 positions,
paired row by row — and the sample numbering runs down the rows either way.
The first version of the emitter had them the wrong way round. A test now
reads the vendor's own template and asserts our pairing equals theirs for all
95 rows, which is what caught it.

Each emitter enforces the limits its own template states, before writing
rather than at upload: Azenta's 500-sample ceiling and semicolon primer
separator; ELIM's 50-character names, its `[A-Za-z0-9_-]` restriction and its
controlled vocabularies; and Berkeley's "LEAVE AT LEAST ONE WELL EMPTY ON
PLATE", which is how that facility confirms plate orientation and is therefore
a constraint, not a preference — a full 96-sample plate is refused.

## What is deliberately not done

- **No primer placement.** Where `CMV_F` lands and how much of the insert it
  reaches is not modelled. The read is aligned by sequence, which needs no
  primer at all. Worth adding when it earns its place (Q24).
- **No chromatogram inspection.** `mixed` is called from the basecaller's own
  ambiguity codes, not by re-examining peak heights in the `DATA` channels.
  The raw traces are parsed and available; reading them properly is a bigger
  job and the ambiguity codes catch the common case.
- **The dual-cassette pcDNA3.1 construct is not built.** The backbone now
  loads with its real sequence and features, which is what that work needs;
  the construct itself is the next piece.
