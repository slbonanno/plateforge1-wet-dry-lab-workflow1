# 0020 — Two promoters, not one ORF

Date: 2026-09-23
Status: accepted
Supersedes: 0016's choice of default vector (the rest of 0016 stands)

## Context

0016 made the default a full IgG1 from a single ORF, heavy and light chain
separated by furin–GSG–T2A. It expresses, and it is one transcript, which is
tidy. It is also the wrong ratio.

2A skipping is incomplete, and the fraction that does skip translates the
downstream ORF at a lower level than the upstream one. IgG assembly wants the
opposite: light chain in **excess** of heavy, conventionally 1.5–2:1. Free
heavy chain that cannot find a partner aggregates and takes assembled yield
down with it. HC-first T2A delivers less light chain than heavy and offers no
handle at all — the ratio is a property of the peptide, not of anything a
designer can change.

## Decision

The default is `pcdna-igg1-dual-vk-v1`: two cassettes, each with its own
promoter and terminator.

```
HC (CMV):  CMV — T7 — Kozak — ATG+leader — [VH] — CH1-hinge-CH2-CH3 — stop — BGH pA
LC (CAG):  CAG — Kozak — ATG+leader — VL — CK — stop — SV40 pA
```

CAG on the light chain is deliberate: it is the stronger promoter, so the
excess falls on the right side. The ratio is now a design parameter — change a
promoter and you change it.

Two smaller choices, both to avoid trouble later: each cassette has its own
start codon, because they are separate transcripts rather than a continued
reading frame; and the two terminators are different, because a repeated
several-hundred-bp element in one plasmid invites recombination.

### `Cassette` is now a first-class thing

`Vector` holds `cassettes`, each a promoter → reading frame → terminator, with
at most one carrying the insert. `upstream` and `downstream` still exist and
still mean what they meant — what precedes and follows the insert *within its
own cassette*. Anything in another cassette is neither, which is the point:
with two promoters the light chain is not downstream of the heavy chain in any
sense a reading frame cares about.

`cloning.verify` now translates **every** cassette. The second cassette cannot
be broken by the insert, but it can be broken by an edit to the vector, and
that is exactly the failure nobody would notice.

`Vector.declared_bp` reports the length of elements whose sequence is actually
declared and counts the ones that are not, rather than returning a plasmid
size that silently omits its own promoters.

## What does not change

**The ordered fragment.** The VH still sits between the heavy chain's leader
and CH1, so the fusion sites are identical (`AGGA` / `GCCA`) and a fragment
designed against the T2A vector is byte-for-byte the fragment designed against
this one. A test asserts it. Swapping vector does not reprint the order.

## Size

~6.8 kb with bacterial selection only, ~8.4 kb with a NeoR cassette for
stables. Only the constant regions and the VH are exact; promoter and backbone
lengths are typical values, and `declared_bp` will not pretend otherwise.

| | bp |
|---|---|
| HC cassette | 2,379–2,412 |
| LC cassette | 2,546 |
| Backbone (ori + AmpR + spacers) | 1,839 |

## The T2A vector stays registered

It is a legitimate format — fewer parts, smaller plasmid, and fine when
absolute yield does not matter. It is not the default, and this record is why.
