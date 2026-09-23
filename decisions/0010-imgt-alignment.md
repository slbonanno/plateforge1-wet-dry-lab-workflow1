# 0010 — Align on ANARCI IMGT numbering, not string position

Date: 2026-09-22
Status: accepted
Closes: Q12

## Context

The alignment figure assumed `sequence_alignment_aa` was a fixed-width IMGT
alignment, so sequences of one germline would share columns. That is false.
Measured on Briney et al. 2019:

| V gene | n | distinct lengths | min | max |
|---|---|---|---|---|
| IGHV3-23 | 16,638 | 44 | 62 | 124 |
| IGHV1-69 | 2,546 | 31 | 70 | 121 |
| IGHV3-53 | 679 | 27 | 67 | 123 |

Reads cover different spans of the domain. Laying them out by string position
puts residues in the wrong columns; the interim guard restricted to the modal
length, which is correct but discards most of the data.

## Decision

Align on `ANARCI_numbering`, which OAS carries per sequence: a per-residue IMGT
position label, grouped by region.

```
{'fwh1': {'15 ': 'P', '16 ': 'G', ...},
 'cdrh3': {'111 ': 'G', '111A': 'I', '112A': 'D', '112 ': 'R', ...},
 'fwh4': {'118 ': 'W', ...}}
```

Each residue already knows its column, so sequences of any length drop into a
shared frame with no aligner. Regions come from the same structure, which also
supplies FR4 — the per-region `fwr*_aa` columns do not.

## Insertion ordering

IMGT inserts into CDR3 outward from the middle: `111, 111A, 111B, … 112C,
112B, 112A, 112`. Insertions after 111 ascend and insertions before 112
descend, so a plain string sort is wrong. `imgt.position_key` implements this.

## The germline reference

`germline_alignment_aa` is aligned to `sequence_alignment_aa` position by
position, and the numbering covers that sequence's non-gap residues in order,
so zipping the two places the germline into the same columns.

When the counts disagree the row is skipped rather than guessed at. A germline
row one residue out of register would make every position look mutated — a
confidently wrong figure is worse than a missing one.

## Consequences

- Ragged real reads render correctly instead of being excluded.
- `min_occupancy` trims insertion columns only a few sequences carry.
- The string-position path remains as a fallback for data without numbering,
  and still restricts to the modal length and says so.
- Nothing else depends on numbering yet, but it is now in the pool, and it is
  the natural basis for position-wise SHM analysis and CDR3 logos.
