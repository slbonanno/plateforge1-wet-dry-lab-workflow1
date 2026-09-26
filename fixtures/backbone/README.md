# Backbone — pcDNA3.1(+), as a SnapGene `.dna` file

`pcDNA3.1.dna`, from Thermo Fisher (Invitrogen) via SnapGene. Ground truth for
`library.snapgene` (rule 6), and the real sequence behind what `library.vector`
previously described only in words.

Verified on load, every time, by `snapgene.check()`:

- **5,428 bp**, circular, ACGT only
- 13 features, including CMV enhancer (235-614), CMV promoter (615-818),
  T7 promoter (863-881), MCS (895-1002), bGH poly(A) (1028-1252),
  f1 ori, SV40 promoter and ori, NeoR/KanR, AmpR, ColE1 ori

A backbone is the one input where a silent error is unrecoverable: a one-base
slip shifts every ORF built on it, survives every test that only looks at the
insert, and surfaces months later as an expression run that produced nothing.
So it is checked against what it is supposed to be before anything uses it,
rather than trusted.

## Format

A flat run of segments: one type byte, a big-endian length, then the payload.

    0   the DNA: one flags byte (bit 0 = circular), then the bases
    3   enzyme recognition list (SnapGene's own UI state; skipped)
    5   primers, as XML
    6   notes, as XML
    8   additional sequence properties, as XML
    9   the header: b"SnapGene" plus three 16-bit version fields
    10  features, as XML

Feature coordinates in the XML are 1-based and inclusive. That is converted
once, on read, so nothing downstream has to remember it.

## Next

The construct this repo actually wants is **CMV-VH and CAG-VL as two
cassettes** (decision 0020) built on this backbone. The backbone is now
loaded with real coordinates, which is what that work needed; the construct
itself is not built yet.
