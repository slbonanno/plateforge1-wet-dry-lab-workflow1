# Sanger traces — three real `.ab1` files from Azenta/GENEWIZ

These are ground truth for the ABIF format. `library.abif` was written against
them (rule 6), and the tests read them rather than anything we generated.

| file | sample | well | mean Q | Q20 bases | what it is |
|---|---|---|---|---|---|
| `194-M13F.ab1` | 194-M13F | F5 | 42.4 | 828 / 929 | a clean forward read |
| `194-M13R.ab1` | 194-M13R | G5 | 32.1 | 719 / 959 | the same insert, reverse |
| `194-3xP3.ab1` | 194-3xP3 | H5 | 11.8 | 54 / 832 | **a failed read** |

They are from an unrelated project and vector; nothing about the sequence
corresponds to this repo's construct. That does not matter, and the set is
better than a matched one would have been:

**The forward/reverse pair is the strongest check in the repo.** Both read the
same molecule from opposite ends, so their overlap must be a perfect reverse
complement. It is: 253 bp, 100% identity, zero mismatches, zero gaps. Nothing
short of a correct parse, quality decode, Mott trim, reverse complement and
aligner produces that, so one test covers the whole chain on real data.

**The failed read is a fixture you cannot fabricate honestly.** It is what a
real dead well looks like — 832 called bases, almost none of them meaning
anything — and it pins the behaviour that matters: trimming returns an empty
window and the well is reported `low_quality` rather than having a sequence
read into it.

## What these files taught us

Opening them changed the design. An ABIF file carries, written by the
sequencer and untouched by downloading:

    TUBE          the well            F5
    CTID / CTNM   plate barcode/name  04A000129809 / 2825298
    RUND1/RUNT1   run start           2024-02-08 00:57:21
    MODL          instrument          3730
    DySN          dye set             Z-BigDyeV3
    User          the vendor          Azenta/GENEWIZ
    PBAS / PCON   bases and Phred scores

So the well never has to be parsed out of a filename, and reruns are ordered
by the instrument's own timestamp rather than by file modification time —
which downloading rewrites, usually backwards. See decision 0025.

`PCON` is a char array whose *bytes* are Phred scores. Decoding it as text
gives every base a plausible-looking quality in the 30s-60s and is meaningless.

## Still worth capturing

- **A genuine rerun pair** — the same sample submitted twice, as the vendor
  named them. Rerun grouping is the weakest part of the module (Q23) and one
  real example would say whether the marker list matches reality.
- **A vendor submission manifest** — the spreadsheet tying sample names to
  wells. It would replace the heuristic entirely.
- **A mixed trace** — two colonies in one well. `mixed` is called from the
  basecaller's ambiguity codes and has never been tested against a real one.
