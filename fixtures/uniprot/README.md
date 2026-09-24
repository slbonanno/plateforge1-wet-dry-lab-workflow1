# UniProt records, as retrieved

Rule 6 applies to reference sequences as much as to file formats: a constant
region typed from memory is a guess, and a guess in a 96-well synthesis order
is expensive. These are the unmodified FASTA records the constant-region
constants in `library/vector.py` were derived from.

| File | Accession | Retrieved | From |
|---|---|---|---|
| `P01857_IGHG1_HUMAN.fasta` | P01857 (SV=2) | 2026-09-23 | `https://rest.uniprot.org/uniprotkb/P01857.fasta` |
| `P01834_IGKC_HUMAN.fasta` | P01834 (SV=2) | 2026-09-23 | `https://rest.uniprot.org/uniprotkb/P01834.fasta` |

## One derivation, and why

P01857 is the **membrane-bound** isoform. Its C-terminus reads
`...HYTQKSLSLSP` followed by `ELQLEESCAEAQDGELDGLWTTITIFITLFLLSVCYSATVTFFKVKWIFSSVVDLKQTIIPDYRNMIGQGA`,
which is the M1/M2 transmembrane tail.

A secreted IgG1 ends `...HYTQKSLSLSPGK`. So `IGHG1_CH1_CH3_AA` is the UniProt
sequence truncated at `KSLSLSP` with `GK` appended. `library.vector` performs
that truncation from the stored record rather than hardcoding the result, and a
test asserts the two agree — so if the record is ever replaced the constant
cannot silently disagree with it.

P01834 is used unmodified.
