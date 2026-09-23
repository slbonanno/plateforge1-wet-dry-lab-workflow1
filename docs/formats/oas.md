# Observed Antibody Space (OAS)

**Direction:** in
**Primary route:** the HuggingFace parquet mirror (`ConvergeBio/oas-unpaired`)
**Fixture:** none yet (Q9) — synthetic units from `library.synth` stand in
**Verified:** mirror schema verified against a real shard, 2026-09-22.
OPIG's own download paths verified **broken** the same day.

## OPIG direct access does not currently work

As of 2026-09-22, `opig.stats.ox.ac.uk` returns **403 for the data files
themselves** while serving directory listings normally:

| Path | Result |
|---|---|
| `/webapps/ngsdb/paired/` | 200, browsable listing |
| `/webapps/ngsdb/paired/Jaffe_2022/csv/` | 200, filenames visible |
| `/webapps/ngsdb/paired/Jaffe_2022/csv/*.csv.gz` | **403** |
| `/webapps/ngsdb/unpaired/` | **403** |
| `/webapps/oas/downloads/` via the search form | **403** (Apache, `naga-www-prod`) |

Confirmed from two unrelated networks and with a browser User-Agent, so it is
not an IP block, a client problem, or a missing login — OAS has no accounts and
the data is CC-BY 4.0. Either an anti-scraping ACL or a misconfiguration on
OPIG's side. **OAS has no API**, so when the files are blocked there is no
other first-party route.

`oas.load()` still reads a `.csv.gz` from a URL or a local path, so this comes
back the moment OPIG unblocks it, or if you obtain a unit some other way.
Contact: `oas_opig@stats.ox.ac.uk`.

## The mirror

`ConvergeBio/oas-unpaired` on HuggingFace: OAS unpaired as parquet, validated
against the original files. 2.07B heavy and 357M light sequences, 4,030 shards,
2.45 TB total — but parquet over range requests means a column projection plus
a row-group limit reads megabytes.

```
data/unpaired_heavy/<Study>/train-000NN.parquet
data/unpaired_light/<Study>/train-000NN.parquet
```

Studies are spelled as citations: `Briney et al., 2019`. `hf.study_files()`
matches case- and punctuation-insensitively so you need not reproduce that
exactly.

### Two differences from a raw OAS data unit

**Metadata is columns, not a header.** A `.csv.gz` unit carries a JSON blob on
line 0; the mirror carries `meta_Run`, `meta_Author`, `meta_Species`,
`meta_Age`, `meta_BSource`, `meta_BType`, `meta_Vaccine`, `meta_Disease`,
`meta_Subject`, `meta_Longitudinal`, `meta_Isotype`, `meta_Chain`, `meta_Link`.
`hf.extract_meta()` folds these back into the dict `normalize()` expects.

**One folder holds every species.** `unpaired_heavy` mixes human, rabbit, mouse
and rhesus. The first row of `Banerjee et al., 2017` is `IGHV1S45*01` — rabbit.
A species filter is therefore not optional, and `hf.load_study()` defaults to
human for that reason.

### Columns confirmed present in a real shard

AIRR standard, plus OAS's own. Everything the pipeline needs survived the
mirroring:

| Column | Note |
|---|---|
| `sequence_alignment_aa` | the amino acid sequence |
| `germline_alignment_aa` | IgBlast's germline, aligned to it — the alignment figure's reference row |
| `v_germline_alignment_aa` | V-region germline only |
| `v_call`, `j_call`, `d_call` | IMGT allele calls, **may be comma-separated ties** |
| `fwr1_aa` … `cdr3_aa` | per-region strings, used for the IMGT region bands |
| `fwr1_start`/`_end` … `fwr4_start`/`_end` | region coordinates, including FR4 |
| `Redundancy` | count after sequence collapse |
| `ANARCI_status` | pipe-delimited flags (see below) |
| `ANARCI_numbering` | per-residue IMGT numbering; the fallback if sequences are ever not column-aligned |
| `aa_hash_hi`, `aa_hash_lo` | added by the mirror, not in OAS |

## Gotchas, all found in real data

**Most `ANARCI_status` flags describe the read, not the molecule.** Measured
on Briney et al. 2019, 19,864 sequences surviving gene and length filters:

| Flag | Count |
|---|---|
| `Shorter than IMGT defined: fw1` | 17,592 |
| `Deletions: 1..15, 73` | 11,471 |
| `Deletions: 1..14, 73` | 3,453 |
| `Shorter than IMGT defined: fw1, fw4` | 2,272 |
| `Missing Conserved Cysteine: 104` | 480 |

Treating every flag as a liability discarded **all 125,743** rows — which is
how this was found. "Shorter than IMGT defined" and low-numbered "Deletions"
are the same artifact: the amplicon starts inside FR1, so the first dozen IMGT
positions are not covered. Rejecting those means rejecting a repertoire on the
basis of primer design.

`oas.COVERAGE_PATTERNS` classifies those as coverage; everything else,
including anything unrecognised, stays a liability. A judgment call worth
knowing: "Deletions" is treated as coverage although a mid-domain deletion
would be real. Position 73 is deleted in ~90% of Briney sequences, which is
systematic rather than 18,000 real events. Revisit if a dataset shows sparse
isolated deletions.

**Ambiguous residues are common and must be filtered.** Real CDR3s come back
looking like `XXXXXXXXXXYW`. `normalize()` counts non-standard residues into
`n_ambiguous` and `Filter(max_ambiguous=0)` drops them by default. Without it
they reach the sampler and end up in a picked library.

**V/J calls are often ties.** `"IGHV3-23*01,IGHV3-23*04"` used as a gene name
produces a silent garbage category. All parsing goes through
`library.germlines`.

**`sequence_alignment_aa` is NOT a fixed-width IMGT alignment.** This was
assumed and is false. Measured on Briney et al. 2019:

| V gene | n | distinct lengths | min | max | column-aligned |
|---|---|---|---|---|---|
| IGHV3-23 | 16,497 | 54 | 62 | 124 | no |
| IGHV1-69 | 2,708 | 38 | 70 | 121 | no |
| IGHV3-53 | 659 | 32 | 67 | 123 | no |

Sequences of one germline do **not** share columns, because reads cover
different spans of the domain. `hf.alignment_report()` measures this, and
`figures.alignment()` restricts to the modal length and states how many it
excluded rather than padding them into the wrong register — correct, but it
discards most of the data.

`ANARCI_numbering` solves this and is now used (decisions/0010). It carries a
per-residue IMGT position per region:

```
{'fwh1': {'15 ': 'P', '16 ': 'G', ...},
 'cdrh3': {'111 ': 'G', '111A': 'I', '112A': 'D', '112 ': 'R', ...},
 'fwh4': {'118 ': 'W', ...}}
```

Every residue knows its column, so ragged reads align with no aligner, and the
regions come from the same structure — including FR4, which the `fwr*_aa`
columns do not provide. Note the insertion ordering: IMGT inserts into CDR3
outward from the middle, `111, 111A, 111B, … 112C, 112B, 112A, 112`, so a
string sort of position labels is wrong.

**Safari corrupts `.csv.gz` downloads** by auto-unzipping. OPIG warns about
this; use another browser or turn the setting off.

## Still to verify

Against real data, once a pool exists: whether `germline_alignment_aa` is ever
empty; the full `ANARCI_status` flag vocabulary and which flags deserve to be
treated as liabilities; whether region columns ever fail to tile the sequence;
and whether the paired mirror (`unpaired_light` is separate — paired units are
a different layout entirely) needs its own adapter.
