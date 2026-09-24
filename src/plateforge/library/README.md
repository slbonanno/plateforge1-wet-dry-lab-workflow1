# library

Sequence acquisition and in-silico construct generation.

**Inputs:** OAS data units (URL or local file), a germline panel, a sample size
**Outputs:** `SEQ` records in the pool, `RUN` ingest artifacts, `LIB` sequence sets

## What exists

| File | Does |
|---|---|
| `germlines.py` | IMGT call parsing (ties, alleles, families) and the scaffold panels |
| `oas.py` | Reads OAS data units: metadata line, streaming, normalization, filtering |
| `pool.py` | Persists the pool — parquet for full rows, SQLite for the filter index |
| `diversity.py` | Diversity-aware sampling and the acceptance criteria |
| `synth.py` | Synthetic OAS-format units for testing and dummy-data generation |
| `hf.py` | Reads the HuggingFace parquet mirror; maps its `meta_*` columns back to an OAS metadata dict |
| `imgt.py` | Maps residues to IMGT positions — classification, not layout (see `msa.py`) |
| `msa.py` | Aligns sequences to their germline with gaps; the germline owns the columns and the numbering |
| `selection.py` | Writes one selection's run bundle: what was fetched, how it was filtered, what was chosen, one alignment per family |
| `germline_db.py` | Germline references by gene name, from IMGT or derived from the pool |
| `figures.py` | Regenerable figures; figure 3 is a visual acceptance test |
| `construct.py` | scFv and VH-only assembly, clone minting, plate layout |
| `vector.py` | What the ordered fragment lands in: the elements the vector already carries, and where the ligation boundaries fall |
| `cloning.py` | How the fragment joins: Golden Gate, Gibson, blunt — one registered strategy each |
| `vendors.py` | Order tables, one registered emitter per vendor, each declaring whether its layout has been verified against a real template |
| `ordering.py` | Composes the four: molecule, vector, strategy, vendor |
| `codon.py` | Reverse translation, GC balancing, BsaI and homopolymer removal |

Two entry points:

```bash
python scripts/demo.py                                   # synthetic, offline
python scripts/fetch_oas.py --list-studies               # what the mirror holds
python scripts/fetch_oas.py --study "Briney et al., 2019"  # real sequences
python scripts/fetch_oas.py --source ~/Downloads/unit.csv.gz   # a local unit
```

Real data comes from the HuggingFace parquet mirror of OAS, because OPIG's own
download paths return 403 and OAS has no API. Install the extra for it:

```bash
pip install -e ".[hf]"
```

`docs/formats/oas.md` records what is blocked, what the mirror's schema looks
like, and the gotchas found in real data.

Generated figures go under `figures/`, which is gitignored and split so a
synthetic figure is never mistaken for a real one:

```
figures/synthetic/   from demo.py
figures/real/        from fetch_oas.py
```

Nothing under `figures/` is committed. Regenerate rather than checking one in —
a stale figure that no longer matches the code is worse than no figure.

## The shape of a run

```python
from plateforge.library import diversity, germlines, oas, pool

# 1. Read one data unit, straight from OAS or from a file already downloaded.
df, meta = oas.load(
    oas.unit_url("Chen_2020", "SRR11937587_1_Heavy_IGHG"),
    oas.Filter(genes=germlines.DEFAULT.genes, cdr3_len_range=(8, 30)),
)

# 2. Add it to the pool. Idempotent — seq_id is a content hash.
run_id = pool.ingest(df, meta)

# 3. Sample a campaign-like set.
picked = diversity.sample(pool.index(), 10, panel=germlines.DEFAULT, seed=42)

# 4. Check it, then freeze it as a LIB artifact.
report = diversity.Spec(min_genes=3).check(picked)
lib_id = pool.make_library("campaign-a", list(picked["seq_id"]),
                           produced_by="library.diversity",
                           parents={run_id: "sampled_from"})
```

Downstream modules receive `lib_id`, never objects from this module.

## Figures

`figures.all_figures()` writes four; `alignment()` and `aa_legend()` are separate
because the alignment needs one germline at a time.

| Figure | Shows |
|---|---|
| 1 germline composition | pool vs sample vs panel target — did the sampler hit quota |
| 2 CDRH3 lengths | per-germline length distributions, sampled members marked |
| 3 identity check | sampler vs random draws; the visual acceptance test |
| 4 V gene usage | what is actually in the pool, grouped by family |
| 5 alignment | per-germline, every residue printed, only divergence from germline coloured, IMGT regions banded above |
| input_summary | germline usage, CDRH3 lengths and an example alignment in one wide figure, for a README |

The alignment view works without running an aligner because OAS ships
`sequence_alignment_aa` IMGT-gapped, so sequences from one germline are already
column-aligned. `normalize()` keeps that as `aa_gapped` alongside the ungapped
`aa_seq`, and keeps OAS's own `germline_alignment_aa` as `germline_aa`. Both
live in parquet only, so fetch them with `pool.fetch(...)` rather than
`pool.index()`.

Gaps are drawn as explicit dashes rather than blank cells, so a gap is
unmistakably a gap. The x-axis numbers the **germline's own residues**: a gap
column is a position some other sequence carries and the germline does not, so
it consumes no germline number. Sequences below the germline are unnumbered.

The reference row is **the germline**, taken from the `germline_alignment_aa`
column OAS ships per sequence — IgBlast's own call, already column-aligned, so
there is no germline table to maintain and nothing asserted that the data did
not provide. `reference_row()` falls back to the observed consensus only when a
data unit carries no germline column, and labels the row so a figure can never
silently imply a germline it did not have.

Every residue letter is printed. Positions matching the germline appear in
muted ink with no fill; only substitutions get a coloured cell, so the sequence
stays readable and divergence is what catches the eye. FR and CDR bands are
drawn above from OAS's own per-region columns (`fwr1_aa`, `cdr1_aa`, ...),
located by searching for each region in the sequence rather than trusting
cumulative lengths, so overlapping or non-tiling regions degrade to "no band"
instead of a wrong one.

Residue colouring is one hue per amino acid, banded by property so families
read as colour neighbourhoods. Twenty categories cannot be separated by colour
alone — that is perception, not palette design — so hues were chosen by
maximising the minimum pairwise OKLab distance under normal, deuteranopic and
protanopic vision within each band (worst pair ≈5.0 ΔE), and the residue letter
is drawn in every cell with ink picked per fill luminance.

## From picked sequences to an order

```bash
python scripts/order.py --pick 94 --reserve A01 H12
```

Writes to `outputs/<clone-set-id>/`:

| File | Contents |
|---|---|
| `clones.csv` | one row per clone: ids, construct, protein, DNA, GC, warnings, well |
| `plate_map.csv` | the 96-well layout as a table |
| `plate_map.txt` | the same as a grid, CDRH3 per well, for reading at the bench |
| `alignment_<gene>.png` | the clones that will be synthesised |

This is where a `SEQ` becomes a `CLN` (decisions/0005). DNA is reverse
translated with human codon usage, balanced to mid-GC, and stripped of the
assembly enzyme's site on both strands and of homopolymer runs — see
decisions/0013 for what real output forced.

Vendor order forms are deliberately absent: those need a real template in
`fixtures/` first (Q15), and Gibson/HiFi adapters need the destination vector
(Q16).

## Germline references

```bash
python scripts/germlines.py --download        # cache IMGT/GENE-DB once
python scripts/germlines.py --panel --show    # resolve the scaffold panel
```

Sources are tried in order and every result says where it came from. `fasta`
and `url` are IMGT and authoritative; `pool` is a per-column consensus of the
germline OAS reports for that gene — no network, but it covers only genes in
the pool and only the span the reads covered, which for 5' truncated amplicons
is not all of FR1. Check `is_derived` before treating one as a reference.

## Reproducibility

Synthetic generation is seeded and must produce identical data in any process.
That is not automatic: Python randomises string hashing per process, so
iterating a `set` while consuming the RNG makes the same seed yield different
data on different machines. `test_generation_is_reproducible_across_processes`
guards it by running generation under several `PYTHONHASHSEED` values and
comparing digests.

## What a run produces

Selecting sequences writes a bundle under `outputs/<LIB id>/` (decision 0015):
`manifest.json` (source, filter funnel, sampler seed, diversity report, code
version), `selected.csv`, `pool_summary.csv`, and one
`alignment_<gene>.png` + `.fasta` + `.aln.txt` per V gene. A re-run mints a new
`LIB` id and gets its own directory, so no run overwrites another's evidence.

`python -c "from plateforge.core import runs; print(runs.describe('LIB-...'))"`
summarises a directory you have found and cannot place.

## Alignments

Sequences are aligned **to their germline**, with gaps (decision 0014). The
germline is the first row and the only numbered one — insertion columns are not
germline positions, and a variant's residues are never numbered. Nothing is
excluded for being short: a 5' truncated read aligns with leading gaps.

Gap penalties follow IMGT regions (decision 0018): frameworks keep the standard
penalty, CDR1/CDR2 are loosened, and CDR3 is very cheap to open — it is the
product of V(D)J recombination and has no germline to be indel-free against.
Without this the aligner pays framework prices to open the junction and buys
mismatches instead.

The built-in aligner needs no install. `mafft` or `clustalo` are used instead
with `--aligner`, if they are on PATH.

```bash
python scripts/fetch_oas.py --study "Briney et al., 2019" \
    --genes IGHV3-23 IGHV3-53 --pick 96 --align-rows 50
```

An IMGT/GENE-DB FASTA on disk is preferred as the reference when given
(`--imgt-fasta`); otherwise the reference is a consensus of the germline calls
OAS ships with these reads, and the figure says so.

## Ordering

What gets ordered is the **VH alone** (decision 0016). The default vector
`pcdna-igg1-dual-vk-v1` expresses a full IgG1 from **two cassettes**, each
chain on its own promoter (decision 0020):

    HC (CMV):  CMV — T7 — Kozak — ATG+leader — [VH] — CH1-hinge-CH2-CH3 — stop — BGH pA
    LC (CAG):  CAG — Kozak — ATG+leader — VL — CK — stop — SV40 pA

CAG on the light chain is deliberate — assembly wants LC in excess of HC, and
a single-ORF T2A construct gives the opposite ratio with no way to tune it.
~6.8 kb total, or ~8.4 kb with a mammalian selection cassette.

Constant regions come from UniProt records stored in `fixtures/uniprot/`, not
typed from memory. `pcdna-igg1-t2a-vk-v1` (one ORF, T2A) and
`pcdna-scfv-vk-v1` (VH–(G4S)3–VL–His6) stay registered.

Golden Gate fusion sites sit in the constant flanks — the last 4 nt of the
leader and the first 4 nt of the linker — so one vector accepts all 96 inserts
even though VH N-termini differ by germline. Every clone's assembled reading
frame is translated and checked before the tables are written.

```bash
python scripts/order.py --list-options      # vectors, strategies, vendors
python scripts/order.py --lib LIB-oas-real-... --barcode PLATE-001
python scripts/order.py --lib LIB-... --vendor idt_eblocks --allow-unverified
```

`clone_*` columns describe the molecule; `order_*` columns describe the parcel.
Vendor layouts that have no real template in `fixtures/` refuse to emit unless
asked, and stamp a caveat file when they do — see `docs/formats/vendors.md`.

## Storage

Per `decisions/0006`, the pool is split: full OAS rows go to parquet
(`bulk` table `oas_pool`), and the columns worth filtering on go to the SQLite
`sequences` table. The pool can grow well past what SQLite would handle without
a migration.

`seq_id` is a content hash of source and sequence, so re-ingesting the same unit
adds nothing. It is still minted once per distinct sequence, in `core.ids`, and
never recomputed downstream.

## Open

- No real OAS fixture yet — everything is tested against synthetic units
  (`docs/formats/oas.md` lists what still needs checking against a real file).
- Paired chains, and reformatting onto scFv / VHH / VH-only / CDRH3 grafts,
  are not built. The panel's fixed light chain (IGKV1-39) is declared but unused.
