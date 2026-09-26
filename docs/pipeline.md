# The pipeline, so far

Two commands. Each writes one self-describing directory under
`$PLATEFORGE_DATA/outputs/`, keyed on the artifact id it produced. Nothing
overwrites anything: a re-run mints a new id and a new directory.

```bash
python scripts/fetch_oas.py --study "Briney et al., 2019" --limit 20000 \
    --genes IGHV3-23 IGHV3-53 --pick 96 --align-rows 50
python scripts/order.py --lib LIB-oas-real-... --barcode PLATE-001 \
    --vendor idt_eblocks --allow-unverified
```

## 1 — Sequences → a selection (`scripts/fetch_oas.py`)

| Stage | Does | Code |
|---|---|---|
| fetch | OAS unpaired heavy, via the HuggingFace parquet mirror (OPIG serves 403; OAS has no API) | `library.hf` |
| normalise | AIRR columns → our schema; V/J calls parsed once | `library.oas`, `library.germlines` |
| filter | named steps with a survival funnel, so "kept 0" is never the only clue | `library.oas.Filter` |
| store | parquet for full rows, SQLite for the filter index | `library.pool` |
| select | greedy farthest-point on normalised CDRH3 edit distance, one CDRH3 per clone | `library.diversity` |
| align | each V gene to its own IMGT germline, gaps relaxed inside CDRs | `library.msa` |
| record | the whole chain, as files | `library.selection`, `core.runs` |

**Out:** `outputs/<LIB id>/` — `manifest.json`, `selected.csv`,
`pool_summary.csv`, `filter_funnel.csv`, `alignment_<gene>.png|.fasta|.aln.txt`,
`input_summary.png`, **`summary.md`** (paste into the README).

## 2 — A selection → orderable DNA (`scripts/order.py`)

| Stage | Does | Code |
|---|---|---|
| layout | A01–H12 column-major, picks reordered so neighbouring wells are similar | `library.diversity.group_for_plate` |
| design | VH only; the vector carries leader, CH1–CH3, T2A, VL, CK | `library.vector` |
| optimise | human codon usage, mid-window GC, BsaI removal, homopolymer breaking | `library.codon` |
| join | Golden Gate fragment: pad–BsaI–N–OH–CDS–OH–N–BsaI–pad | `library.cloning` |
| verify | the **assembled ORF** translated and checked, not the fragment alone | `library.cloning.verify` |
| emit | our table always; a vendor layout only when it has a real template | `library.vendors` |

**Out:** `outputs/<CLN set id>/` — `clones.csv`, `plate_map.csv|.txt`,
`order_generic.csv`, `order_<vendor>.csv` (+ `CAVEAT.txt`), `vector.json`,
`metadata_table.png`, `alignment_<gene>.png`, **`summary.md`**.

## 3 — Plates through the lab (`plateforge.assay`)

The clone is tracked; plates are containers it passes through; steps are
registry entries (decision 0019).

```python
from plateforge.assay import plates, steps, protocol

plan = steps.get("harvest_supernatant").plan(culture, volume_ul=500.0)
steps.commit(plan)
steps.record(plan, operator="shivan", lots=["LOT-42"],
             deviations="column 7 short")

plates.trace(igg.plate_id, "D07")        # every well it descends from
plates.layout_matches(ordered, igg)      # assert the map survived
```

The chain: `transfect → harvest_supernatant → add_beads → immunoprecipitate →
elute_neutralise → magnet_transfer → quantify_bli → normalise`. `plan` has no
side effects; `record` takes what happened. A human step has the same shape as
a robot step.

## 4 — Sequencing the clones (`scripts/sequencing.py`)

```bash
python scripts/sequencing.py order  --clones .../clones.csv --vendor genewiz --primer CMV_F
python scripts/sequencing.py verify --clones .../clones.csv --traces ~/Downloads/seq/
python scripts/sequencing.py inspect
```

| Stage | Does | Code |
|---|---|---|
| order | fills a vendor's real order form; enforces its stated limits before upload | `emit.sequencing` |
| read | `.ab1` traces: bases, Phred scores, and the well/plate/run-time the sequencer wrote | `library.abif` |
| trim | Mott's algorithm, so 5' noise is not aligned as mutations | `library.abif.mott_trim` |
| screen | each read against all 96 expected **inserts**, to catch a plate swap | `library.dna.screen` |
| align | banded affine alignment, local at both ends | `library.dna.align_to` |
| call | `exact` / `silent` / `missense` / `indel` / `frameshift` / `mixed` / `no_insert` / `wrong_clone` / `low_quality` | `library.sanger` |
| reruns | grouped by sample name, ordered by the instrument's timestamp; superseded, never deleted | `library.sanger.resolve_reruns` |

Every format here is verified against a real file (decision 0025): three
Azenta `.ab1` traces, the pcDNA3.1(+) SnapGene map, and the three vendors'
order-form templates. **Out:** `outputs/<RUN id>/` — `calls.csv`,
`usable_wells.csv`, `reruns.csv`, and the filled order form.

## 5 — Reads → hits (`plateforge.assay.hits`)

```python
from plateforge.assay import hits

result = hits.call(paired, "ratio_and_z", blanks=empty_wells)
result.verdict      # callable | degraded | uncallable
result.reasons      # sentences, not codes
result.hits
```

A hit is a well that beats **its own** control well. The caller is allowed to
refuse, and refusing is right about 10% of the time (decision 0024): on a
plate whose top decile is under 10× its background, no threshold works, and
the honest output is "re-read this sooner" rather than a short confident
list. Saturation is flagged, never refused on — measurement said that one the
other way round from the obvious guess.

```bash
python scripts/call_hits.py --pairs 600 --seed 17
```

Re-measures every threshold and writes the tables plus two worked plates.
**Out:** `outputs/<RUN id>/` — `callers.csv`, `by_verdict.csv`,
`by_dynamic.csv`, `by_saturation.csv`, `per_plate.parquet`,
`calls_callable.png`, `calls_uncallable.png`, `pair_*.png`.

## 6 — Worklists (`plateforge.emit`)

```python
from plateforge.emit import worklists
worklists.available()                                   # and how well we know each
worklists.emit(plan.transfer_frame, "tecan_evo", allow_unverified=True)
```

Every emitter declares `verified` / `documented` / `sketch` (decision 0021).
Nothing is verified but our own format. Which instrument this lab has lives in
`core.lab`, read from `plateforge-data/lab.json`.

## Blocked on real files

| | Needs |
|---|---|
| INTEGRA worklists | one exported CSV of each import type → `fixtures/integra/` |
| Octet BLI | one quantitation export + its plate map → `fixtures/octet/` |
| Vendor order forms | one upload template each → `fixtures/<vendor>/` |
| Protocol numbers | real volumes, incubations, wash counts (Q19) |
| Gen5 metadata | one real multi-sheet workbook → `fixtures/gen5/` (Q20) |
| Hit thresholds | a plate with a known positive control (Q21) |
| Rerun grouping | a vendor submission manifest → `fixtures/sanger/` (Q23) |

## Reading an unexplained directory

```python
from plateforge.core import runs
print(runs.describe("LIB-..."))
```

Every `manifest.json` records the git SHA that produced it.
