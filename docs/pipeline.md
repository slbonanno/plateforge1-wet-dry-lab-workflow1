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

## 4 — Worklists (`plateforge.emit`)

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

## Reading an unexplained directory

```python
from plateforge.core import runs
print(runs.describe("LIB-..."))
```

Every `manifest.json` records the git SHA that produced it.
