# fixtures

Real, unmodified example files from external systems. These are ground truth
for parsers and emitters, and they let the full test suite run without any
instrument present.

Do not hand-edit files in here. If a file needs anonymising, note that in the
corresponding `docs/formats/` entry and keep the original structure intact.

## What is still missing, and what each file would unblock

This is the shopping list. Everything below is currently an assumption, a
refusal, or a number calibrated on simulated data. Each row names the one
file that turns it into something measured. Capture instructions are in the
per-directory READMEs.

| Capture | Where it goes | What it unblocks | Question |
|---|---|---|---|
| **One Gen5 workbook** — a real multi-plate scan, saved as .xlsx exactly as Gen5 writes it, nothing edited | `fixtures/gen5/` | Which sheet is which plate, which read is 450 nm, what the timestamps mean. The grid finder works by shape and needs no fixture; everything *around* the grid is unparsed | Q20 |
| **A plate with a known positive control**, target and control, read at the usual development time | `fixtures/gen5/` | `hits.DEFAULTS["min_signal"] = 0.2` — the weakest number in the repo, reached for rather than measured | Q21 |
| **The same plate read twice**, at the normal time and again over-developed | `fixtures/gen5/` | Tests the saturation finding (0024) directly instead of inferring it from the simulator's own optics. Also settles Q22: whether a second read is the cheap way to rank saturated wells |
| **INTEGRA CSVs**, one of each import type the VIALAB software accepts, plus its Excel template | `fixtures/integra/` | `emit.worklists.integra_vialab`, which currently refuses and says so | Q7 |
| **An Octet quantitation export** plus the plate map that went with it | `fixtures/octet/` | Parsing BLI concentrations back onto wells, so `normalise` runs on real numbers | Q7 |
| **Vendor upload templates** — IDT eBlocks, GenScript, Genewiz, whichever we actually order from | `fixtures/<vendor>/` | `library.vendors`; the IDT layout is a two-column guess that refuses without `allow_unverified=True` | Q15 |
| **A genuine rerun pair** — one sample submitted twice, named as the vendor named them | `fixtures/sanger/` | Rerun grouping is a heuristic on the sample name; one real example would say whether the marker list matches reality | Q23 |
| **A vendor submission manifest** — the sheet tying sample names to wells | `fixtures/sanger/` | Would replace that heuristic entirely | Q23 |
| **A few thousand real OAS rows**, as downloaded | `fixtures/oas/` | The library suite runs entirely on synthetic OAS-format units. Real ingest works; the fixture is what stops it silently breaking | Q9 |
| **Real protocol numbers** — bead and buffer volumes, incubation times, wash counts | a protocol JSON, not a fixture | `assay.protocol` defaults, which are declared placeholders today | Q19 |

Now closed, and worth noting because they show what a single real file buys:
the three Azenta `.ab1` traces, the pcDNA3.1(+) SnapGene map, and the three
vendors' order-form templates. Between them they turned the ABIF reader, the
backbone, and all three order-form emitters from guesses into things checked
against the real thing — and caught two bugs that no amount of reasoning had
(decision 0025).

None of these block building. They block *claiming*: until each arrives, the
thing it would verify is marked `heuristic`, `documented` or `SIMULATED`
rather than presented as known.
