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
| `figures.py` | Regenerable figures; figure 3 is a visual acceptance test |

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
