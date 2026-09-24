"""End-to-end run on SYNTHETIC data, writing figures.

    python scripts/demo.py [--out DIR] [--n 4000] [--pick 10]

Output goes to figures/synthetic/ by default, kept apart from
figures/real/ so a synthetic figure is never mistaken for a real one.
The whole figures/ tree is gitignored; the few committed figures the
README uses live in docs/figures/.

Everything here works offline. Swap `synth.make_unit(...)` for an OAS URL
(`oas.unit_url("Chen_2020", "SRR11937587_1_Heavy_IGHG")`) to run it for real.
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures/synthetic",
                    help="where to write figures (gitignored)")
    ap.add_argument("--n", type=int, default=4000, help="synthetic sequences to generate")
    ap.add_argument("--pick", type=int, default=10, help="how many to sample")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--scratch", action="store_true",
                    help="use a throwaway data dir instead of PLATEFORGE_DATA")
    args = ap.parse_args()

    if args.scratch or not os.environ.get("PLATEFORGE_DATA"):
        os.environ["PLATEFORGE_DATA"] = tempfile.mkdtemp(prefix="plateforge-demo-")
        print(f"data dir: {os.environ['PLATEFORGE_DATA']}  (throwaway)")

    from plateforge.core import artifacts
    from plateforge.library import diversity, figures, germlines, oas, pool, synth

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    unit = synth.make_unit(Path(os.environ["PLATEFORGE_DATA"]) / "synthetic_unit.csv.gz",
                           n=args.n, seed=7)
    df, meta = oas.load(unit, oas.Filter(chain="H", cdr3_len_range=(8, 30)))
    run_id = pool.ingest(df, meta)
    print(f"\ningest {run_id}")
    print(f"  scanned {meta['_rows_scanned']}, kept {meta['_rows_kept']}")
    print("  " + json.dumps(pool.stats()["by_gene"]))

    idx = pool.index()
    picked = diversity.sample(idx, args.pick, panel=germlines.DEFAULT, seed=args.seed)
    lib_id = pool.make_library("demo", list(picked["seq_id"]),
                               produced_by="library.diversity",
                               params={"n": args.pick, "panel": germlines.DEFAULT.name,
                                       "seed": args.seed},
                               parents={run_id: "sampled_from"})
    print(f"\nlibrary {lib_id}")
    print(picked[["seq_id", "v_gene", "cdr3_aa", "cdr3_len"]].to_string(index=False))

    report = diversity.Spec(min_genes=len(germlines.DEFAULT.genes)).check(picked)
    print("\ndiversity report")
    for key, val in report.items():
        if isinstance(val, dict):
            print(f"  {'PASS' if val['pass'] else 'FAIL'}  {key}: "
                  f"{val['value']} (limit {val['limit']})")
    print(f"  OVERALL {'PASS' if report['pass'] else 'FAIL'}")
    print(f"\nlineage up from library: {artifacts.lineage(lib_id, 'up')}")

    made = figures.all_figures(idx, picked, out, panel=germlines.DEFAULT)

    # The alignment view needs the gapped sequence, which lives in parquet.
    full = pool.fetch(list(idx["seq_id"]),
                      columns=["seq_id", "v_gene", "j_gene", "aa_gapped", "germline_aa"]
                              + oas.REGION_COLUMNS)
    busiest = idx["v_gene"].value_counts().index[0]
    made.extend(figures.alignments_per_gene(full, out, prefix="fig5_alignment",
                                            max_rows=50).values())
    made.append(figures.aa_legend(out / "fig5b_aa_legend.png"))

    print("\nfigures:")
    for p in made:
        print(f"  {p}")


if __name__ == "__main__":
    main()
