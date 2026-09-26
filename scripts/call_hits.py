"""Choose a hit-calling rule by measuring it, on SIMULATED plates.

    python scripts/call_hits.py                    # the default, ~600 pairs
    python scripts/call_hits.py --pairs 2000 --seed 7

Every number in `decisions/0024-hit-calling.md` comes out of this script. It
is here so that a claim in a patch can be re-run on the machine that owns the
data rather than believed, and so that changing a threshold means re-running
it rather than arguing.

What it produces, in one run bundle:

  callers.csv         mean precision / recall / F1 per caller
  by_verdict.csv      the same, split by what the plate itself said
  by_dynamic.csv      F1 against top-decile-over-background, ungated
  by_saturation.csv   F1 against how much of the plate hit the ceiling
  per_plate.parquet   one row per plate per caller, everything above
  calls_*.png         two example plates, with the calls drawn on them

Nothing here touches a real measurement. The plates are synthetic and every
output says so.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from plateforge.assay import elisa, figures, hits
from plateforge.core import ids, runs

BANDS = [0, 6, 12, 30, np.inf]
BAND_LABELS = ["<6x", "6-12x", "12-30x", ">30x"]
SAT_BANDS = [-0.001, 0.02, 0.10, 0.20, 1.0]
SAT_LABELS = ["<2%", "2-10%", "10-20%", ">20%"]


def plate_context(pairs) -> pd.DataFrame:
    """The two plate properties a verdict could plausibly rest on.

    Computed here, outside `hits`, so that the thresholds in `hits.DEFAULTS`
    are answerable to a measurement rather than circular.
    """
    rows = []
    for pair in pairs:
        truth = pair.truth()
        target = truth[hits.TARGET].astype(float)
        background = float(np.median(np.sort(target)[:max(3, len(target) // 10)]))
        top = float(target.nlargest(max(1, len(target) // 10)).mean())
        rows.append({
            "pair_id": pair.pair_id,
            "dynamic_range": top / background if background > 0 else np.inf,
            "saturated_fraction": float((target >= hits.DEFAULTS["saturation_od"]).mean()),
            "binders": int(truth["is_binder"].sum()),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=600)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--clones", type=int, default=96)
    ap.add_argument("--label", default="hit-calling")
    args = ap.parse_args()

    clone_ids = [f"CLN-SIM-{i:04d}" for i in range(args.clones)]
    print(f"simulating {args.pairs} pairs x {args.clones} clones")
    pairs = elisa.simulate_many(clone_ids, args.pairs, seed=args.seed)
    context = plate_context(pairs)

    print("scoring every caller, twice: refusing and not refusing")
    gated = hits.evaluate(pairs, respect_verdict=True)
    ungated = hits.evaluate(pairs, respect_verdict=False)

    callers = hits.summarise(gated)

    # The verdict table is computed for ONE caller, not pooled over all four.
    # Pooled, every plate is counted four times and the totals read as plate
    # counts while being row counts -- which is exactly the sort of quiet
    # inflation this script exists to prevent.
    reference = "ratio_over_blank"
    one = ungated[ungated["caller"] == reference]
    grouped = one.groupby("verdict")
    by_verdict = grouped[["precision", "recall", "f1"]].mean().round(3)
    by_verdict["tp"] = grouped["tp"].sum()
    by_verdict["fp"] = grouped["fp"].sum()
    by_verdict["plates"] = grouped["pair_id"].nunique()

    # The banding is computed ungated and on plates that had at least one
    # binder: a plate with nothing to find scores 0 for reasons that have
    # nothing to do with whether it could have been read.
    detail = ungated.merge(context, on="pair_id")
    live = detail[(detail["binders"] > 0) &
                  (detail["caller"] == reference)].copy()
    live["range_band"] = pd.cut(live["dynamic_range"], BANDS, labels=BAND_LABELS)
    live["sat_band"] = pd.cut(live["saturated_fraction"], SAT_BANDS,
                              labels=SAT_LABELS)

    def band(column: str) -> pd.DataFrame:
        grouped = live.groupby(column, observed=True)
        out = grouped[["precision", "recall", "f1"]].mean().round(3)
        out["plates"] = grouped.size()
        return out

    by_dynamic, by_saturation = band("range_band"), band("sat_band")

    print("\n== callers, with refusals respected ==")
    print(callers.to_string())
    print(f"\n== what the verdict was worth ({reference}, nobody refusing) ==")
    print(by_verdict.to_string())
    print(f"\n== F1 by dynamic range ({reference}, plates with >=1 binder) ==")
    print(by_dynamic.to_string())
    print(f"\n== F1 by saturated fraction ({reference}, >=1 binder) ==")
    print(by_saturation.to_string())

    bundle = runs.Bundle(ids.mint_stamped("RUN", args.label), "hit_calling",
                         label=args.label)
    bundle.section("simulation", {
        "SIMULATED": True,
        "warning": "synthetic plates; these numbers calibrate a rule, they "
                   "are not an assay result",
        "pairs": args.pairs, "clones": args.clones, "seed": args.seed,
    })
    bundle.section("thresholds_in_force", hits.DEFAULTS)
    bundle.section("reference_caller", {
        "caller": reference,
        "why": "the verdict and banding tables are computed for one caller, "
               "not pooled over four, so a plate is counted once",
    })
    bundle.section("headline", {
        "best_f1_caller": str(callers.index[0]),
        "refused_fraction": round(float(
            (gated["verdict"] == hits.UNCALLABLE).mean()), 3),
        "f1_callable": float(by_verdict.loc[hits.CALLABLE, "f1"]),
        "f1_degraded": float(by_verdict.loc[hits.DEGRADED, "f1"])
        if hits.DEGRADED in by_verdict.index else None,
        "f1_uncallable_if_not_refused": float(
            by_verdict.loc[hits.UNCALLABLE, "f1"])
        if hits.UNCALLABLE in by_verdict.index else None,
    })

    for name, table in (("callers", callers), ("by_verdict", by_verdict),
                        ("by_dynamic", by_dynamic),
                        ("by_saturation", by_saturation)):
        bundle.write_table(f"{name}.csv", table.reset_index(),
                           role=f"{name.replace('_', ' ')}, simulated")
    detail.to_parquet(bundle.dir / "per_plate.parquet", index=False)
    bundle.add(bundle.dir / "per_plate.parquet",
               role="one row per plate per caller")

    # Two worked examples: a plate the caller could read, and one it refused.
    verdicts = gated[gated["caller"] == "ratio_over_blank"] \
        .set_index("pair_id")["verdict"].to_dict()
    for wanted in (hits.CALLABLE, hits.UNCALLABLE):
        pair = next((p for p in pairs if verdicts.get(p.pair_id) == wanted), None)
        if pair is None:
            continue
        callset = hits.call(pair.truth(), "ratio_over_blank")
        made = figures.hit_calls(pair, callset,
                                 bundle.dir / f"calls_{wanted}.png")
        bundle.add(made, role=f"an example {wanted} plate, calls drawn on it")
        made = figures.elisa_pair(pair, bundle.dir / f"pair_{wanted}.png")
        bundle.add(made, role=f"the same {wanted} pair, target vs control")

    bundle.close()
    print(f"\nrun bundle {bundle.dir}")
    for f in sorted(p.name for p in bundle.dir.iterdir()):
        print(f"  {f}")


if __name__ == "__main__":
    main()
