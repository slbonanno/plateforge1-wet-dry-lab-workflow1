"""Generate SIMULATED paired ELISA data. Not real measurements.

    python scripts/simulate_elisa.py --pairs 600
    python scripts/simulate_elisa.py --pairs 2000 --seed 7 --register 3

One pair is one experiment: the same 96 clones read on a target-antigen plate
and on a control-antigen plate (HLA-His), position for position. A clone is a
hit because its target signal exceeds its control signal, not because it is
high on its own -- which is why the pair, and not the plate, is the unit.

Writes a run bundle: the reads, the latent truth, per-plate features for a
classifier, a baseline, and a figure. Everything is stamped `simulated`.

`--register N` also puts N pairs through the real plate lineage, so the
ingest path is exercised rather than assumed.
"""
from __future__ import annotations

import argparse
import json

import pandas as pd

from plateforge.assay import elisa, plates
from plateforge.core import runs


def baseline(features: pd.DataFrame) -> dict:
    """How separable target and control are, by plate pattern alone.

    Reported because the honest answer is 'it depends on the plate'. A panel
    with no hits has a target plate that IS a control plate, and no model can
    or should tell them apart. A classifier that scores well on those is
    reading something it should not.
    """
    try:
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.model_selection import GroupKFold, cross_val_predict
    except ImportError:
        return {"skipped": "scikit-learn not installed; pip install -e '.[ml]'"}

    # A full plate has no blank wells, so its blank-derived features are NaN
    # by construction, not by accident. -1 is outside every real range, so a
    # tree can split on "this plate had no blanks" rather than being handed a
    # value it cannot use.
    X = features.drop(columns=[c for c in ("pair_id", "scenario", "plate_kind",
                                           "is_target", "n_binders")
                               if c in features.columns]).fillna(-1.0)
    y = features["is_target"].astype(int)
    pred = cross_val_predict(GradientBoostingClassifier(random_state=0), X, y,
                             groups=features["pair_id"], cv=GroupKFold(5))
    frame = features.assign(correct=(pred == y))
    by_scenario = frame.groupby("scenario")["correct"].agg(["mean", "size"])
    buckets = pd.cut(frame["n_binders"], [-1, 0, 2, 10, 30, 96])
    by_hits = frame.groupby(buckets, observed=True)["correct"].agg(["mean", "size"])
    return {
        "model": "GradientBoostingClassifier, 5-fold grouped by pair",
        "accuracy": round(float((pred == y).mean()), 4),
        "by_scenario": {k: round(float(v), 4)
                        for k, v in by_scenario["mean"].items()},
        "by_binder_count": {str(k): round(float(v), 4)
                            for k, v in by_hits["mean"].items()},
        "note": "accuracy near chance on plates with no binders is correct, "
                "not a failure: those two plates are the same experiment.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=600)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--clones", type=int, default=96)
    ap.add_argument("--scenario", help="force one scenario for every pair")
    ap.add_argument("--register", type=int, default=0,
                    help="also push this many pairs through the plate store")
    ap.add_argument("--label", default="elisa-sim")
    args = ap.parse_args()

    clone_ids = [f"CLN-SIM-{i:04d}" for i in range(args.clones)]
    print(f"simulating {args.pairs} pairs x {args.clones} clones "
          f"({args.pairs * args.clones * 2:,} wells)")
    pairs = elisa.simulate_many(clone_ids, args.pairs, seed=args.seed,
                                scenario=args.scenario)

    reads = pd.concat([p.long() for p in pairs], ignore_index=True)
    truth = pd.concat([p.truth() for p in pairs], ignore_index=True)
    features = elisa.feature_table(pairs)
    binders = {p.pair_id: sum(b.is_binder for b in p.behaviours.values())
               for p in pairs}
    features["n_binders"] = features["pair_id"].map(binders)

    from plateforge.core import ids
    bundle = runs.Bundle(ids.mint_stamped("RUN", args.label),
                         "simulated_elisa", label=args.label)
    bundle.section("simulation", {
        "SIMULATED": True,
        "warning": "synthetic data; never present as measurements",
        "pairs": args.pairs, "clones": args.clones, "seed": args.seed,
        "wells": int(len(reads)),
        "scenario_forced": args.scenario,
        "scenario_counts": features.drop_duplicates("pair_id")["scenario"]
                                   .value_counts().to_dict(),
        "assumed": elisa.ASSUMED,
        "hit_rate_mix": elisa.HIT_RATE_MIX,
        "control_binder_rate": elisa.CONTROL_BINDER_RATE,
        "anti_tag_rate": elisa.ANTI_TAG_RATE,
    })
    bundle.section("binders", {
        "per_plate_median": float(pd.Series(binders).median()),
        "per_plate_max": int(pd.Series(binders).max()),
        "plates_with_none": int((pd.Series(binders) == 0).sum()),
    })

    reads.to_parquet(bundle.dir / "reads.parquet", index=False)
    bundle.add(bundle.dir / "reads.parquet", role="every well, both plates")
    truth.to_parquet(bundle.dir / "truth.parquet", index=False)
    bundle.add(bundle.dir / "truth.parquet", role="paired wells + latent truth")
    features.to_parquet(bundle.dir / "plate_features.parquet", index=False)
    bundle.add(bundle.dir / "plate_features.parquet",
               role="one row per plate, for the classifier")
    bundle.write_table("reads_sample.csv", reads.head(2000),
                       role="the first 2000 rows, readable without parquet")

    print("\nfitting a baseline classifier...")
    report = baseline(features)
    bundle.section("baseline", report)
    if "accuracy" in report:
        print(f"  accuracy {report['accuracy']:.3f}")
        for key, value in sorted(report["by_binder_count"].items()):
            print(f"    binders {key:<10} {value:.3f}")

    made = elisa.heatmap(pairs, bundle.dir / "example_plates.png", n=4)
    bundle.add(made, role="four example pairs")

    if args.register:
        print(f"\nregistering {args.register} pair(s) through the plate store")
        for pair in pairs[:args.register]:
            for kind, frame in (("target", pair.target),
                                ("control", pair.control)):
                plate = plates.create("assay", 96, label=f"elisa-{kind}",
                                      is_virtual=True,
                                      notes=f"SIMULATED {pair.pair_id} "
                                            f"{pair.scenario}")
                plates.fill(plate, frame, volume_ul=100.0)
                plates.save(plate)
                plates.save_reads(plate.plate_id, frame, channel=elisa.CHANNEL,
                                  units="OD450")
        bundle.section("registered", {"pairs": args.register})

    bundle.write_text("summary.md", summary_markdown(bundle, report, binders),
                      role="README block")
    bundle.close()

    print(f"\nrun bundle {bundle.dir}")
    for f in sorted(p.name for p in bundle.dir.iterdir()):
        print(f"  {f}")


def summary_markdown(bundle, report: dict, binders: dict) -> str:
    counts = pd.Series(binders)
    lines = [
        "### Simulated ELISA data",
        "",
        "**Synthetic. Not measurements.** Generated by `library`-style "
        "simulation so the analysis path can be built before a real plate "
        "reader export exists.",
        "",
        f"- **Design.** {len(binders)} paired experiments. Each pair is the "
        "same 96 clones on a target-antigen plate and a control-antigen "
        "plate (HLA-His), position for position. A clone is a hit when its "
        "target signal beats its own control signal.",
        f"- **Binders.** median {counts.median():.0f} per plate, max "
        f"{counts.max()}, and {int((counts == 0).sum())} plates with none — "
        "a panel that yields nothing is common and is included deliberately.",
        "- **Also modelled.** Clones that bind the control antigen, and "
        "anti-tag clones that light up both plates; expression failures; "
        "TMB over- and under-development; poor washing; weak coating; edge "
        "evaporation; reader gradients; pipetting scatter.",
        "- **Optics.** OD develops toward the reader ceiling, so an "
        "over-developed plate saturates rather than brightening — which is "
        "why it loses its hits.",
    ]
    if "accuracy" in report:
        lines += [
            f"- **Baseline.** Telling a target plate from a control plate by "
            f"signal pattern alone: {report['accuracy']:.0%} accuracy "
            "(gradient boosting, grouped 5-fold). Near chance on plates with "
            "no binders, which is the correct answer — those two plates are "
            "the same experiment.",
        ]
    lines += [
        "",
        "| File | Holds |",
        "|---|---|",
        f"| `{bundle.dir.name}/reads.parquet` | every well of every plate |",
        f"| `{bundle.dir.name}/truth.parquet` | paired wells plus the latent truth |",
        f"| `{bundle.dir.name}/plate_features.parquet` | one row per plate, for training |",
        f"| `{bundle.dir.name}/example_plates.png` | four pairs, as heatmaps |",
        f"| `{bundle.dir.name}/manifest.json` | every parameter, and the code version |",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
