"""Regenerate the figures committed under `docs/figures/`.

    python scripts/readme_figures.py

`docs/figures/` is the only place figures are committed, and the README points
at it. Everything else lands under `plateforge-data/outputs/` and is
gitignored. Mixing the two is the failure mode CLAUDE.md names as the worst
one available: a synthetic figure presented as real.

So this script draws exactly two kinds of figure, and keeps them separable:

**Simulated figures** are reproducible anywhere, from a fixed seed, and carry
SIMULATED on their face. The ELISA pair and the hit-call map are these.

**The pipeline map is not.** Its numbers are how many sequences *this machine*
actually pulled, filtered and ordered, so it is built from the run bundles in
your own data root. If those bundles are not there, the figure is **skipped**
rather than drawn from plausible-looking numbers -- which is the whole point.
Run `scripts/fetch_oas.py` and `scripts/order.py` first, then re-run this.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from plateforge.assay import elisa, figures, hits
from plateforge.core import paths

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "docs" / "figures"
SEED = 17


def newest(prefix: str) -> Path | None:
    root = paths.data_root() / "outputs"
    if not root.exists():
        return None
    found = sorted((d for d in root.iterdir()
                    if d.is_dir() and d.name.startswith(prefix)
                    and (d / "manifest.json").exists()),
                   key=lambda d: d.name)
    return found[-1] if found else None


def simulated_figures() -> list[Path]:
    """A pair worth looking at: full plate, callable, and not trivially easy."""
    clone_ids = [f"CLN-SIM-{i:04d}" for i in range(96)]
    pairs = elisa.simulate_many(clone_ids, 120, seed=SEED)

    best, best_score = None, -1.0
    for pair in pairs:
        truth = pair.truth()
        if len(truth) < 96:                     # a full plate reads clearest
            continue
        callset = hits.call(truth, "ratio_over_blank")
        if callset.verdict != hits.CALLABLE:
            continue
        binders = int(truth["is_binder"].sum())
        if not 4 <= binders <= 14:              # enough to see, few enough to count
            continue
        score = hits.score_calls(callset.calls, truth)
        # Deliberately NOT the best-scoring plate. A README figure showing a
        # perfect call is an advertisement; one with a miss and a false call
        # on it is the honest picture of what this does.
        interesting = score["fn"] + score["fp"]
        if 1 <= interesting <= 4 and score["f1"] > best_score:
            best, best_score = (pair, callset), score["f1"]

    if best is None:
        raise SystemExit("no suitable pair found; widen the search")
    pair, callset = best
    made = [
        figures.elisa_pair(pair, OUT / "fig7_elisa_pair.png"),
        figures.hit_calls(pair, callset, OUT / "fig8_hit_calls.png"),
    ]
    print(f"  from {pair.pair_id} ({pair.scenario}), seed {SEED}")
    return made


def pipeline_figure() -> Path | None:
    """Only from real bundles. No bundles, no figure."""
    lib, clones = newest("LIB-"), newest("CLN-")
    # Naming which bundle is missing matters. "no LIB/CLN run bundle" sent
    # someone looking for a LIB bundle that was sitting right there, three
    # lines above its own name in the same output.
    wanted = [("LIB selection bundle", lib, "scripts/fetch_oas.py"),
              ("CLN ordering bundle", clones, "scripts/order.py")]
    missing = [(what, how) for what, found, how in wanted if found is None]
    if missing:
        print(f"  SKIPPED — no {' and no '.join(w for w, _ in missing)} in "
              f"{paths.data_root() / 'outputs'}.")
        print("  Run " + " and ".join(how for _, how in missing)
              + ", then re-run this. The figure is not drawn from made-up "
                "counts.")
        return None

    lib_manifest = json.loads((lib / "manifest.json").read_text())
    clone_manifest = json.loads((clones / "manifest.json").read_text())
    source = lib_manifest.get("source", {})
    funnel = lib_manifest.get("funnel") or lib_manifest.get("filter_funnel") or []
    selection = lib_manifest.get("selection", {})

    stages = [{"name": "OAS study", "unit": "sequences",
               "count": source.get("rows_fetched") or source.get("limit"),
               "note": str(source.get("study", "")).strip()}]
    for step in funnel if isinstance(funnel, list) else []:
        stages.append({"name": str(step.get("step", "filter")),
                       "unit": "sequences", "count": step.get("kept"),
                       "note": str(step.get("why", "")).strip()})
    stages.append({"name": "selected", "unit": "clones",
                   "count": selection.get("n_selected") or selection.get("pick"),
                   "note": "diversity-aware sampling"})
    stages.append({"name": "ordered", "unit": "fragments",
                   "count": (clone_manifest.get("design", {}) or {}).get("n_clones"),
                   "note": "one 96-well plate"})
    stages = [s for s in stages if s.get("count")]
    if len(stages) < 2:
        print("  SKIPPED — the bundles are there but do not carry counts.")
        return None

    made = figures.pipeline_map(stages, OUT / "fig9_pipeline.png",
                                title="What one run did, on this machine")
    print(f"  from {lib.name} and {clones.name}")
    return made


def input_summary() -> Path | None:
    """Copy the real input-summary figure out of the newest selection bundle.

    This is the one committed figure drawn from *real* OAS data rather than
    from simulation, which is why it is copied rather than regenerated: it
    belongs to a specific run, and the manifest beside it says which one.
    No bundle, no figure -- the same rule as the pipeline map.
    """
    bundle = newest("LIB-")
    if bundle is None:
        print("  SKIPPED -- no LIB run bundle in "
              f"{paths.data_root() / 'outputs'}. Run scripts/fetch_oas.py.")
        return None
    source = bundle / "input_summary.png"
    if not source.exists():
        print(f"  SKIPPED -- {bundle.name} has no input_summary.png")
        return None
    target = OUT / "input_summary.png"
    target.write_bytes(source.read_bytes())
    print(f"  from {bundle.name}")
    return target


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-pipeline", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    print("real input summary (from a selection bundle):")
    made = []
    got = input_summary()
    if got:
        made.append(got)

    print("simulated figures (stamped SIMULATED, fixed seed):")
    made += simulated_figures()
    if not args.skip_pipeline:
        print("pipeline figure (real counts only):")
        got = pipeline_figure()
        if got:
            made.append(got)

    print("\nwrote:")
    for path in made:
        print(f"  {path.relative_to(REPO)}  ({path.stat().st_size // 1024} KB)")
    print("\nREADME.md is yours — these are just on disk for you to point at.")


if __name__ == "__main__":
    main()
