"""Reproduce every result this repo claims, on this machine.

    python scripts/run_all.py              # the fast path, a few minutes
    python scripts/run_all.py --full       # real OAS download included

Everything I assert in a patch should be checkable here rather than taken on
trust. Each stage prints what it produced and where, and the whole run writes
a bundle recording the code version, the timings and the headline numbers --
so `plateforge-data/outputs/` on your machine holds the record of the work,
not a transcript somewhere else.

Stages that need the network are skipped unless asked for, so the default run
works on a plane.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from plateforge.core import ids, paths, runs, version

REPO = Path(__file__).resolve().parents[1]


def run(label: str, argv: list[str]) -> dict:
    print(f"\n=== {label} " + "=" * max(0, 60 - len(label)))
    started = time.time()
    done = subprocess.run([sys.executable] + argv, cwd=REPO,
                          capture_output=True, text=True)
    elapsed = round(time.time() - started, 1)
    tail = done.stdout.strip().splitlines()[-18:]
    print("\n".join(tail))
    if done.returncode != 0:
        print(done.stderr.strip()[-2000:])
    print(f"--- {label}: {'ok' if done.returncode == 0 else 'FAILED'} "
          f"in {elapsed}s")
    return {"label": label, "argv": argv, "ok": done.returncode == 0,
            "seconds": elapsed,
            "tail": tail if done.returncode == 0 else done.stderr[-2000:]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="include the stages that download real OAS data")
    ap.add_argument("--pairs", type=int, default=400)
    ap.add_argument("--pick", type=int, default=96)
    ap.add_argument("--skip-tests", action="store_true")
    args = ap.parse_args()

    print(f"plateforge  code {version.code_version()}")
    print(f"data root   {paths.data_root()}")

    stages = []
    if not args.skip_tests:
        stages.append(run("test suite", ["-m", "pytest", "-q"]))
    stages.append(run("synthetic pipeline", ["scripts/demo.py"]))
    stages.append(run("simulated ELISA",
                      ["scripts/simulate_elisa.py", "--pairs", str(args.pairs),
                       "--register", "2"]))
    stages.append(run("hit calling",
                      ["scripts/call_hits.py", "--pairs", str(args.pairs)]))
    stages.append(run("sequencing fixtures", ["scripts/sequencing.py", "inspect"]))
    if args.full:
        stages.append(run("real OAS selection", [
            "scripts/fetch_oas.py", "--study", "Briney et al., 2019",
            "--limit", "20000", "--genes", "IGHV3-23", "IGHV3-53",
            "--pick", str(args.pick), "--align-rows", "50"]))

    bundle = runs.Bundle(ids.mint_stamped("RUN", "run-all"), "reproduction",
                         label="run_all")
    bundle.section("environment", {
        "code_version": version.code_version(),
        "python": sys.version.split()[0],
        "data_root": str(paths.data_root()),
        "full": args.full,
    })
    bundle.section("stages", [{k: v for k, v in s.items() if k != "tail"}
                              for s in stages])
    bundle.write_text("transcript.txt", "\n\n".join(
        f"=== {s['label']} ({'ok' if s['ok'] else 'FAILED'}, {s['seconds']}s)\n"
        + ("\n".join(s["tail"]) if isinstance(s["tail"], list) else s["tail"])
        for s in stages), role="what each stage printed")
    bundle.close()

    failed = [s["label"] for s in stages if not s["ok"]]
    print(f"\nrecord written to {bundle.dir}")
    print("all stages passed" if not failed else f"FAILED: {', '.join(failed)}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
