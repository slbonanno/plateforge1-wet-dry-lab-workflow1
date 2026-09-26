"""Order Sanger sequencing for a clone plate, and verify what comes back.

Two directions, one script.

    # going out: fill a vendor's order form from a clone table
    python scripts/sequencing.py order --clones plateforge-data/outputs/CLN-.../clones.csv \
        --vendor genewiz --primer CMV_F

    # coming back: judge a directory of .ab1 files against those clones
    python scripts/sequencing.py verify --clones .../clones.csv --traces ~/Downloads/seq/

    # neither: look at the real fixture traces and the backbone
    python scripts/sequencing.py inspect

Every run writes a bundle under `plateforge-data/outputs/`, so the record of
what was submitted and what came back lives on the machine that did it.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from plateforge.core import ids, runs
from plateforge.emit import sequencing as forms
from plateforge.library import abif, dna, sanger, snapgene
from plateforge.library import vector as vectors

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "fixtures"


def load_clones(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = {"clone_id", "well"} - set(frame.columns)
    if missing:
        raise SystemExit(f"{path}: missing column(s) {', '.join(sorted(missing))}")
    return frame


# --- ordering ----------------------------------------------------------------

def do_order(args) -> None:
    frame = load_clones(args.clones)
    frame = frame.assign(length_bp=args.plasmid_bp,
                         concentration_ng_ul=args.concentration)
    form = forms.order(frame, args.vendor, primer=args.primer,
                       allow_issues=args.force)

    bundle = runs.Bundle(ids.mint_stamped("RUN", "seq-order"), "sequencing_order",
                         label=f"seq-order-{args.vendor}")
    bundle.section("order", {
        "vendor": args.vendor, "primer": args.primer,
        "samples": form.n_samples, "confidence": form.confidence,
        "template": form.template, "notes": form.notes, "issues": form.issues,
        "source": str(args.clones),
    })
    written = form.write(bundle.dir, stem="order")
    bundle.add(written, role=f"{args.vendor} order form, ready to upload")
    bundle.write_table("submitted.csv",
                       frame[["well", "clone_id"]].assign(primer=args.primer),
                       role="what went on the plate, in our own format")
    bundle.close()

    print(f"{args.vendor}: {form.n_samples} samples, primer {args.primer}")
    for note in form.notes:
        print(f"  - {note}")
    for issue in form.issues:
        print(f"  ! {issue}")
    print(f"\nrun bundle {bundle.dir}")
    print(f"  {written.name}")


# --- verifying ---------------------------------------------------------------

def do_verify(args) -> None:
    frame = load_clones(args.clones)
    column = next((c for c in ("cds_dna", "insert_dna", "dna", "cds")
                   if c in frame.columns), None)
    if column is None:
        raise SystemExit(
            f"{args.clones}: no column holding the ordered DNA "
            "(looked for cds_dna, insert_dna, dna, cds)")

    vector = vectors.get(args.vector) if args.vector else None
    expectations = sanger.expectations_from(frame, vector, dna_column=column)

    paths = sorted(Path(args.traces).glob("*.ab1"))
    if not paths:
        raise SystemExit(f"no .ab1 files under {args.traces}")
    traces = abif.read_all(paths)
    print(f"{len(traces)} trace(s) from {args.traces}")

    calls, report = sanger.verify(traces, expectations, by=args.match,
                                  prefer=args.prefer)
    counts = sanger.summarise(calls)
    good = sanger.usable_wells(calls)

    print("\n== verdicts ==")
    print(counts.to_string(index=False))
    print(f"\n{len(good)} of {len(expectations)} wells usable")
    if report is not None and len(report.superseded):
        print(f"\n{len(report.superseded)} trace(s) superseded by a rerun")
    if report is not None and report.unresolved:
        print(f"{len(report.unresolved)} trace(s) a human needs to look at")

    trouble = calls[~calls["verdict"].isin(sanger.USABLE)]
    if len(trouble):
        print("\n== wells to look at ==")
        print(trouble[["well", "clone_id", "verdict", "notes"]]
              .to_string(index=False, max_colwidth=80))

    bundle = runs.Bundle(ids.mint_stamped("RUN", "seq-verify"),
                         "sequencing_verification", label="seq-verify")
    bundle.section("verification", {
        "traces": len(traces), "wells_expected": len(expectations),
        "wells_usable": len(good), "matched_by": args.match,
        "vector": args.vector or "none (insert is its own context)",
        "verdicts": counts.set_index("verdict")["wells"].to_dict(),
    })
    bundle.write_table("calls.csv", calls, role="one row per trace, judged")
    bundle.write_table("usable_wells.csv", pd.DataFrame({"well": good}),
                       role="the wells worth carrying forward")
    if report is not None:
        bundle.write_table("reruns.csv", report.table,
                           role="which trace was used for each sample, and why")
        bundle.section("reruns", {
            "kept": len(report.keep), "superseded": len(report.superseded),
            "unresolved": len(report.unresolved),
            "warning": "grouping is by sample name and is a heuristic; see "
                       "decision 0025 and check reruns.csv",
        })
    bundle.close()
    print(f"\nrun bundle {bundle.dir}")


# --- inspecting the fixtures -------------------------------------------------

def do_inspect(_args) -> None:
    """What the real files in `fixtures/` actually contain."""
    print("== backbone ==")
    plasmid = snapgene.read(FIXTURES / "backbone" / "pcDNA3.1.dna")
    print(f"  {plasmid.name}: {plasmid.summary()}")
    problems = snapgene.check(plasmid, 5428, ("CMV promoter", "MCS", "AmpR"))
    print(f"  checks: {'ok' if not problems else problems}")
    for feature in plasmid.features[:6]:
        print(f"    {feature.start:5d}-{feature.end:5d}  {feature.name}")

    print("\n== traces ==")
    traces = abif.read_all(sorted((FIXTURES / "sanger").glob("*.ab1")))
    print(pd.DataFrame([t.summary() for t in traces]).to_string(index=False))

    good = [t for t in traces if t.bases_at_least(20) > 300]
    if len(good) >= 2:
        first, second = good[0], good[1]
        a_start, a_end = abif.mott_trim(first.quality)
        b_start, b_end = abif.mott_trim(second.quality)
        aligned = dna.align_to(
            dna.reverse_complement(second.sequence[b_start:b_end]),
            first.sequence[a_start:a_end])
        print(f"\n{first.sample} against reverse-complement of {second.sample}:")
        print(f"  {aligned.length} bp overlap, identity {aligned.identity:.4f}, "
              f"{aligned.mismatches} mismatches, {aligned.gaps} gaps")
        print("  (two reads of one molecule from opposite ends; a perfect "
              "overlap is the parser, the trim and the aligner all agreeing)")

    print("\n== order forms ==")
    print(forms.available().to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    order = sub.add_parser("order", help="fill a vendor order form")
    order.add_argument("--clones", required=True)
    order.add_argument("--vendor", default=forms.DEFAULT,
                       choices=sorted(forms.FORMS.keys()))
    order.add_argument("--primer", default="CMV_F")
    order.add_argument("--plasmid-bp", type=int, default=6800)
    order.add_argument("--concentration", type=float, default=60.0)
    order.add_argument("--force", action="store_true",
                       help="write the form even if it breaks the vendor's limits")
    order.set_defaults(func=do_order)

    verify = sub.add_parser("verify", help="judge .ab1 files against a clone table")
    verify.add_argument("--clones", required=True)
    verify.add_argument("--traces", required=True)
    verify.add_argument("--vector", default="pcdna-igg1-dual-vk-v1")
    verify.add_argument("--match", default="well", choices=("well", "sample"))
    verify.add_argument("--prefer", default="later", choices=("later", "best"))
    verify.set_defaults(func=do_verify)

    inspect = sub.add_parser("inspect", help="what the real fixtures contain")
    inspect.set_defaults(func=do_inspect)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
