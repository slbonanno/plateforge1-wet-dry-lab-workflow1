"""Turn picked sequences into orderable DNA, a plate map and an order table.

    python scripts/order.py --lib LIB-oas-real-...        # from a saved pick
    python scripts/order.py --lib LIB-... --vendor idt_eblocks --allow-unverified
    python scripts/order.py --list-options                # vectors, strategies, vendors

Produces, under outputs/<clone-set-id>/:

    clones.csv          one row per clone: what is ordered AND what is
                        expressed, the vector, the fusion sites, every warning
    plate_map.csv       96-well layout, wells grouped so neighbours are similar
    plate_map.txt       the same as a grid, for reading
    metadata_table.png  the clone record as a figure, first 10 wells
    summary.md          a README block: vector, cloning, codon optimisation,
                        how the assembled ORF was verified
    order_generic.csv   our own order table -- columns we defined
    order_<vendor>.csv  the vendor's layout, only with --allow-unverified until
                        a real template is in fixtures/ (Q15)
    vector.json         the insertion context these fragments were designed for
    alignment_<gene>.png  the clones that will be synthesised

What is ordered is the VH alone: the working vector already carries the Kozak,
the secretion leader, the (G4S)3 linker, the light chain and the tag. See
library/vector.py and decision 0016.
"""
from __future__ import annotations

import argparse
import os
import sys


import json                                                               # noqa: E402

from plateforge.core import artifacts, paths, wells                      # noqa: E402
from plateforge.library import (cloning, construct, diversity, figures,  # noqa: E402
                                germline_db, germlines, ordering, pool,
                                vector as vectors, vendors)


def grid(plated, plate_format=96, label="cdr3_aa"):
    """A plate as a human reads it.

    Clone ids share a long prefix, so truncating from the front prints the
    same string in every well. The CDRH3 is what actually distinguishes one
    well from the next at the bench; the id's unique tail is the fallback.
    """
    n_rows, n_cols = wells.dims(plate_format)
    if label in plated.columns:
        text = plated[label].astype(str)
    else:
        text = plated["clone_id"].str.rsplit("-", n=1).str[-1]
    cell = 12                      # fixed, so long CDRH3s cannot break the grid
    def fit(value: str) -> str:
        return value if len(value) <= cell else value[:cell - 1] + "\u2026"

    by_well = {w: fit(v) for w, v in zip(plated["well"], text)}

    out = ["     " + "".join(f"{c:>{cell + 1}d}" for c in range(1, n_cols + 1))]
    for r in range(1, n_rows + 1):
        row = wells.row_label(r)
        cells = [f"{by_well.get(f'{row}{c:02d}', '-'):>{cell + 1}}"
                 for c in range(1, n_cols + 1)]
        out.append(f"{row:>3}  " + "".join(cells))
    out.append("")
    out.append(f"{len(plated)} clones, {label} shown, truncated to {cell} "
               f"characters. Empty wells marked '-'. clones.csv has the full values.")
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", help="an existing LIB artifact of picked sequences")
    ap.add_argument("--study", help="pick fresh from this study instead")
    ap.add_argument("--pick", type=int, default=96)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--format", default="scfv", choices=["scfv", "vh_only"])
    ap.add_argument("--linker", default="G4S3")
    ap.add_argument("--orientation", default="VH-VL", choices=["VH-VL", "VL-VH"])
    ap.add_argument("--enzyme", default="BsaI")
    ap.add_argument("--vector", default="pcdna-igg1-dual-vk-v1",
                    help="insertion context; --list-options shows them")
    ap.add_argument("--strategy", default="golden_gate",
                    choices=sorted(cloning.STRATEGIES.keys()))
    ap.add_argument("--vendor", default="generic",
                    help="order table layout; --list-options shows them")
    ap.add_argument("--allow-unverified", action="store_true",
                    help="emit a vendor layout that has no real template in "
                         "fixtures/ yet. It will be stamped as unverified.")
    ap.add_argument("--min-order-length", type=int, default=300,
                    help="pad short fragments up to this (IDT eBlocks minimum)")
    ap.add_argument("--no-group", action="store_true",
                    help="lay the plate out in pick order rather than "
                         "grouping similar clones together")
    ap.add_argument("--barcode", help="plate barcode to stamp on every well")
    ap.add_argument("--list-options", action="store_true")
    ap.add_argument("--reserve", nargs="*", default=[],
                    help="wells to keep free, e.g. A01 H12")
    ap.add_argument("--imgt-fasta", help="IMGT FASTA for the light chain")
    args = ap.parse_args()

    if args.list_options:
        print("vectors:")
        for name in vectors.VECTORS:
            print("  " + vectors.get(name).describe().replace("\n", "\n  "))
        print("\ncloning strategies:")
        for name in cloning.STRATEGIES:
            meta = cloning.STRATEGIES.meta(name)
            print(f"  {name:<14} enzyme={meta.get('enzyme')} "
                  f"scarless={meta.get('scarless')}")
        print("\nvendors:")
        print(vendors.available().to_string(index=False))
        return

    # --- the sequences ----------------------------------------------------
    if args.lib:
        picked = pool.library_members(args.lib)
        if picked.empty:
            sys.exit(f"no members in {args.lib}")
        parent = args.lib
    else:
        idx = pool.index()
        if idx.empty:
            sys.exit("the pool is empty; run scripts/fetch_oas.py first")
        picked = diversity.sample(idx, args.pick, seed=args.seed)
        parent = pool.make_library("order", list(picked["seq_id"]),
                                   produced_by="library.diversity",
                                   params={"n": args.pick, "seed": args.seed})
    print(f"{len(picked)} sequences from {parent}")

    # --- design -----------------------------------------------------------
    vec = vectors.get(args.vector)
    print("\n" + vec.describe())
    problems = vec.check_fusion_sites()
    for problem in problems:
        print(f"  ! {problem}")

    plated = ordering.design(
        picked, vector_name=args.vector, strategy=args.strategy,
        min_order_length=args.min_order_length, enzymes=(args.enzyme,),
        fmt=args.format, group_by_similarity=not args.no_group,
        reserved=args.reserve, barcode=args.barcode)

    set_id = construct.register(plated, parents={parent: "designed_from"},
                                params={"format": args.format,
                                        "vector": args.vector,
                                        "strategy": args.strategy,
                                        "enzyme": args.enzyme,
                                        "min_order_length": args.min_order_length,
                                        "grouped": not args.no_group})
    out = paths.output_dir(set_id)

    # --- report -----------------------------------------------------------
    rep = ordering.report(plated)
    print(f"\nclone set {set_id}")
    for key, val in rep.items():
        print(f"  {key}: {val}")
    flagged = plated[plated["synthesis_warnings"].notna()]
    for _, r in flagged.head(5).iterrows():
        print(f"    {r['clone_id']}: {r['synthesis_warnings']}")
    broken = plated[plated["assembly_problems"].notna()]
    if len(broken):
        print(f"\n  {len(broken)} clone(s) do not assemble cleanly:")
        for _, r in broken.head(5).iterrows():
            print(f"    {r['clone_id']}: {r['assembly_problems']}")

    # --- write ------------------------------------------------------------
    plated["clone_id_tail"] = plated["clone_id"].str.rsplit("-", n=1).str[-1]
    plated.to_csv(out / "clones.csv", index=False)
    plated[["well", "clone_id", "seq_id", "v_gene", "cdr3_aa",
            "order_length", "order_gc"]].to_csv(out / "plate_map.csv", index=False)
    (out / "plate_map.txt").write_text(grid(plated) + "\n")
    (out / "plate_map_ids.txt").write_text(grid(plated, label="clone_id_tail") + "\n")
    (out / "vector.json").write_text(json.dumps(vec.as_dict(), indent=2) + "\n")

    wanted = {"generic", args.vendor}
    emitted = {}
    for name in sorted(wanted):
        try:
            order = vendors.emit(plated, name,
                                 allow_unverified=args.allow_unverified)
        except ValueError as exc:
            print(f"\n  {name}: {exc}")
            continue
        order.table.to_csv(out / f"order_{name}.csv", index=False)
        emitted[name] = order.as_dict()
        if not order.template_verified:
            (out / f"order_{name}.CAVEAT.txt").write_text(
                "\n".join(order.caveats) + "\n")
        if len(order.issues):
            order.issues.to_csv(out / f"order_{name}_issues.csv", index=False)
            print(f"\n  {name}: {order.issues['name'].nunique()} clone(s) "
                  f"outside the product limits — see order_{name}_issues.csv")
    (out / "order_summary.json").write_text(
        json.dumps({"report": rep, "vendors": emitted}, indent=2) + "\n")

    figures.metadata_table(plated, 10, out / "metadata_table.png")
    (out / "summary.md").write_text(
        ordering.summary_markdown(plated, args.vector, set_id, args.strategy))

    full = pool.fetch(list(plated["seq_id"]),
                      columns=["seq_id", "v_gene", "j_gene", "aa_gapped",
                               "germline_aa"])
    # One alignment per V gene on the plate, not just the busiest: a plate
    # with three germlines on it needs three pictures to be described.
    made = list(figures.alignments_per_gene(full, out, max_rows=50).values())
    figures.stamp(out, run_id=set_id, note="clone set")

    print(f"\nREADME block ready to paste: {out / 'summary.md'}")
    print(f"\nlineage: {artifacts.lineage(set_id, 'up')}")
    print(f"\nwritten to {out}")
    for f in sorted(out.iterdir()):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
