"""Turn picked sequences into orderable constructs, a plate map and an alignment.

    python scripts/order.py --lib LIB-oas-real-...        # from a saved pick
    python scripts/order.py --study "Briney et al., 2019" --pick 96

Produces, under outputs/<clone-set-id>/:

    clones.csv      the 2B metadata table, one row per clone
    plate_map.csv   96-well layout, as the vendor ships it
    plate_map.txt   the same as a grid, for reading
    alignment.png   the clones that will be synthesised

Vendor order forms are deliberately absent: those need a real example of each
vendor's template in fixtures/ first (Q15), and Gibson/HiFi adapters need the
destination vector (Q16).
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PLATEFORGE_DATA", os.path.expanduser("~/plateforge-data"))

from plateforge.core import artifacts, paths, wells                      # noqa: E402
from plateforge.library import (construct, diversity, figures, germline_db,  # noqa: E402
                                germlines, oas, pool)


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
    ap.add_argument("--reserve", nargs="*", default=[],
                    help="wells to keep free, e.g. A01 H12")
    ap.add_argument("--imgt-fasta", help="IMGT FASTA for the light chain")
    args = ap.parse_args()

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

    # --- the light chain, and where it came from --------------------------
    vl = None
    if construct.FORMATS.meta(args.format)["needs_light_chain"]:
        cache = paths.raw_dir("imgt") / "imgt_reference_aa.fasta"
        fasta = args.imgt_fasta or (str(cache) if cache.exists() else None)
        vl = construct.light_chain_for(fasta_path=fasta)
        if vl is None:
            sys.exit(
                f"could not resolve {germlines.DEFAULT.light_chain}.\n"
                "Either run  python scripts/germlines.py --download\n"
                "or ingest light chains:  python scripts/fetch_oas.py "
                "--chain light --any-gene --study \"<a study>\"")
        print(f"light chain: {vl.summary()}")
        if vl.is_derived:
            print("  WARNING: this is a consensus inferred from the pool, not "
                  "IMGT. It covers only the span the reads covered, which for "
                  "5' truncated amplicons is not all of FR1.")

    # --- design -----------------------------------------------------------
    clones = construct.build(
        picked, vl.sequence if vl else "", fmt=args.format,
        linker=args.linker, orientation=args.orientation,
        enzymes=(args.enzyme,),
        germline_source=vl.source if vl else None)

    plated = construct.to_plate(clones, 96, reserved=args.reserve)
    set_id = construct.register(plated, parents={parent: "designed_from"},
                                params={"format": args.format,
                                        "linker": args.linker,
                                        "orientation": args.orientation,
                                        "enzyme": args.enzyme,
                                        "light_chain_source": vl.source if vl else None})
    out = paths.output_dir(set_id)

    # --- report -----------------------------------------------------------
    print(f"\nclone set {set_id}")
    print(f"  construct: {args.format}"
          + (f", {args.orientation}, {args.linker}" if args.format == "scfv" else ""))
    print(f"  protein length: {plated['aa_length'].min()}-{plated['aa_length'].max()} aa")
    if "gc" in plated:
        print(f"  GC: {plated['gc'].min():.1%}-{plated['gc'].max():.1%}"
              f" (median {plated['gc'].median():.1%})")
        print(f"  {args.enzyme} sites removed: {int(plated['sites_removed'].sum())}")
        flagged = plated[plated["synthesis_warnings"].notna()]
        if len(flagged):
            print(f"  {len(flagged)} clones carry synthesis warnings:")
            for _, r in flagged.head(5).iterrows():
                print(f"    {r['clone_id']}: {r['synthesis_warnings']}")
        else:
            print("  no synthesis warnings")

    plated["clone_id_tail"] = plated["clone_id"].str.rsplit("-", n=1).str[-1]
    plated.to_csv(out / "clones.csv", index=False)
    plated[["well", "clone_id", "seq_id", "v_gene", "cdr3_aa",
            "dna_length", "gc"]].to_csv(out / "plate_map.csv", index=False)
    (out / "plate_map.txt").write_text(grid(plated) + "\n")
    (out / "plate_map_ids.txt").write_text(grid(plated, label="clone_id_tail") + "\n")

    full = pool.fetch(list(plated["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa",
                               "anarci_numbering"])
    gene = plated["v_gene"].value_counts().index[0]
    figures.alignment(full, gene, out / f"alignment_{gene}.png", max_rows=24)

    print(f"\nlineage: {artifacts.lineage(set_id, 'up')}")
    print(f"\nwritten to {out}")
    for f in sorted(out.iterdir()):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
