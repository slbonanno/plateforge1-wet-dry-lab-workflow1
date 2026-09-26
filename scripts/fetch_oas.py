"""Ingest real OAS sequences into the pool.

    python scripts/fetch_oas.py --list-studies
    python scripts/fetch_oas.py --study "Briney et al., 2019" --limit 20000
    python scripts/fetch_oas.py --source ~/Downloads/some_unit.csv.gz

Three ways in, in order of preference:

  --study    the HuggingFace parquet mirror of OAS (default, and the only
             route that currently works: OPIG serves 403 for its own .csv.gz
             files, and OAS has no API)
  --source   a local .csv.gz data unit, or a URL to one
  --list-studies   what the mirror holds

Figures go to figures/real/, kept apart from the synthetic output of
scripts/demo.py. The whole figures/ tree is gitignored.
"""
from __future__ import annotations

import argparse
import sys


from plateforge.library import (diversity, figures, germlines, msa, oas,  # noqa: E402
                                pool, selection)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", help="study name in the HuggingFace mirror")
    ap.add_argument("--source", help="a local .csv.gz data unit, or a URL to one")
    ap.add_argument("--list-studies", action="store_true")
    ap.add_argument("--chain", default="heavy", choices=["heavy", "light"])
    ap.add_argument("--species", default="human",
                    help="mirror only; 'any' keeps every species")
    ap.add_argument("--limit", type=int, default=20000,
                    help="stop after this many kept sequences (0 = no cap)")
    ap.add_argument("--shards", type=int, default=2,
                    help="mirror only; how many parquet shards to touch")
    ap.add_argument("--row-groups", type=int, default=1,
                    help="mirror only; row groups per shard (keeps the read small)")
    ap.add_argument("--pick", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="figures/real")
    ap.add_argument("--any-gene", action="store_true",
                    help="do not restrict to the panel germlines")
    ap.add_argument("--allow-duplicate-cdr3", action="store_true",
                    help="permit the same CDRH3 under different V genes "
                         "(off by default: two wells, one molecule)")
    ap.add_argument("--explain", action="store_true",
                    help="print the filter funnel even when rows survive")
    ap.add_argument("--keep-liabilities", action="store_true",
                    help="keep sequences ANARCI flagged (length notes are kept either way)")
    ap.add_argument("--genes", nargs="*", default=None,
                    help="V genes to keep, e.g. --genes IGHV3-23 IGHV3-53. "
                         "Default is the panel; --any-gene keeps everything.")
    ap.add_argument("--align-rows", type=int, default=50,
                    help="sequences per family in the alignment figures")
    ap.add_argument("--aligner", default="reference",
                    choices=["reference", "mafft", "clustalo"],
                    help="reference = built-in germline-anchored aligner (no "
                         "install needed). The others are used if on PATH.")
    ap.add_argument("--imgt-fasta",
                    help="IMGT/GENE-DB AA FASTA on disk; used as the alignment "
                         "reference instead of the OAS germline consensus")
    args = ap.parse_args()

    if args.list_studies:
        from plateforge.library import hf
        for name in hf.list_studies(args.chain):
            print(name)
        return

    if not args.study and not args.source:
        sys.exit("give --study (mirror) or --source (local file or URL); "
                 "--list-studies shows what the mirror has")

    filt = oas.Filter(
        chain="H" if args.chain == "heavy" else "L",
        genes=(None if (args.any_gene or args.chain == "light")
               else (args.genes or germlines.DEFAULT.genes)),
        cdr3_len_range=(8, 30),
        min_redundancy=None,
        exclude_liabilities=not args.keep_liabilities,
    )

    if args.study:
        from plateforge.library import hf
        df, meta = hf.load_study(
            args.study, filt, chain=args.chain,
            species=None if args.species == "any" else args.species,
            limit=args.limit or None, max_shards=args.shards,
            row_groups_per_shard=args.row_groups)
    else:
        df, meta = oas.load(args.source, filt, limit=args.limit or None)

    meta["_source_kind"] = "mirror" if args.study else "local file"
    print(f"\nspecies={meta.get('Species')}  chain={meta.get('Chain')}  "
          f"btype={meta.get('BType')}  isotype={meta.get('Isotype')}")
    print(f"disease={meta.get('Disease')}  author={meta.get('Author')}")
    print(f"scanned {meta['_rows_scanned']:,}   kept {meta['_rows_kept']:,}")

    funnel = meta.get("_funnel")
    if args.explain and funnel is not None:
        print("\nfilter funnel (first chunk):")
        print(funnel.to_string(index=False))

    if not len(df):
        print(f"\nspecies filter dropped {meta.get('_species_dropped', 0):,} rows")
        if funnel is not None:
            print("\nfilter funnel (first chunk) -- the first step that drops "
                  "everything is the culprit:")
            print(funnel.to_string(index=False))
        sys.exit("\nNothing survived. Relax whichever step above emptied it: "
                 "--any-gene, --keep-liabilities, --species any.")

    print("\nANARCI flags in what was kept:")
    counts = oas.flag_counts(df)
    print(counts.head(8).to_string() if len(counts) else "  (none)")

    drift = pool.schema_drift()
    if drift:
        print(f"\nthe stored pool is missing {len(drift)} column(s) this version "
              f"produces: {', '.join(drift[:6])}"
              + (" ..." if len(drift) > 6 else ""))
        print("re-ingesting those rows so they pick the columns up "
              "(seq_id is a content hash, so a plain re-run would skip them).")

    run_id = pool.ingest(df, meta, refresh=bool(drift))
    idx = pool.index()
    print(f"\ningest {run_id}\npool now holds {len(idx):,} sequences")
    print(idx["v_gene"].value_counts().head(12).to_string())

    full = pool.fetch(list(idx["seq_id"]),
                      columns=["seq_id", "v_gene", "j_gene", "aa_gapped", "germline_aa"]
                              + oas.REGION_COLUMNS)
    if args.aligner != "reference" and args.aligner not in msa.available_backends():
        print(f"\n{args.aligner} is not on PATH; using the built-in aligner")
        args.aligner = "reference"

    # The pool is cumulative across ingests, so it holds genes that earlier
    # runs brought in. Sampling the whole index would draw from those too --
    # and the manifest would then record a gene filter the selection visibly
    # violates, which reads as a broken filter rather than as a wide pool.
    # Restrict the selection to the same genes this run asked for, and say how
    # much of the pool that leaves eligible.
    eligible = idx
    if filt.genes:
        eligible = idx[idx["v_gene"].isin(filt.genes)]
        if len(eligible) < len(idx):
            print(f"\nthe pool is cumulative: {len(idx):,} sequences held, "
                  f"{len(eligible):,} match this run's gene filter "
                  f"({', '.join(filt.genes)}). Selecting from those.")
        if eligible.empty:
            sys.exit("no pooled sequences match the requested genes")

    picked = diversity.sample(eligible, args.pick, seed=args.seed,
                              unique_cdr3=not args.allow_duplicate_cdr3)
    print("\n" + picked[["seq_id", "v_gene", "cdr3_aa", "cdr3_len"]].to_string(index=False))
    if len(picked) < args.pick:
        print(f"\nasked for {args.pick}, got {len(picked)}: the pool ran out of "
              "distinct CDRH3s. Widen the filters or lower --pick; padding with "
              "repeats would mean ordering the same molecule twice.")

    report = diversity.Spec(min_genes=1).check(picked)
    print("\ndiversity report")
    for key, val in report.items():
        if isinstance(val, dict):
            print(f"  {'PASS' if val['pass'] else 'FAIL'}  {key}: "
                  f"{val['value']} (limit {val['limit']})")

    lib_id = pool.make_library("oas-real", list(picked["seq_id"]),
                               produced_by="library.diversity",
                               params={"source": meta["_source_ref"], "seed": args.seed},
                               parents={run_id: "sampled_from"})
    print(f"\nlibrary {lib_id}")

    # The alignments describe the sequences that were SELECTED -- these are
    # the ones marched through the downstream steps. Aligning the whole pool
    # would be a picture of the repertoire, not of this order.
    picked_rows = full[full["seq_id"].isin(set(picked["seq_id"]))]

    # Figures first, so the summary panel can go into the run bundle too.
    made = figures.all_figures(idx, picked, args.out)
    gene = picked["v_gene"].value_counts().index[0]
    made.append(figures.aa_legend(f"{args.out}/aa_legend.png"))
    made.append(figures.summary_panel(idx, picked, picked_rows, gene=gene,
                                      out=f"{args.out}/input_summary.png",
                                      max_alignment_rows=args.align_rows))
    # The per-gene alignments also go here, not only into the run bundle.
    # They used to live only here; moving them left an orphan PNG at the old
    # path that no later run overwrote, so it sat there looking current.
    made.extend(figures.alignments_per_gene(picked_rows, args.out,
                                            max_rows=args.align_rows,
                                            backend=args.aligner,
                                            fasta_path=args.imgt_fasta).values())
    bundle = selection.run(
        lib_id=lib_id, pool_index=eligible, picked=picked,
        alignment_rows=picked_rows, meta=meta,
        n_requested=args.pick, seed=args.seed,
        filter_spec=filt, funnel=funnel, diversity_report=report,
        max_alignment_rows=args.align_rows, backend=args.aligner,
        fasta_path=args.imgt_fasta,
        extra_files={"input_summary.png": f"{args.out}/input_summary.png"},
    )

    figures.stamp(args.out, run_id=lib_id, note="selection run")
    orphans = figures.stale(args.out, made)
    if orphans:
        print(f"\n{len(orphans)} figure(s) in {args.out} were NOT written by "
              f"this run and are stale:")
        for p in orphans[:8]:
            print(f"  {p.name}")
        print("  delete them, or `rm -rf figures/` before re-running.")

    print(f"\nrun bundle {bundle.dir}")
    print(f"README block ready to paste: {bundle.dir / 'summary.md'}")
    for f in sorted(p.name for p in bundle.dir.iterdir()):
        print(f"  {f}")
    print("\nfigures:")
    for p in made:
        print(f"  {p}")


if __name__ == "__main__":
    main()
