"""Resolve germline reference sequences by gene name.

    python scripts/germlines.py --download                 # fetch IMGT, cache it
    python scripts/germlines.py IGHV3-23 IGKV1-39          # resolve, any source
    python scripts/germlines.py --fasta ~/imgt_aa.fasta IGKV1-39
    python scripts/germlines.py --panel                    # the whole scaffold panel

Two kinds of answer, and the difference matters:

  fasta/url  IMGT/GENE-DB. Authoritative and complete.
  pool       consensus of the germline OAS reports for that gene. No network,
             but covers only genes in the pool and only the span reads covered.

The cache is written under PLATEFORGE_DATA/raw/imgt/ so a download happens once.
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PLATEFORGE_DATA", os.path.expanduser("~/plateforge-data"))

from plateforge.core import paths                                   # noqa: E402
from plateforge.library import germline_db as gdb, germlines        # noqa: E402

CACHE_NAME = "imgt_reference_aa.fasta"


def cache_path():
    return paths.raw_dir("imgt") / CACHE_NAME


def download() -> int:
    import urllib.request
    dest = cache_path()
    print(f"fetching {gdb.IMGT_GENEDB_URL}")
    try:
        req = urllib.request.Request(gdb.IMGT_GENEDB_URL,
                                     headers={"User-Agent": "plateforge"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except Exception as exc:                                        # noqa: BLE001
        print(f"failed: {type(exc).__name__}: {exc}")
        print("\nIMGT may be unreachable from here, as OPIG's own downloads were.")
        print("Download the amino acid reference by hand from")
        print("  https://www.imgt.org/vquest/refseqh.html")
        print(f"and save it to {dest}, then re-run without --download.")
        return 1
    records = gdb.parse_imgt_fasta(text)
    dest.write_text(text)
    print(f"saved {len(text):,} bytes, {len(records):,} genes -> {dest}")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="*", help="gene names, e.g. IGHV3-23")
    ap.add_argument("--download", action="store_true", help="fetch and cache IMGT")
    ap.add_argument("--fasta", help="an IMGT FASTA to use instead of the cache")
    ap.add_argument("--panel", action="store_true",
                    help="resolve the scaffold panel and its light chain")
    ap.add_argument("--show", action="store_true", help="print the sequences")
    args = ap.parse_args()

    if args.download:
        sys.exit(download())

    genes = list(args.genes)
    if args.panel:
        panel = germlines.DEFAULT
        genes += panel.genes + ([panel.light_chain] if panel.light_chain else [])
    if not genes:
        sys.exit("name some genes, or use --panel. --download fetches IMGT first.")

    fasta = args.fasta or (str(cache_path()) if cache_path().exists() else None)
    if not fasta:
        print("no IMGT cache found; falling back to pool-derived germlines.\n"
              "Run --download, or pass --fasta, for authoritative references.\n")

    report = gdb.coverage_report(genes, fasta_path=fasta)
    print(report.to_string(index=False))

    derived = report[report["derived"] == True]                     # noqa: E712
    if len(derived):
        print(f"\n{len(derived)} of {len(report)} are pool-derived consensus, not "
              "IMGT. They cover only the span the reads covered.")
    missing = report[~report["resolved"]]
    if len(missing):
        print(f"\nunresolved: {', '.join(missing['gene'])}")

    if args.show:
        for gene in genes:
            found = gdb.resolve(gene, fasta_path=fasta)
            if found:
                print(f"\n>{found.gene} [{found.source}]")
                print(found.sequence)


if __name__ == "__main__":
    main()
