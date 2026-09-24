"""One selection run, recorded end to end.

The shape of the work is always the same:

    get sequences  ->  filter to a pool  ->  select N  ->  align each family

and the thing that must survive is not any one of those steps but the join
between them. `run()` performs the last step and writes the whole chain into a
`core.runs` bundle: where the sequences came from, what the filter did to them,
how the N were chosen, and one germline-anchored alignment per V gene.

Everything here is a recording of decisions made elsewhere -- the sampler and
the filter are not re-implemented, they are passed in and asked to explain
themselves.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core import runs
from . import figures, msa


def run(*, lib_id: str, pool_index: pd.DataFrame, picked: pd.DataFrame,
        alignment_rows: pd.DataFrame, meta: dict,
        n_requested: int, seed: int,
        filter_spec=None, funnel: pd.DataFrame | None = None,
        diversity_report: dict | None = None,
        genes: list[str] | None = None,
        max_alignment_rows: int = 50,
        backend: str = "reference",
        fasta_path=None,
        extra_files: dict[str, str | Path] | None = None) -> runs.Bundle:
    """Write the bundle for one selection. Returns it, already closed."""
    b = runs.Bundle(lib_id, "sequence_selection", label=meta.get("_source_ref", ""))

    b.section("source", {
        "kind": meta.get("_source_kind", "oas"),
        "ref": meta.get("_source_ref"),
        "study": meta.get("Author"),
        "species": meta.get("Species"),
        "chain": meta.get("Chain"),
        "btype": meta.get("BType"),
        "isotype": meta.get("Isotype"),
        "disease": meta.get("Disease"),
        "scanned": meta.get("_rows_scanned"),
        "kept": meta.get("_rows_kept"),
        "shards": meta.get("_shards"),
    })

    if filter_spec is not None:
        b.section("filter", getattr(filter_spec, "as_dict", lambda: str(filter_spec))())
    if funnel is not None and len(funnel):
        b.write_table("filter_funnel.csv", funnel, role="how the pool was narrowed")
        b.section("filter", {"funnel": funnel.to_dict("records")})

    b.section("pool", {
        "n": int(len(pool_index)),
        "by_gene": {str(k): int(v) for k, v in
                    pool_index["v_gene"].value_counts().items()},
    })

    b.section("selection", {
        "n_requested": int(n_requested),
        "n_selected": int(len(picked)),
        "seed": int(seed),
        "sampler": "library.diversity.sample "
                   "(greedy farthest-point on normalised CDRH3 edit distance)",
        "by_gene": {str(k): int(v) for k, v in
                    picked["v_gene"].value_counts().items()},
    })
    if diversity_report:
        b.section("diversity", diversity_report)

    keep = [c for c in ("seq_id", "v_gene", "j_gene", "cdr3_aa", "cdr3_len",
                        "aa_gapped", "germline_aa") if c in picked.columns]
    b.write_table("selected.csv", picked[keep], role="the sequences this run chose")
    summary = (pool_index.groupby("v_gene")
               .agg(n=("seq_id", "size"),
                    cdr3_len_min=("cdr3_len", "min"),
                    cdr3_len_max=("cdr3_len", "max"),
                    cdr3_len_median=("cdr3_len", "median"))
               .reset_index())
    b.write_table("pool_summary.csv", summary, role="what they were drawn from")

    genes = genes or list(picked["v_gene"].dropna().value_counts().index)
    aligned = {}
    for gene in genes:
        rows = alignment_rows[alignment_rows["v_gene"] == gene]
        if not len(rows):
            continue
        png = figures.alignment(rows, gene, b.dir / f"alignment_{gene}.png",
                                max_rows=max_alignment_rows, backend=backend,
                                fasta_path=fasta_path)
        b.add(png, role=f"{gene} aligned to germline")
        for ext in (".fasta", ".aln.txt"):
            side = png.with_suffix(ext)
            if side.exists():
                b.add(side, role=f"{gene} alignment, {ext.lstrip('.')}")
        m, _ = figures.build_msa(rows, gene, max_rows=max_alignment_rows,
                                 backend=backend, fasta_path=fasta_path)
        if m is not None:
            aligned[gene] = {
                "n_sequences": m.n_sequences,
                "columns": m.width,
                "germline_length": sum(1 for x in m.numbers if x is not None),
                "insertion_columns": sum(1 for x in m.numbers if x is None),
                "reference": m.ref_source,
                "backend": m.backend,
                "median_identity_to_germline": round(
                    float(pd.Series([m.identity(r) for r in m.rows]).median()), 4),
            }
    b.section("alignments", aligned)
    b.section("aligner", {
        "requested": backend,
        "available": msa.available_backends(),
        "numbering": "germline residue index; insertion columns are unnumbered "
                     "and only the reference row is numbered",
    })

    for name, path in (extra_files or {}).items():
        p = Path(path)
        if p.exists():
            (b.dir / name).write_bytes(p.read_bytes())
            b.add(b.dir / name, role="figure")

    b.close()
    # The README block is written last, from the manifest that was just
    # closed, so it can never describe a run other than this one.
    from ..core import runs as _runs
    b.write_text("summary.md", summary_markdown(_runs.read_manifest(lib_id)),
                 role="README block describing this run")
    b.close()
    return b


# --- the run, in a paragraph someone can paste ------------------------------

def summary_markdown(manifest: dict) -> str:
    """A short README block describing how this run's sequences were obtained.

    Generated from the manifest rather than written by hand, because a
    hand-written one goes stale the first time the query changes and nobody
    notices. Every number here came out of the run it describes.
    """
    s = manifest.get("sections", {})
    src, filt = s.get("source", {}), s.get("filter", {})
    sel, aln = s.get("selection", {}), s.get("alignments", {})
    obj = manifest.get("obj_id", "")
    d = f"outputs/{obj}"

    def n(x, dash="—"):
        return f"{int(x):,}" if isinstance(x, (int, float)) else dash

    genes = filt.get("genes") or ["all"]
    lo_hi = filt.get("cdr3_len_range")
    by_gene = ", ".join(f"{g} {v}" for g, v in (sel.get("by_gene") or {}).items())

    lines = [
        "### Input sequences",
        "",
        f"- **Source.** OAS (Observed Antibody Space), unpaired "
        f"{str(src.get('chain') or 'heavy').lower()} chain, "
        f"`{src.get('species', '?')}`, study `{src.get('study') or src.get('ref', '?')}`."
        + (" Read via the HuggingFace parquet mirror — OPIG's own download "
           "paths return 403 and OAS has no API."
           if src.get("kind") == "mirror" else
           f" Read from a local data unit (`{src.get('ref', '?')}`)."),
        f"- **Filter.** {n(src.get('scanned'))} rows scanned → "
        f"{n(src.get('kept'))} kept. Productive only; V gene in "
        f"{', '.join(genes)}; CDRH3 "
        + (f"{lo_hi[0]}–{lo_hi[1]} aa" if lo_hi else "any length")
        + "; ANARCI liabilities dropped (5′ truncation notes are not "
          "liabilities); no ambiguous residues. Per-step funnel in "
          f"`{d}/filter_funnel.csv`.",
        f"- **Selection.** {n(sel.get('n_selected'))} of {n(s.get('pool', {}).get('n'))} "
        f"chosen by greedy farthest-point sampling on normalised CDRH3 edit "
        f"distance (seed {sel.get('seed', '?')}), one CDRH3 per clone — two "
        f"wells of the same molecule is not a screen. "
        + (f"By germline: {by_gene}." if by_gene else ""),
        f"- **Alignment.** Each V gene aligned to its own germline with gaps "
        f"(Needleman–Wunsch, BLOSUM62, affine). Germline is the top row and "
        f"the only numbered one; colour marks departure from it.",
    ]
    if aln:
        for gene, a in aln.items():
            lines.append(
                f"  - {gene}: {a.get('n_sequences')} seqs, {a.get('columns')} "
                f"columns ({a.get('germline_length')} germline + "
                f"{a.get('insertion_columns')} insertion), median identity to "
                f"germline {a.get('median_identity_to_germline')}. "
                f"Reference: {a.get('reference')}.")
    lines += [
        "",
        "| File | Holds |",
        "|---|---|",
        f"| `{d}/manifest.json` | the whole query: source, filter, seed, code version |",
        f"| `{d}/selected.csv` | the {n(sel.get('n_selected'))} sequences, in rank order |",
        f"| `{d}/pool_summary.csv` | what they were drawn from, per germline |",
        f"| `{d}/filter_funnel.csv` | rows surviving each filter step |",
        f"| `{d}/alignment_<gene>.png` `.fasta` `.aln.txt` | one alignment per germline |",
        "",
        f"Produced by code `{manifest.get('code_version', '?')}` on "
        f"{(manifest.get('finished_utc') or '')[:10]}. Re-running mints a new "
        f"`LIB` id and a new directory; nothing here is overwritten.",
        "",
        "### From these sequences to an order",
        "",
        f"`python scripts/order.py --lib {obj}` turns the selection into "
        "constructs and a plate, under `outputs/<CLN set id>/`:",
        "",
        "| Table | One row per | Holds |",
        "|---|---|---|",
        "| `clones.csv` | clone | `clone_id`, `seq_id`, construct format, "
        "linker, protein and DNA sequence, length, GC, Type IIs sites removed, "
        "synthesis warnings, germline and its source, `well`, `plate_barcode` |",
        "| `plate_map.csv` | well | `well` (A01, column-major), `clone_id`, "
        "`cdr3_aa`, `dna_length`, `gc` — the layout as a vendor ships it |",
        "| `plate_map.txt` | — | the same grid, for reading at the bench |",
        "| `order_generic.csv` | well | our own order table: name, sequence, "
        "length, GC, well, vector, cloning strategy |",
        "| `order_<vendor>.csv` | well | the vendor's own layout. Emitted only "
        "with `--allow-unverified` until a real upload template is in "
        "`fixtures/` (Q15), and stamped with a caveat file when it is |",
        "| `vector.json` | — | the insertion context these fragments were "
        "designed for, including the two fusion sites |",
        "",
        "What is ordered is the **VH alone**: the vector already carries the "
        "Kozak, the secretion leader, the (G4S)3 linker, the IGKV1-39/IGKJ1 "
        "light chain and the His6 tag. Every vendor table is derived from "
        "`clones.csv` rather than regenerated, so what is ordered and what is "
        "recorded cannot drift apart.",
    ]
    return "\n".join(lines) + "\n"
