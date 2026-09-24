"""From picked sequences to the DNA you actually put in a cart.

Four separable decisions, deliberately kept separable:

    what molecule   library.construct   scFv, VH-only, ...
    where it goes   library.vector      what the vector already supplies
    how it joins    library.cloning     Golden Gate, Gibson, blunt
    who makes it    library.vendors     the order table's shape and limits

This module only composes them. Every one can be swapped without touching the
others, which is the whole point -- a different cloning strategy is a
different argument here, not a rewrite.

The important distinction it maintains:

    insert      what is ORDERED. With the working vector, the VH alone.
    expressed   what the cell MAKES. leader-VH-(G4S)3-VL-His6, because the
                vector supplies everything but the VH.

Conflating those two is how a plate of fragments gets ordered that cannot
express anything, so both are carried as columns and the assembled reading
frame is translated and checked for every clone before the table is written.
"""
from __future__ import annotations

import json

import pandas as pd

from ..core import ids
from . import cloning, codon, construct, diversity, vector as vectors


def design(picked: pd.DataFrame, *,
           vector_name: str = "pcdna-igg1-dual-vk-v1",
           strategy: str = "golden_gate",
           min_order_length: int = 300,
           enzymes: tuple[str, ...] = ("BsaI",),
           sequence_column: str = "aa_seq",
           fmt: str = "igg1",
           germline_source: str | None = None,
           group_by_similarity: bool = True,
           plate_format: int = 96,
           plate_order: str = "column",
           reserved: list[str] | None = None,
           barcode: str | None = None,
           **strategy_kwargs) -> pd.DataFrame:
    """One row per clone: the molecule, the fragment, the well.

    `group_by_similarity` reorders the picks so neighbouring wells hold
    similar clones. It changes which well a clone lands in and nothing else --
    the set was already chosen by the sampler and is not touched here.
    """
    vec = vectors.get(vector_name)
    five, three = vec.fusion_sites()
    # Where each clone's germline reference came from, resolved once per gene.
    sources: dict[str, str | None] = {}
    if group_by_similarity:
        picked = diversity.group_for_plate(picked)

    rows = []
    for _, seq in picked.iterrows():
        vh = str(seq[sequence_column])
        gene = seq.get("v_gene")
        if gene not in sources:
            from . import germline_db
            found = germline_db.resolve(str(gene)) if gene else None
            sources[gene] = (f"{found.source}:{found.allele}"
                             if found and found.allele else
                             (found.source if found else None))
        opt = codon.optimize(vh, enzymes=enzymes)
        fragment = cloning.design(opt.dna, vec, strategy,
                                  min_length=min_order_length,
                                  **strategy_kwargs)
        orf = vec.orf(opt.dna)
        expressed = codon.translate(orf).rstrip("*")
        chains = {name: len(p.rstrip("*"))
                  for name, p in vec.proteins(opt.dna).items()}
        problems = cloning.verify(fragment, opt.dna, vec, expected_protein=vh)

        rows.append({
            "clone_id": ids.mint_stamped("CLN", fmt),
            "seq_id": seq["seq_id"],
            "construct": fmt,
            "v_gene": seq.get("v_gene"),
            "j_gene": seq.get("j_gene"),
            "cdr3_aa": seq.get("cdr3_aa"),
            "cdr3_len": seq.get("cdr3_len"),
            # what is ordered
            "insert_name": vec.insert,
            "insert_aa_seq": vh,
            "insert_aa_length": len(vh),
            "insert_dna_seq": opt.dna,
            "order_dna_seq": fragment.dna,
            "order_length": fragment.length,
            "order_gc": round(fragment.gc, 4),
            # how it joins
            "vector": vec.name,
            "cloning_strategy": fragment.strategy,
            "enzyme": vec.enzyme if fragment.strategy == "golden_gate" else None,
            "fusion_site_5p": five,
            "fusion_site_3p": three,
            "cloneable": fragment.cloneable,
            # what is expressed -- the clone's own protein and the ORF that
            # encodes it, which is the fragment PLUS everything the vector
            # already carries. clone_* is the molecule; order_* is the parcel.
            "clone_aa_seq": expressed,
            "aa_length": len(expressed),
            "clone_dna_seq": orf,
            "dna_length": len(orf),
            "gc": round(codon.gc_fraction(orf), 4),
            "linker": None,
            "orientation": None,
            "cassettes_json": json.dumps(chains, sort_keys=True),
            "n_chains": len(chains),
            "light_chain": f"{vectors.VL_V_GENE}/{vectors.VL_J_GENE}",
            # how it was made
            "germline": seq.get("v_gene"),
            "germline_source": germline_source or sources.get(gene),
            "gc_swaps": opt.gc_swaps,
            "sites_removed": opt.sites_removed,
            "synthesis_warnings": "; ".join(opt.warnings + fragment.warnings) or None,
            "assembly_problems": "; ".join(problems) or None,
        })

    clones = pd.DataFrame(rows)
    return construct.to_plate(clones, plate_format=plate_format,
                              order=plate_order, reserved=reserved,
                              barcode=barcode)


def report(clones: pd.DataFrame) -> dict:
    """What a human should be told before this becomes an order."""
    def count(column):
        return int(clones[column].notna().sum()) if column in clones else 0

    lengths = clones["order_length"] if "order_length" in clones else pd.Series(dtype=int)
    gc = clones["order_gc"] if "order_gc" in clones else pd.Series(dtype=float)
    return {
        "clones": int(len(clones)),
        "vector": clones["vector"].iloc[0] if len(clones) else None,
        "cloning_strategy": clones["cloning_strategy"].iloc[0] if len(clones) else None,
        "order_bp_min": int(lengths.min()) if len(lengths) else None,
        "order_bp_max": int(lengths.max()) if len(lengths) else None,
        "order_bp_total": int(lengths.sum()) if len(lengths) else None,
        "gc_min": round(float(gc.min()), 3) if len(gc) else None,
        "gc_max": round(float(gc.max()), 3) if len(gc) else None,
        "with_synthesis_warnings": count("synthesis_warnings"),
        "with_assembly_problems": count("assembly_problems"),
        "not_cloneable": int((~clones["cloneable"]).sum())
                         if "cloneable" in clones else 0,
    }


def summary_markdown(clones: pd.DataFrame, vector_name: str | None = None,
                     set_id: str = "", strategy: str | None = None) -> str:
    """A README block for one order: how the DNA was made and verified.

    Generated from the clone table, like the selection block, so the numbers
    in it are the run's own.
    """
    from . import vector as vectors

    vec = vectors.get(vector_name or (clones["vector"].iloc[0] if len(clones)
                                      else "pcdna-igg1-t2a-vk-v1"))
    rep = report(clones)
    d = f"outputs/{set_id}" if set_id else "outputs/<CLN set id>"
    five, three = vec.fusion_sites()
    enzyme = vec.enzyme
    warned = rep["with_synthesis_warnings"]

    single = len(vec.cassettes) == 1
    if single:
        shape = ("Full IgG1, **one ORF** — heavy and light chain split at "
                 "translation by furin–GSG–T2A. One transcript, so the chain "
                 "ratio is a property of 2A skipping and cannot be tuned.")
    else:
        shape = (f"Full IgG1, **{len(vec.cassettes)} cassettes** — each chain "
                 f"on its own promoter and terminator, so the chain ratio is "
                 f"a design parameter. Assembly wants light chain in excess of "
                 f"heavy (~1.5–2:1); free heavy chain aggregates and takes "
                 f"yield with it.")
    lines = [
        "### Constructs and synthesis",
        "",
        f"- **Format.** {shape} Only the **{vec.insert}** is ordered.",
        f"- **Vector.** `{vec.name}` ({vec.backbone}):",
        "",
    ]
    for c in vec.cassettes:
        parts = " → ".join(f"**[{e.name}]**" if e.kind == "insert" else e.name
                           for e in c.elements)
        lines.append(f"      {c.name}: {parts}")
    lines += [
        "",
        f"- **Ligation.** {rep['cloning_strategy']}, {enzyme} (Type IIs). "
        f"Fusion sites `{five}` (5′) / `{three}` (3′), taken from the constant "
        f"flanks — not the VH edges, because VH N-termini differ by germline "
        f"(EVQL vs QVQL) and all 96 must share one vector. Overhangs are "
        f"copies of vector sequence, so the junction is scarless.",
        f"- **Codon optimisation.** `library.codon`: human codon usage, GC "
        f"targeted to the middle of the window (not merely inside it, which "
        f"parks on the boundary), internal {enzyme} sites removed by "
        f"synonymous swap, homopolymer runs >6 broken. Repair is scored by "
        f"severity, so splitting one long run into two shorter ones counts as "
        f"progress.",
        f"- **Result.** {rep['clones']} fragments, {rep['order_bp_min']}–"
        f"{rep['order_bp_max']} bp ({rep['order_bp_total']:,} bp total), GC "
        f"{rep['gc_min']:.0%}–{rep['gc_max']:.0%}. "
        + (f"{warned} carry synthesis warnings." if warned
           else "No synthesis warnings."),
        f"- **Verification.** Not the fragment alone — the **assembled ORF**. "
        f"For every clone the vector's upstream + VH + downstream is "
        f"concatenated and translated, and checked for frame, internal stop "
        f"codons, and that the designed VH is actually present in the product "
        f"({rep['with_assembly_problems']} failures). A fragment that is "
        f"correct on its own and wrong once ligated is only visible at the "
        f"join.",
        f"- **Padding.** Fragments are padded outside the {enzyme} sites to a "
        f"300 bp floor (IDT eBlocks minimum) with fixed filler carrying no "
        f"Type IIs site. The coding sequence is never touched.",
        "",
        "### Where the sequences came from",
        "",
        f"- **VH.** OAS, selected upstream — see the selection run's "
        f"`summary.md`. Germline references: "
        f"{', '.join(sorted(set(str(s) for s in clones['germline_source'].dropna())))}.",
        f"- **VL.** IMGT `{vectors.VL_V_GENE}*01` + `{vectors.VL_JUNCTION}` + "
        f"`{vectors.VL_J_GENE}*01`, from the germline tables bundled with "
        f"`anarci`. `{vectors.VL_JUNCTION}` is the CDR-L3 junction: it is **not** "
        f"germline-encoded by either gene — a V gene stops at the conserved "
        f"Cys+2 and a J gene starts mid-CDR3 — so it is stored as its own "
        f"constant rather than hidden inside a 107-aa string. It is the "
        f"germline-configured IGKV1-39/IGKJ1 CDR-L3 used in most synthetic "
        f"libraries; swap it for a new vector entry if your library differs.",
        f"- **Constant regions.** UniProt `{vectors.IGHG1_ACCESSION}` (IGHG1) "
        f"and `{vectors.IGKC_ACCESSION}` (IGKC). The records are stored "
        f"unmodified in `fixtures/uniprot/`; the one derivation — P01857 is the "
        f"membrane-bound isoform, so the secreted C-terminus is truncated at "
        f"`KSLSLSP` + `GK` — is re-derived from the file by a test, not "
        f"hardcoded.",
        f"- **Leader.** Mouse Ig-kappa secretion leader, on both chains. "
        f"Without it the antibody is not secreted.",
        "",
        "| File | Holds |",
        "|---|---|",
        f"| `{d}/clones.csv` | one row per clone, {len(clones.columns)} columns: "
        "`clone_*` is the molecule, `order_*` is the parcel |",
        f"| `{d}/plate_map.csv` `.txt` | 96-well layout, A01–H12 column-major, "
        "neighbouring wells grouped by similarity |",
        f"| `{d}/order_generic.csv` | our own order table |",
        f"| `{d}/order_<vendor>.csv` | the vendor's layout (see its CAVEAT file) |",
        f"| `{d}/vector.json` | the insertion context, including both fusion sites |",
        f"| `{d}/metadata_table.png` | the clone record, as a figure |",
        f"| `{d}/alignment_<gene>.png` | the clones that will be synthesised |",
    ]
    return "\n".join(lines) + "\n"
