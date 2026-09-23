"""Turning picked sequences into orderable constructs.

A picked sequence is an observation. A construct is a molecule you can order:
a defined architecture, a definite amino acid sequence, and DNA that a vendor
will actually synthesise. This is where a SEQ becomes a CLN (decisions/0005) --
a clone exists once a sequence is realised in a construct.

Adding an architecture is a new function plus one registration; nothing here
moves. The formats this covers now are scFv and VH-only; VHH and CDRH3
grafting are the same shape of problem and not yet written.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from ..core import artifacts, ids, registry, stores, wells
from . import codon, germline_db, germlines

FORMATS = registry.Registry("construct format")

LINKERS = {
    "G4S3": "GGGGSGGGGSGGGGS",      # 15 residues, the classic scFv linker
    "G4S4": "GGGGSGGGGSGGGGSGGGGS",  # 20 residues, less prone to diabody formation
    "Whitlow": "GSTSGSGKPGSGEGSTKG",  # 218 linker
}


@dataclass
class Construct:
    """One designed molecule, before it is given a clone id."""
    seq_id: str
    fmt: str
    protein: str
    parts: dict[str, str] = field(default_factory=dict)
    notes: str = ""


@FORMATS.register("scfv", needs_light_chain=True)
def scfv(vh: str, vl: str, linker: str = "G4S3", orientation: str = "VH-VL") -> Construct:
    """VH and VL joined by a flexible linker.

    Orientation and linker are explicit because both change expression and
    neither has a universally right answer.
    """
    link = LINKERS.get(linker, linker)
    if orientation == "VH-VL":
        protein = vh + link + vl
        parts = {"first": vh, "linker": link, "second": vl}
    elif orientation == "VL-VH":
        protein = vl + link + vh
        parts = {"first": vl, "linker": link, "second": vh}
    else:
        raise ValueError(f"unknown orientation {orientation!r}")
    return Construct(seq_id="", fmt="scfv", protein=protein, parts=parts,
                     notes=f"{orientation}, {linker} ({len(link)} aa linker)")


@FORMATS.register("vh_only", needs_light_chain=False)
def vh_only(vh: str, vl: str | None = None, **_) -> Construct:
    """The heavy variable domain alone -- a domain antibody, or a VHH scaffold."""
    return Construct(seq_id="", fmt="vh_only", protein=vh, parts={"vh": vh},
                     notes="heavy variable domain only")


def light_chain_for(panel: germlines.Panel | None = None, **resolve_kwargs):
    """Resolve the panel's fixed light chain, keeping its provenance.

    Returns the Germline, so a caller can see whether it came from IMGT or was
    inferred from the pool -- a 5' truncated consensus is not something to
    build 96 orderable constructs on without knowing.
    """
    panel = panel or germlines.DEFAULT
    if not panel.light_chain:
        return None
    return germline_db.resolve(panel.light_chain, **resolve_kwargs)


def build(picked: pd.DataFrame, vl_sequence: str, fmt: str = "scfv",
          linker: str = "G4S3", orientation: str = "VH-VL",
          enzymes: tuple[str, ...] = ("BsaI",),
          sequence_column: str = "aa_seq",
          germline_source: str | None = None,
          optimize_dna: bool = True) -> pd.DataFrame:
    """Design a construct per picked sequence and mint a clone id for each.

    Returns a frame with one row per clone, carrying everything an order form
    and a plate map need. Nothing is written to a store here; `register()`
    does that, so designing and persisting stay separable.
    """
    maker = FORMATS.get(fmt)
    rows = []
    for _, seq in picked.iterrows():
        vh = str(seq[sequence_column])
        built = maker(vh, vl_sequence, linker=linker, orientation=orientation)
        row = {
            "clone_id": ids.mint_stamped("CLN", fmt),
            "seq_id": seq["seq_id"],
            "construct": fmt,
            "v_gene": seq.get("v_gene"),
            "cdr3_aa": seq.get("cdr3_aa"),
            "linker": linker if built.parts.get("linker") else None,
            "orientation": orientation if fmt == "scfv" else None,
            "clone_aa_seq": built.protein,
            "aa_length": len(built.protein),
            "germline": seq.get("v_gene"),
            "germline_source": germline_source,
            "notes": built.notes,
        }
        if optimize_dna:
            opt = codon.optimize(built.protein, enzymes=enzymes)
            row.update({
                "clone_dna_seq": opt.dna,
                "dna_length": len(opt.dna),
                "gc": round(opt.gc, 4),
                "sites_removed": opt.sites_removed,
                "gc_swaps": opt.gc_swaps,
                "synthesis_warnings": "; ".join(opt.warnings) or None,
            })
        rows.append(row)
    return pd.DataFrame(rows)


def to_plate(clones: pd.DataFrame, plate_format: int = 96,
             order: str = "column", reserved: list[str] | None = None,
             barcode: str | None = None) -> pd.DataFrame:
    """Lay clones out on a plate, in the order a vendor ships them.

    Reserved wells are skipped rather than overwritten -- a control well that
    silently receives a clone is the kind of error that is only found at the
    bench.
    """
    reserved = [wells.normalize(w, plate_format) for w in (reserved or [])]
    available = [w for w in wells.all_wells(plate_format, order) if w not in reserved]
    if len(clones) > len(available):
        raise ValueError(
            f"{len(clones)} clones will not fit {len(available)} free wells on a "
            f"{plate_format}-well plate ({len(reserved)} reserved)")

    barcode = barcode or f"PLATE-{ids.stamp()}"
    laid = clones.reset_index(drop=True).copy()
    laid["well"] = available[:len(laid)]
    laid["plate_barcode"] = barcode
    laid["plate_format"] = plate_format
    laid["row"] = [w[0] for w in laid["well"]]
    laid["col"] = [int(w[1:]) for w in laid["well"]]
    return laid


def register(clones: pd.DataFrame, *, produced_by: str = "library.construct",
             parents: dict[str, str] | None = None,
             params: dict | None = None) -> str:
    """Persist clones and record a CLN-set artifact. Returns the artifact id."""
    conn = stores.connect("library")
    cols = [c for c in CLONE_COLUMNS if c in clones.columns]
    rows = clones[cols].copy()
    rows["created_at"] = ids.now_iso()
    rows.to_sql("clones", conn, if_exists="append", index=False)
    conn.commit()

    set_id = ids.mint_stamped("LIB", "clones")
    artifacts.register(
        set_id, "clone_set", produced_by, store="library",
        params=params or {}, parents=parents or {},
        meta={"n": len(clones),
              "construct": clones["construct"].iloc[0] if len(clones) else None},
    )
    conn.executemany(
        "INSERT OR REPLACE INTO library_members (lib_id, seq_id, rank) VALUES (?,?,?)",
        [(set_id, s, i) for i, s in enumerate(clones["seq_id"])])
    conn.commit()
    return set_id


CLONE_COLUMNS = [
    "clone_id", "seq_id", "construct", "v_gene", "cdr3_aa", "linker",
    "orientation", "clone_aa_seq", "aa_length", "clone_dna_seq", "dna_length",
    "gc", "sites_removed", "gc_swaps", "synthesis_warnings", "germline",
    "germline_source", "well", "plate_barcode", "notes",
]
