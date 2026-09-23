"""Construct design: reverse translation, scFv assembly, plate layout."""
import re

import pandas as pd
import pytest

from plateforge.core import stores, wells
from plateforge.library import codon, construct, diversity, oas, pool, synth

VL = ("DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGT"
      "DFTLTISSLQPEDFATYYCQQSYSTPLTFGQGTKVEIK")
VH = ("EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFT"
      "ISRDNSKNTLYLQMNSLRAEDTAVYYCAKARWGYFDYWGQGTLVTVSS")


def _fresh(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()


# --- reverse translation ---------------------------------------------------

def test_round_trip_preserves_the_protein():
    for protein in (VH, VL, VH + "GGGGSGGGGSGGGGS" + VL):
        out = codon.optimize(protein)
        assert codon.translate(out.dna) == protein
        assert len(out.dna) == 3 * len(protein)


def test_bsai_is_removed_from_both_strands():
    out = codon.optimize(VH + VL, enzymes=("BsaI",))
    assert "GGTCTC" not in out.dna
    assert "GAGACC" not in out.dna, "the reverse complement is cut too"


def test_multiple_enzymes_avoided_together():
    out = codon.optimize(VH + VL, enzymes=("BsaI", "BsmBI"))
    for site in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG"):
        assert site not in out.dna, site


def test_sites_to_avoid_includes_reverse_complements():
    sites = codon.sites_to_avoid(("BsaI",))
    assert "GGTCTC" in sites and "GAGACC" in sites


def test_a_literal_site_can_be_passed():
    assert "AAGCTT" in codon.sites_to_avoid(("AAGCTT",))


def test_reverse_complement():
    assert codon.reverse_complement("GGTCTC") == "GAGACC"
    assert codon.reverse_complement("ACGT") == "ACGT"


def test_gc_lands_mid_window_not_on_the_edge():
    out = codon.optimize(VH + "GGGGSGGGGSGGGGS" + VL)
    assert 0.45 <= out.gc <= 0.55, f"GC {out.gc:.3f} should sit near the middle"
    assert not out.warnings


def test_naive_optimisation_would_be_gc_rich():
    """Documents why GC balancing exists: the top human codon is GC-rich."""
    naive = "".join(codon.HUMAN_USAGE[aa][0][0] for aa in VH)
    assert codon.gc_fraction(naive) > 0.65


def test_homopolymer_runs_are_broken():
    # Prolines and threonines together are what produced CCCCCC in real output.
    protein = "MPPPPTTTPPPPKKKKEEEE" + VH
    out = codon.optimize(protein)
    assert codon.homopolymers(out.dna) == []
    assert codon.translate(out.dna) == protein


def test_homopolymer_detection():
    assert codon.homopolymers("AACCCCCCGG") == ["CCCCCC"]
    assert codon.homopolymers("ACGTACGT") == []


def test_unknown_residues_are_reported_not_silently_dropped():
    out = codon.optimize("EVQLXXXLES")
    assert not out.clean
    assert "NNN" in out.dna
    assert any("unknown" in str(u) for u in out.unresolved)


def test_custom_usage_table_is_honoured():
    only = {"A": [("GCG", 1.0)], "M": [("ATG", 1.0)]}
    out = codon.optimize("MAAA", usage=only, gc_window=(0.0, 1.0))
    assert out.dna == "ATGGCGGCGGCG"


# --- scFv assembly ---------------------------------------------------------

def test_scfv_orientation_and_linker():
    made = construct.scfv(VH, VL, linker="G4S3", orientation="VH-VL")
    assert made.protein == VH + construct.LINKERS["G4S3"] + VL
    assert made.protein.startswith(VH) and made.protein.endswith(VL)

    flipped = construct.scfv(VH, VL, orientation="VL-VH")
    assert flipped.protein.startswith(VL) and flipped.protein.endswith(VH)


def test_linker_length_matters():
    short = construct.scfv(VH, VL, linker="G4S3")
    long_ = construct.scfv(VH, VL, linker="G4S4")
    assert len(long_.protein) - len(short.protein) == 5


def test_a_literal_linker_can_be_passed():
    made = construct.scfv(VH, VL, linker="AAAAA")
    assert "AAAAA" in made.protein


def test_unknown_orientation_rejected():
    with pytest.raises(ValueError, match="orientation"):
        construct.scfv(VH, VL, orientation="sideways")


def test_vh_only_needs_no_light_chain():
    made = construct.vh_only(VH, None)
    assert made.protein == VH
    assert not construct.FORMATS.meta("vh_only")["needs_light_chain"]
    assert construct.FORMATS.meta("scfv")["needs_light_chain"]


# --- build and plate -------------------------------------------------------

def _picked(tmp_path, n=12, seed=7):
    src = synth.make_unit(tmp_path / "u.csv.gz", n=4000, seed=seed)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    return diversity.sample(pool.index(), n, seed=1)


def test_build_mints_one_clone_per_sequence(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    picked = _picked(tmp_path)
    clones = construct.build(picked, VL)
    assert len(clones) == len(picked)
    assert clones["clone_id"].nunique() == len(clones)
    assert clones["clone_id"].str.startswith("CLN-").all()
    assert set(clones["seq_id"]) == set(picked["seq_id"])
    stores.close_all()


def test_every_clone_dna_round_trips(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    clones = construct.build(_picked(tmp_path), VL)
    for _, row in clones.iterrows():
        assert codon.translate(row["clone_dna_seq"]) == row["clone_aa_seq"]
    stores.close_all()


def test_no_clone_carries_the_assembly_enzyme_site(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    clones = construct.build(_picked(tmp_path), VL, enzymes=("BsaI",))
    for dna in clones["clone_dna_seq"]:
        assert "GGTCTC" not in dna and "GAGACC" not in dna
        assert not re.search(r"(A{6,}|C{6,}|G{6,}|T{6,})", dna)
    stores.close_all()


def test_plate_layout_is_column_major_and_padded(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    clones = construct.build(_picked(tmp_path, n=12), VL)
    plated = construct.to_plate(clones, 96)
    assert list(plated["well"][:3]) == ["A01", "B01", "C01"]
    assert all(len(w) == 3 for w in plated["well"]), "wells are zero padded"


def test_reserved_wells_are_skipped(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    clones = construct.build(_picked(tmp_path, n=12), VL)
    plated = construct.to_plate(clones, 96, reserved=["A01", "B01"])
    assert "A01" not in set(plated["well"])
    assert "B01" not in set(plated["well"])
    assert plated["well"].iloc[0] == "C01"
    stores.close_all()


def test_overfull_plate_refuses_rather_than_dropping_clones(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    clones = construct.build(_picked(tmp_path, n=96), VL)
    with pytest.raises(ValueError, match="will not fit"):
        construct.to_plate(clones, 96, reserved=["A01", "H12"])
    stores.close_all()


def test_full_plate_fits_exactly(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    clones = construct.build(_picked(tmp_path, n=96), VL)
    plated = construct.to_plate(clones, 96)
    assert len(plated) == 96
    assert plated["well"].nunique() == 96
    assert set(plated["well"]) == set(wells.all_wells(96))
    stores.close_all()


def test_register_persists_and_links(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.core import artifacts
    picked = _picked(tmp_path, n=8)
    lib = pool.make_library("src", list(picked["seq_id"]),
                            produced_by="test")
    plated = construct.to_plate(construct.build(picked, VL), 96)
    set_id = construct.register(plated, parents={lib: "designed_from"})

    conn = stores.connect("library")
    n = conn.execute("SELECT COUNT(*) FROM clones").fetchone()[0]
    assert n == len(plated)
    assert lib in artifacts.lineage(set_id, "up")
    assert artifacts.get(set_id).meta["construct"] == "scfv"
    stores.close_all()


def test_clone_ids_survive_a_second_build(tmp_path, monkeypatch):
    """Clones are minted, not derived, so designing twice gives distinct ids
    rather than silently colliding (decisions/0005)."""
    _fresh(tmp_path, monkeypatch)
    picked = _picked(tmp_path, n=6)
    first = construct.build(picked, VL)
    second = construct.build(picked, VL)
    assert not set(first["clone_id"]) & set(second["clone_id"])
    assert list(first["seq_id"]) == list(second["seq_id"])
    stores.close_all()
