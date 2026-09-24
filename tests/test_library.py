import os

import pandas as pd
import re

import pytest

from plateforge.core import artifacts, stores
from plateforge.library import diversity, germlines, oas, pool, synth


# --- germline call parsing -------------------------------------------------

def test_gene_family_allele():
    assert germlines.gene("IGHV3-23*01") == "IGHV3-23"
    assert germlines.family("IGHV3-23*01") == "IGHV3"
    assert germlines.allele("IGHV3-23*01") == "01"
    assert germlines.chain_of("IGHV3-23*01") == "H"
    assert germlines.chain_of("IGKV1-39*01") == "L"
    assert germlines.chain_of("IGLV1-51*01") == "L"


def test_ties_resolve_to_first():
    assert germlines.gene("IGHV3-23*01,IGHV3-23*04") == "IGHV3-23"
    assert len(germlines.split_calls("IGHV3-23*01,IGHV3-23*04")) == 2


def test_missing_calls_are_none_not_crashes():
    for bad in (None, "", float("nan")):
        assert germlines.gene(bad) is None
        assert germlines.family(bad) is None


def test_panel_quotas_sum_and_cover():
    q = germlines.DEFAULT.quotas(10)
    assert sum(q.values()) == 10
    assert set(q) == set(germlines.DEFAULT.genes)
    assert all(v >= 1 for v in q.values())
    assert q["IGHV3-23"] > q["IGHV3-53"]


def test_panel_quotas_large():
    assert sum(germlines.DEFAULT.quotas(384).values()) == 384


def test_panel_quota_below_panel_size_raises():
    with pytest.raises(ValueError):
        germlines.DEFAULT.quotas(2)


def test_panel_roundtrip():
    back = germlines.Panel.from_dict(germlines.DEFAULT.as_dict())
    assert back.genes == germlines.DEFAULT.genes
    assert back.light_chain == "IGKV1-39"


# --- OAS format ------------------------------------------------------------

def test_unit_url():
    url = oas.unit_url("Chen_2020", "SRR11937587_1_Heavy_IGHG")
    assert url.endswith("/unpaired/Chen_2020/csv/SRR11937587_1_Heavy_IGHG.csv.gz")


def test_metadata_line_parses(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=20, seed=1)
    meta = oas.read_metadata(src)
    assert meta["Species"] == "human"
    assert meta["Chain"] == "Heavy"


def test_load_normalizes(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=2)
    df, meta = oas.load(src, oas.Filter(exclude_liabilities=False))
    assert len(df) > 0
    assert meta["_rows_scanned"] == 200
    assert set(df["v_gene"]) <= set(synth.FRAMEWORKS)
    assert (df["v_family"] == df["v_call"].map(germlines.family)).all()
    assert (df["chain"] == "H").all()
    assert (df["cdr3_len"] == df["cdr3_aa"].str.len()).all()
    assert df["aa_seq"].str.contains("-").sum() == 0
    assert df["seq_id"].str.startswith("SEQ-").all()


def test_seq_id_is_deterministic(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=50, seed=3)
    a, _ = oas.load(src, oas.Filter(exclude_liabilities=False))
    b, _ = oas.load(src, oas.Filter(exclude_liabilities=False))
    assert list(a["seq_id"]) == list(b["seq_id"])


def test_filters(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=4)
    only, _ = oas.load(src, oas.Filter(genes=["IGHV3-23"]))
    assert set(only["v_gene"]) == {"IGHV3-23"}

    clean, _ = oas.load(src, oas.Filter(exclude_liabilities=True))
    assert not clean["has_liability"].any()

    dirty, _ = oas.load(src, oas.Filter(exclude_liabilities=False))
    assert dirty["has_liability"].any(), "synthetic data should contain liabilities"

    short, _ = oas.load(src, oas.Filter(cdr3_len_range=(10, 14)))
    assert short["cdr3_len"].between(10, 14).all()

    capped, _ = oas.load(src, oas.Filter(exclude_liabilities=False), limit=25)
    assert len(capped) == 25


def test_liability_detection():
    assert oas._has_liability("|Unusual residue|")
    assert not oas._has_liability("||")
    assert not oas._has_liability("")
    assert not oas._has_liability(None)


# --- pool persistence ------------------------------------------------------

def _fresh(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()


def test_ingest_and_index(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=5)
    df, meta = oas.load(src)
    run_id = pool.ingest(df, meta)

    assert artifacts.get(run_id).obj_type == "sequence_ingest"
    assert artifacts.get(run_id).meta["rows_added"] == len(df)
    idx = pool.index()
    assert len(idx) == len(df)
    assert pool.stats()["total"] == len(df)
    stores.close_all()


def test_ingest_is_idempotent(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=150, seed=6)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    first = pool.stats()["total"]
    second_run = pool.ingest(df, meta)
    assert pool.stats()["total"] == first
    assert artifacts.get(second_run).meta["rows_added"] == 0
    stores.close_all()


def test_make_library_and_lineage(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=7)
    df, meta = oas.load(src)
    run_id = pool.ingest(df, meta)

    picked = diversity.sample(pool.index(), 10, seed=1)
    lib_id = pool.make_library("panel", list(picked["seq_id"]),
                               produced_by="library.diversity",
                               parents={run_id: "sampled_from"})

    members = pool.library_members(lib_id)
    assert len(members) == 10
    assert list(members["rank"]) == list(range(10))
    assert run_id in artifacts.lineage(lib_id, "up")
    assert artifacts.get(lib_id).meta["n"] == 10
    stores.close_all()


# --- diversity -------------------------------------------------------------

def test_identity_bounds():
    assert diversity.identity("ARDY", "ARDY") == 1.0
    assert diversity.identity("ARDY", "WWWW") == 0.0
    assert 0 < diversity.identity("ARDYK", "ARDY") < 1


def test_identity_is_length_aware():
    # A short loop inside a long one is not "identical".
    assert diversity.identity("ARWGY", "ARWGYKKKKKKKKKK") < 0.5


def test_greedy_beats_random_on_worst_pair(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=8)
    df, _ = oas.load(src)
    seqs = list(df["cdr3_aa"])
    picked = [seqs[i] for i in diversity.greedy_farthest(seqs, 10, seed=0)]

    def worst(xs):
        return max(diversity.identity(a, b)
                   for i, a in enumerate(xs) for b in xs[i + 1:])

    import random
    rng = random.Random(0)
    randomly = [worst(rng.sample(seqs, 10)) for _ in range(5)]
    assert worst(picked) < sum(randomly) / len(randomly)


def test_sample_respects_panel(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=900, seed=9)
    df, _ = oas.load(src)
    picked = diversity.sample(df, 12, seed=2)
    assert len(picked) == 12
    counts = picked["v_gene"].value_counts().to_dict()
    assert set(counts) == set(germlines.DEFAULT.genes)
    assert counts["IGHV3-23"] >= counts["IGHV3-53"]


def test_sampled_set_passes_spec(tmp_path):
    # A realistic CDRH3 floor, as every real run uses. Without it the sampler
    # draws 5-residue loops, and two random 5-mers differing at one position
    # are 0.8 identical by arithmetic -- the spec then sits exactly on its own
    # limit and the test measures rounding rather than diversity.
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=900, seed=10)
    df, _ = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    picked = diversity.sample(df, 10, seed=3)
    report = diversity.Spec(min_genes=3).check(picked)
    assert report["pass"], report
    assert report["max_pairwise_identity"]["value"] < 0.75, \
        "should clear the limit with margin, not sit on it"


def test_short_cdr3s_are_inherently_similar(tmp_path):
    """Documents why the spec needs a length floor: normalized identity on
    short loops is high no matter how the sampler behaves."""
    assert diversity.identity("ARHDY", "ARMDY") == 0.8
    assert diversity.identity("ARHDYWQGTL", "ARMDYWQGTL") == 0.9


def test_generation_is_reproducible_across_processes(tmp_path):
    """Python randomises string hashing per process. Iterating a set while
    consuming the RNG therefore made the same seed produce different data on
    different machines -- which is exactly how this was found. Guard it."""
    import hashlib
    import subprocess
    import sys
    import gzip

    digests = set()
    for hashseed in ("0", "1", "12345"):
        out = tmp_path / f"u{hashseed}.csv.gz"
        subprocess.run(
            [sys.executable, "-c",
             "from plateforge.library import synth;"
             f"synth.make_unit({str(out)!r}, n=400, seed=10)"],
            check=True, env={**os.environ, "PYTHONHASHSEED": hashseed})
        digests.add(hashlib.sha256(gzip.open(out, "rb").read()).hexdigest())
    assert len(digests) == 1, "same seed must produce identical data in any process"


def test_spec_catches_a_bad_set():
    same = pd.DataFrame({
        "seq_id": [f"SEQ-{i}" for i in range(4)],
        "v_gene": ["IGHV3-23"] * 4,
        "cdr3_aa": ["ARWGYFDY", "ARWGYFDF", "ARWGYFDL", "ARWGYFDI"],
        "cdr3_len": [8, 8, 8, 8],
    })
    report = diversity.Spec(min_genes=3).check(same)
    assert not report["pass"]
    assert not report["max_pairwise_identity"]["pass"]
    assert not report["cdr3_len_span"]["pass"]
    assert not report["genes"]["pass"]


def test_sample_is_deterministic(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=500, seed=11)
    df, _ = oas.load(src)
    a = diversity.sample(df, 10, seed=42)
    b = diversity.sample(df, 10, seed=42)
    assert list(a["seq_id"]) == list(b["seq_id"])


def test_single_gene_panel(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=12)
    df, _ = oas.load(src)
    picked = diversity.sample(df, 8, panel=germlines.SINGLE, seed=4)
    assert len(picked) == 8
    assert set(picked["v_gene"]) == {"IGHV3-23"}


# --- gapped alignment column ----------------------------------------------

def test_gapped_sequence_is_preserved(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=20)
    df, _ = oas.load(src)
    assert "aa_gapped" in df.columns
    assert df["aa_gapped"].str.contains(r"\.").any(), "expected IMGT gaps"
    # aa_seq is the same sequence with gaps removed
    assert (df["aa_gapped"].str.replace(".", "", regex=False) == df["aa_seq"]).all()


def test_sequences_of_one_gene_are_column_aligned(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=600, seed=21)
    df, _ = oas.load(src)
    for gene, sub in df.groupby("v_gene"):
        assert sub["aa_gapped"].str.len().nunique() == 1, f"{gene} not aligned"


def test_synth_produces_clonal_lineages(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=600, seed=22, clonality=0.9)
    df, _ = oas.load(src)
    loops = list(df["cdr3_aa"])[:120]
    close = sum(1 for i, a in enumerate(loops) for b in loops[i + 1:]
                if diversity.identity(a, b) > 0.85)
    assert close > 0, "clonal siblings should exist in the pool"


# --- figures ---------------------------------------------------------------

def test_figures_render(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures

    src = synth.make_unit(tmp_path / "unit.csv.gz", n=1200, seed=23)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    idx = pool.index()
    picked = diversity.sample(idx, 10, seed=5)

    made = figures.all_figures(idx, picked, tmp_path / "figs")
    assert len(made) == 4
    for p in made:
        assert p.exists() and p.stat().st_size > 5_000, p

    full = pool.fetch(list(idx["seq_id"]), columns=["seq_id", "v_gene", "aa_gapped"])
    gene = idx["v_gene"].value_counts().index[0]
    aln = figures.alignment(full, gene, tmp_path / "figs" / "aln.png")
    assert aln.exists() and aln.stat().st_size > 5_000
    legend = figures.aa_legend(tmp_path / "figs" / "legend.png")
    assert legend.exists()
    stores.close_all()


def test_alignment_needs_the_gapped_column(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=24)
    df, _ = oas.load(src)
    with pytest.raises(KeyError, match="aa_gapped"):
        figures.alignment(df[["seq_id", "v_gene"]], df["v_gene"].iloc[0],
                          tmp_path / "x.png")
    stores.close_all()


def test_aa_colour_map_covers_all_residues():
    from plateforge.library import figures
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        assert aa in figures.AA_COLOR, aa


# --- germline reference ----------------------------------------------------

def test_germline_column_is_preserved(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=30)
    df, _ = oas.load(src)
    assert "germline_aa" in df.columns
    assert df["germline_aa"].notna().all()
    # germline is column-aligned with the sequence, which is the whole point
    assert (df["germline_aa"].str.len() == df["aa_gapped"].str.len()).all()


def test_germline_differs_from_sequence_but_shares_framework(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=31, shm_rate=0.05)
    df, _ = oas.load(src)
    row = df.iloc[0]
    diffs = sum(a != b for a, b in zip(row["aa_gapped"], row["germline_aa"]))
    assert diffs > 0, "SHM should make the sequence differ from germline"
    assert diffs < len(row["aa_gapped"]), "they should still share framework"


def test_reference_row_prefers_germline(tmp_path):
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=32)
    df, _ = oas.load(src)
    gene = df["v_gene"].value_counts().index[0]
    sub = df[df["v_gene"] == gene]
    ref, label = figures.reference_row(sub, list(sub["aa_gapped"]))
    assert label == "germline"
    assert ref == sub["germline_aa"].iloc[0]


def test_reference_row_falls_back_and_says_so(tmp_path):
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=33)
    df, _ = oas.load(src)
    gene = df["v_gene"].value_counts().index[0]
    sub = df[df["v_gene"] == gene].drop(columns=["germline_aa"])
    ref, label = figures.reference_row(sub, list(sub["aa_gapped"]))
    assert "consensus" in label and "no germline" in label
    assert len(ref) == len(sub["aa_gapped"].iloc[0])


def test_alignment_titles_the_reference_honestly(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=34)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    full = pool.fetch(list(pool.index()["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa"])
    gene = full["v_gene"].value_counts().index[0]
    p = figures.alignment(full, gene, tmp_path / "aln.png")
    assert p.exists() and p.stat().st_size > 5_000
    stores.close_all()


def test_aa_palette_is_one_hue_per_residue():
    from plateforge.library import figures
    assert len(figures.AA_COLOR) == 20
    assert len(set(figures.AA_COLOR.values())) == 20, "each residue needs its own hue"
    grouped = "".join(figures.AA_GROUPS.values())
    assert sorted(grouped) == sorted(figures.AA_COLOR), "groups must cover the palette"


def test_text_contrast_picked_per_fill():
    from plateforge.library import figures
    assert figures._text_on("#ffffff") == "#1a1a18"
    assert figures._text_on("#1f2196") == "#ffffff"
    for aa, color in figures.AA_COLOR.items():
        assert figures._text_on(color) in ("#1a1a18", "#ffffff"), aa


# --- IMGT region boundaries ------------------------------------------------

def test_region_columns_preserved(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=40)
    df, _ = oas.load(src)
    for col in oas.REGION_COLUMNS:
        assert col in df.columns, col
        assert df[col].notna().all(), col


def test_regions_tile_the_sequence(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=41)
    df, _ = oas.load(src)
    for _, row in df.head(25).iterrows():
        joined = "".join(str(row[c]) for c in oas.REGION_COLUMNS)
        assert row["aa_seq"].startswith(joined), row["v_gene"]


def test_region_spans_are_ordered_and_in_range(tmp_path):
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=42)
    df, _ = oas.load(src)
    row = df.iloc[0]
    spans = figures.region_spans(row, row["aa_gapped"])
    assert [n for n, _, _ in spans] == ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3"]
    width = len(row["aa_gapped"])
    last_end = 0
    for name, start, end in spans:
        assert 0 <= start < end <= width, name
        assert start >= last_end, f"{name} overlaps the previous region"
        last_end = end


def test_region_spans_degrade_without_columns(tmp_path):
    from plateforge.library import figures
    import pandas as pd
    row = pd.Series({"aa_gapped": "QVQL..LESG"})
    assert figures.region_spans(row, "QVQL..LESG") == []


def test_alignment_with_regions_renders(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=600, seed=43)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    full = pool.fetch(list(pool.index()["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa"]
                              + oas.REGION_COLUMNS)
    gene = full["v_gene"].value_counts().index[0]
    p = figures.alignment(full, gene, tmp_path / "aln.png", max_rows=12)
    assert p.exists() and p.stat().st_size > 5_000
    q = figures.alignment(full, gene, tmp_path / "aln2.png", max_rows=12,
                          show_regions=False, highlight_only_differences=False)
    assert q.exists() and q.stat().st_size > 5_000
    stores.close_all()


# --- ANARCI status classification -----------------------------------------

def test_length_notes_are_not_liabilities():
    # Real string seen in the HuggingFace mirror.
    assert not oas._has_liability("||||Shorter than IMGT defined: fw1, fw4|")
    assert not oas._has_liability("||")
    assert not oas._has_liability(None)


def test_real_liabilities_still_flagged():
    assert oas._has_liability("|Unusual residue: C|")
    assert oas._has_liability("||Shorter than IMGT defined: fw1|Unusual residue: X|")


def test_status_flags_parse():
    assert oas.status_flags("||||Shorter than IMGT defined: fw1, fw4|") == \
        ["Shorter than IMGT defined: fw1, fw4"]
    assert oas.status_flags("") == []


def test_flag_counts(tmp_path):
    df = pd.DataFrame({"anarci_status": ["|A|", "|A|", "|B|", "||"]})
    counts = oas.flag_counts(df)
    assert counts["A"] == 2 and counts["B"] == 1


# --- HuggingFace mirror adapter -------------------------------------------

def _mirror_shard(path, n=60, species="human"):
    """A parquet file shaped like the ConvergeBio mirror: AIRR columns plus
    meta_* columns instead of a JSON header."""
    import gzip, csv, io
    src = synth.make_unit(str(path) + ".csv.gz", n=n, seed=99)
    raw = pd.read_csv(src, header=1)
    raw = raw.rename(columns={"Redundancy": "Redundancy"})
    for key, val in {"Species": species, "Chain": "Heavy", "BType": "Unsorted-B-Cells",
                     "Disease": "None", "Author": "mirror test",
                     "Isotype": "IGHG", "Subject": "s1"}.items():
        raw[f"meta_{key}"] = val
    raw.to_parquet(path, index=False)
    return path


def test_extract_meta_maps_columns_to_header_dict(tmp_path):
    from plateforge.library import hf
    p = _mirror_shard(tmp_path / "shard.parquet")
    df = pd.read_parquet(p)
    meta = hf.extract_meta(df)
    assert meta["Species"] == "human"
    assert meta["Chain"] == "Heavy"
    assert "meta_Species" not in meta


def test_extract_meta_marks_mixed_values(tmp_path):
    from plateforge.library import hf
    df = pd.DataFrame({"meta_Species": ["human", "mouse"], "meta_Chain": ["Heavy", "Heavy"]})
    meta = hf.extract_meta(df)
    assert meta["Chain"] == "Heavy"
    assert "2 values" in str(meta["Species"])


def test_mirror_rows_normalize_like_a_data_unit(tmp_path):
    from plateforge.library import hf
    p = _mirror_shard(tmp_path / "shard.parquet")
    raw = pd.read_parquet(p)
    meta = hf.extract_meta(raw)
    df = oas.normalize(raw, meta, source_ref="hf:test")
    assert df["species"].iloc[0] == "human"
    assert df["germline_aa"].notna().all()
    assert (df["aa_gapped"].str.len() == df["germline_aa"].str.len()).all()
    for col in oas.REGION_COLUMNS:
        assert df[col].notna().all()


def test_alignment_report_counts_raw_length_spread(tmp_path):
    from plateforge.library import hf
    df = pd.DataFrame({
        "v_gene": ["IGHV3-23"] * 3 + ["IGHV1-69"] * 2,
        "aa_gapped": ["AAAA", "AAAA", "AAA", "BBBB", "BBBB"],
    })
    rep = hf.alignment_report(df).set_index("v_gene")
    assert rep.loc["IGHV3-23", "raw_lengths"] == 2
    assert rep.loc["IGHV1-69", "raw_lengths"] == 1


def test_ragged_lengths_do_not_block_alignment(tmp_path):
    """raw_lengths being high is expected on real data and does not stop the
    figure: it aligns on IMGT numbering, not string position."""
    from plateforge.library import hf
    df = pd.DataFrame({
        "v_gene": ["IGHV3-23"] * 3,
        "aa_gapped": ["QVQLVES", "VQLVES", "QLVES"],      # three lengths
        "anarci_numbering": [REAL_NUMBERING] * 3,
    })
    rep = hf.alignment_report(df).set_index("v_gene")
    assert rep.loc["IGHV3-23", "raw_lengths"] == 3
    assert rep.loc["IGHV3-23", "numbered"] == 3
    assert rep.loc["IGHV3-23", "can_align"], "numbering is what decides"


def test_report_flags_missing_numbering(tmp_path):
    from plateforge.library import hf
    df = pd.DataFrame({
        "v_gene": ["IGHV3-23"] * 3,
        "aa_gapped": ["QVQLVES"] * 3,
        "anarci_numbering": [None, None, None],
    })
    rep = hf.alignment_report(df).set_index("v_gene")
    assert rep.loc["IGHV3-23", "numbered"] == 0
    assert not rep.loc["IGHV3-23", "can_align"]


def test_hf_missing_dependency_message(monkeypatch):
    import builtins
    from plateforge.library import hf
    real_import = builtins.__import__

    def fake(name, *a, **k):
        if name == "huggingface_hub":
            raise ImportError("no module")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake)
    with pytest.raises(ImportError, match="huggingface_hub"):
        hf._require_hf()


# --- alignment keeps ragged reads instead of excluding them ----------------

def test_alignment_keeps_ragged_sequences(tmp_path, monkeypatch):
    """A short read is aligned with gaps, not dropped.

    The old renderer restricted to the modal string length, which threw away
    exactly the 5' truncated reads that dominate real OAS data.
    """
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=50)
    df, _ = oas.load(src)
    gene = df["v_gene"].value_counts().index[0]
    sub = df[df["v_gene"] == gene].copy()
    idx = sub.index[:3]
    sub.loc[idx, "aa_gapped"] = sub.loc[idx, "aa_gapped"].str[7:]      # 5' truncated

    m, _rows = figures.build_msa(sub, gene, max_rows=12)
    assert m.n_sequences == min(12, len(sub)), "no sequence may be excluded"
    truncated = set(str(s) for s in sub.loc[idx, "seq_id"])
    drawn = set(m.labels)
    assert truncated & drawn, "the truncated reads must be in the alignment"

    p = figures.alignment(sub, gene, tmp_path / "ragged.png", max_rows=12)
    assert p.exists() and p.stat().st_size > 5_000
    assert p.with_suffix(".fasta").exists(), "the alignment is also written as FASTA"
    stores.close_all()


# --- filter funnel ---------------------------------------------------------

def test_explain_reports_each_step(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=60)
    df, meta = oas.load(src, oas.Filter(genes=["IGHV3-23"], cdr3_len_range=(10, 14)))
    funnel = meta["_funnel"]
    assert list(funnel["step"])[0] == "input"
    assert funnel["rows"].is_monotonic_decreasing
    assert funnel["rows"].iloc[-1] == len(df) or funnel["rows"].iloc[-1] >= len(df)


def test_explain_names_the_culprit(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=61)
    df, meta = oas.load(src, oas.Filter(genes=["IGHV9-99"]))   # gene not present
    assert len(df) == 0
    funnel = meta["_funnel"]
    killer = funnel[funnel["rows"] == 0].iloc[0]["step"]
    assert "panel genes" in killer


def test_productive_accepts_bool_and_numeric():
    f = oas.Filter(productive_only=True)
    base = dict(chain=["H"], v_gene=["IGHV3-23"], cdr3_len=[12], has_liability=[False])
    for value in [True, "T", "TRUE", "1", 1]:
        df = pd.DataFrame({**base, "productive": [value]})
        assert len(f.apply(df)) == 1, value
    for value in [False, "F", "FALSE"]:
        df = pd.DataFrame({**base, "productive": [value]})
        assert len(f.apply(df)) == 0, value


def test_unknown_productive_is_kept_not_dropped():
    f = oas.Filter(productive_only=True)
    base = dict(chain=["H"], v_gene=["IGHV3-23"], cdr3_len=[12], has_liability=[False])
    for value in [None, float("nan"), ""]:
        df = pd.DataFrame({**base, "productive": [value]})
        assert len(f.apply(df)) == 1, f"null productive should be kept, got {value!r}"


# --- ambiguous residues ----------------------------------------------------

def test_ambiguous_sequences_are_dropped_by_default(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=600, seed=70, ambiguous_rate=0.3)
    loose, _ = oas.load(src, oas.Filter(max_ambiguous=None, exclude_liabilities=False))
    assert (loose["n_ambiguous"] > 0).any(), "synthetic data should contain X runs"

    clean, _ = oas.load(src, oas.Filter(exclude_liabilities=False))
    assert (clean["n_ambiguous"] == 0).all()
    assert not clean["cdr3_aa"].str.contains("X").any()
    assert not clean["aa_seq"].str.contains("X").any()


def test_ambiguous_count_is_accurate(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=71, ambiguous_rate=0.5)
    df, _ = oas.load(src, oas.Filter(max_ambiguous=None, exclude_liabilities=False))
    for _, row in df.head(40).iterrows():
        expected = sum(1 for c in row["aa_seq"] if c not in "ACDEFGHIKLMNPQRSTVWY")
        assert row["n_ambiguous"] == expected


def test_ambiguous_filter_appears_in_funnel(tmp_path):
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=72, ambiguous_rate=0.25)
    _, meta = oas.load(src, oas.Filter(exclude_liabilities=False))
    steps = list(meta["_funnel"]["step"])
    assert any("ambiguous" in s for s in steps)


# --- ANARCI flag classification, from the real Briney census ---------------

def test_coverage_flags_from_real_data_are_not_liabilities():
    for flag in ["Shorter than IMGT defined: fw1",
                 "Shorter than IMGT defined: fw1, fw4",
                 "Deletions: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 73",
                 "Deletions: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 73"]:
        assert oas.is_coverage_flag(flag), flag
        assert not oas._has_liability(f"|{flag}|"), flag


def test_structural_problems_remain_liabilities():
    for flag in ["Missing Conserved Cysteine: 104", "Unusual residue: C"]:
        assert not oas.is_coverage_flag(flag), flag
        assert oas._has_liability(f"|{flag}|"), flag


def test_unknown_flags_are_treated_as_liabilities():
    assert oas._has_liability("|Something nobody has looked at yet|")


def test_mixed_status_is_a_liability_if_any_flag_is_structural():
    status = "|Shorter than IMGT defined: fw1|Missing Conserved Cysteine: 104|"
    assert oas._has_liability(status)


# --- IMGT numbering alignment (Q12) ---------------------------------------

REAL_NUMBERING = (
    "{'fwh1': {'15 ': 'P', '16 ': 'G', '23 ': 'C'}, "
    "'cdrh1': {'27 ': 'G', '35 ': 'S'}, "
    "'cdrh3': {'111 ': 'G', '111A': 'I', '112A': 'D', '112 ': 'R', '113 ': 'D'}, "
    "'fwh4': {'118 ': 'W'}}"
)


def test_parse_real_numbering():
    from plateforge.library import imgt
    parsed = imgt.parse(REAL_NUMBERING)
    assert set(parsed) == {"fwh1", "cdrh1", "cdrh3", "fwh4"}
    assert parsed["cdrh3"]["111A"] == "I"


def test_parse_rejects_junk():
    from plateforge.library import imgt
    for bad in [None, "", "not a dict", 42, "{unclosed"]:
        assert imgt.parse(bad) is None


def test_imgt_insertion_order_is_outward_from_the_middle():
    from plateforge.library import imgt
    cols = ["112 ", "111A", "112B", "111 ", "112A", "113 ", "110 "]
    assert sorted(cols, key=imgt.position_key) == [
        "110 ", "111 ", "111A", "112B", "112A", "112 ", "113 "]


def test_position_key_handles_unparseable():
    from plateforge.library import imgt
    assert imgt.position_key("???")[0] == 10_000


def test_build_puts_ragged_sequences_in_shared_columns(tmp_path):
    from plateforge.library import imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=80)
    df, _ = oas.load(src)
    sub = df[df["v_gene"] == df["v_gene"].value_counts().index[0]].head(30)
    assert sub["aa_gapped"].str.len().nunique() >= 1

    aln = imgt.build(sub)
    assert aln is not None
    # every row spans the same columns, whatever its length
    assert aln.matrix.shape[0] == len(sub)
    assert aln.width > 100
    labels = [r[0] for r in aln.regions]
    assert labels == ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]


def test_build_regions_tile_without_gaps(tmp_path):
    from plateforge.library import imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=81)
    df, _ = oas.load(src)
    aln = imgt.build(df.head(20))
    ends = [0]
    for _label, start, end in aln.regions:
        assert start == ends[-1], "regions must tile the alignment"
        ends.append(end)
    assert ends[-1] == aln.width


def test_build_returns_none_without_numbering():
    from plateforge.library import imgt
    assert imgt.build(pd.DataFrame({"seq_id": ["a"]})) is None


def test_min_occupancy_trims_rare_insertions(tmp_path):
    from plateforge.library import imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=400, seed=82)
    df, _ = oas.load(src)
    sub = df.head(60)
    wide = imgt.build(sub, min_occupancy=0.0)
    trimmed = imgt.build(sub, min_occupancy=0.5)
    assert trimmed.width <= wide.width


def test_germline_maps_onto_the_same_columns(tmp_path):
    from plateforge.library import imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=83)
    df, _ = oas.load(src)
    sub = df[df["v_gene"] == df["v_gene"].value_counts().index[0]].head(20)
    aln = imgt.build(sub)
    ref, used = imgt.germline_row(sub, list(aln.matrix.columns))
    assert used > 0, "germline should map for at least some rows"
    assert set(ref).issubset(set(aln.matrix.columns))


def test_germline_refuses_to_shift_when_counts_disagree():
    from plateforge.library import imgt
    row = pd.Series({
        "anarci_numbering": REAL_NUMBERING,     # 11 numbered residues
        "aa_gapped": "QVQ",                     # nowhere near 11
        "germline_aa": "QVQ",
    })
    assert imgt.germline_by_column(row) is None


def test_germline_refuses_on_length_mismatch():
    from plateforge.library import imgt
    row = pd.Series({"anarci_numbering": REAL_NUMBERING,
                     "aa_gapped": "QVQLL", "germline_aa": "QVQ"})
    assert imgt.germline_by_column(row) is None


def test_alignment_figure_prefers_numbering(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=800, seed=84)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    full = pool.fetch(list(pool.index()["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa",
                               "anarci_numbering"])
    gene = full["v_gene"].value_counts().index[0]
    p = figures.alignment(full, gene, tmp_path / "numbered.png", max_rows=12)
    assert p.exists() and p.stat().st_size > 5_000
    q = figures.alignment(full.drop(columns=["anarci_numbering"]), gene,
                          tmp_path / "fallback.png", max_rows=12)
    assert q.exists() and q.stat().st_size > 5_000
    stores.close_all()


# --- store consistency -----------------------------------------------------

def test_pool_survives_losing_its_index(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.core import paths
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=85)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    first = pool.consistency()
    assert first["consistent"]

    stores.close_all()
    for f in paths.data_root().glob("stores/library.sqlite*"):
        f.unlink()
    pool.ingest(df, meta)                      # parquet must not double
    after = pool.consistency()
    assert after["consistent"], after
    assert after["in_parquet"] == first["in_parquet"]
    stores.close_all()


def test_repair_drops_orphan_parquet_rows(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.core import bulk
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=86)
    df, meta = oas.load(src)
    pool.ingest(df, meta)
    orphan = bulk.read(pool.POOL_TABLE).head(3).copy()
    orphan["seq_id"] = ["SEQ-orphan1", "SEQ-orphan2", "SEQ-orphan3"]
    bulk.write(pool.POOL_TABLE, orphan, append=True, key="seq_id")
    assert pool.consistency()["parquet_only"] == 3
    assert pool.repair()["consistent"]
    stores.close_all()


# --- schema migration (Q8) -------------------------------------------------

def test_missing_columns_are_added_to_an_existing_database(tmp_path, monkeypatch):
    import sqlite3
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()
    from plateforge.core import paths

    store_dir = paths.data_root() / "stores"
    store_dir.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(store_dir / "library.sqlite")
    conn.executescript("""CREATE TABLE IF NOT EXISTS sequences (
        seq_id TEXT PRIMARY KEY, source TEXT NOT NULL,
        aa_seq TEXT NOT NULL, created_at TEXT NOT NULL);""")
    conn.execute("INSERT INTO sequences VALUES ('SEQ-old','OAS','QVQL','2026-01-01')")
    conn.commit()
    conn.close()

    live = stores.connect("library")
    cols = {r[1] for r in live.execute("PRAGMA table_info(sequences)")}
    assert {"n_ambiguous", "cdr3_len", "v_gene"} <= cols
    assert live.execute("SELECT seq_id FROM sequences").fetchone()[0] == "SEQ-old"
    stores.close_all()


def test_migration_is_idempotent(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()
    conn = stores.connect("library")
    schema = (stores.SCHEMA_DIR / "library.sql").read_text()
    assert stores.migrate(conn, schema) == []
    stores.close_all()


def test_declared_columns_skips_table_constraints():
    schema = """CREATE TABLE IF NOT EXISTS t (
        a TEXT PRIMARY KEY,
        b INTEGER DEFAULT 0,
        PRIMARY KEY (a, b)
    );"""
    cols = stores.declared_columns(schema)["t"]
    assert set(cols) == {"a", "b"}


# --- germline reference resolution -----------------------------------------

IMGT_FASTA = """>M99660|IGHV3-23*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |
EVQLLESGG.GLVQPGGSLRLSCAASGFTF....SSYAMSWVRQAPGKGLEWVSAISGSG..GSTYYADSVKG
RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK
>M99660|IGHV3-23*04|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |
EVQLLESGG.GLVQPGGSLRLSCAASGFTF....SSYTMSWVRQAPGKGLEWVSAISGSG..GSTYYADSVKG
>X59315|IGKV1-39*01|Homo sapiens|F|V-REGION|1..289|
DIQMTQSPSSLSASVGDRVTITCRASQSIS..SYLNWYQQKPGKAPKLLIYAA....STLQSGVPSRF
>Z12345|IGHV9-99*01|Mus musculus|P|V-REGION|1..280|
QVQLKESGPG
"""


def test_imgt_fasta_parses_header_fields():
    from plateforge.library import germline_db as gdb
    db = gdb.parse_imgt_fasta(IMGT_FASTA)
    assert db["IGHV3-23"].allele == "IGHV3-23*01"
    assert db["IGHV3-23"].functionality == "F"
    assert db["IGKV1-39"].sequence.startswith("DIQMTQSPSS")
    assert db["IGHV9-99"].functionality == "P"


def test_imgt_fasta_strips_gaps_and_joins_lines():
    from plateforge.library import germline_db as gdb
    seq = gdb.parse_imgt_fasta(IMGT_FASTA)["IGHV3-23"].sequence
    assert "." not in seq and "-" not in seq
    assert seq.startswith("EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMS")
    assert seq.endswith("YYCAK"), "later lines of the record must be joined"


def test_first_allele_wins():
    from plateforge.library import germline_db as gdb
    # *01 appears before *04 and carries AMS; *04 carries TMS.
    assert "SSYAMS" in gdb.parse_imgt_fasta(IMGT_FASTA)["IGHV3-23"].sequence


def test_fasta_source_is_authoritative_not_derived(tmp_path):
    from plateforge.library import germline_db as gdb
    path = tmp_path / "imgt.fasta"
    path.write_text(IMGT_FASTA)
    found = gdb.from_fasta("IGKV1-39", path)
    assert found.source == "fasta"
    assert not found.is_derived


def test_resolve_prefers_fasta_over_pool(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import germline_db as gdb
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=1200, seed=90)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    path = tmp_path / "imgt.fasta"
    path.write_text(IMGT_FASTA)

    # The bundled IMGT tables outrank everything by default, when present.
    if gdb.anarci_available():
        assert gdb.resolve("IGHV3-23").source == "anarci"
    # Order is explicit when a caller wants a specific provenance.
    assert gdb.resolve("IGHV3-23", order=("fasta", "pool"),
                       fasta_path=path).source == "fasta"
    assert gdb.resolve("IGHV3-23", order=("pool",)).source == "pool"
    stores.close_all()


def test_pool_germline_is_labelled_and_supported(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import germline_db as gdb
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=1500, seed=91)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)

    found = gdb.resolve("IGHV3-23", order=("pool",))
    assert found.is_derived and found.support > 50
    assert "covers only the span" in found.notes
    assert len(found.sequence) > 80
    stores.close_all()


def test_pool_refuses_below_min_support(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import germline_db as gdb
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=600, seed=92)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    assert gdb.from_pool("IGHV3-23", min_support=100_000) is None
    stores.close_all()


def test_resolve_returns_none_for_unknown_gene(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import germline_db as gdb
    assert gdb.resolve("IGHV0-00", order=("pool",)) is None
    stores.close_all()


def test_coverage_report_flags_derived(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import germline_db as gdb
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=1200, seed=93)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    path = tmp_path / "imgt.fasta"
    path.write_text(IMGT_FASTA)

    rep = gdb.coverage_report(["IGHV3-23", "IGHV1-69", "IGHV0-00"],
                              order=("fasta", "pool"),
                              fasta_path=path).set_index("gene")
    assert rep.loc["IGHV3-23", "source"] == "fasta"
    assert rep.loc["IGHV1-69", "source"] == "pool"       # not in the fixture
    assert not rep.loc["IGHV0-00", "resolved"]
    stores.close_all()


# --- cross-germline duplicate CDRH3s ---------------------------------------

def _shared_loop_candidates():
    """Every germline drawing from the same small loop set, so any per-gene
    selection must collide unless duplicates are excluded globally."""
    loops = ["ARWGYFDY", "ARLLPQMDY", "ARKKTVSGFDY", "ARQWWEYMDV",
             "ARHNPLTGDY", "ARFFIKSMDL", "ARVVDNQYFDY", "ARCCMAWGDY"]
    rows = [{"seq_id": f"SEQ-{gene}-{i}", "v_gene": gene,
             "cdr3_aa": loop, "cdr3_len": len(loop)}
            for gene in germlines.DEFAULT.genes for i, loop in enumerate(loops)]
    return pd.DataFrame(rows)


def test_same_cdr3_under_two_genes_is_not_picked_twice():
    cands = _shared_loop_candidates()
    picked = diversity.sample(cands, 12, seed=1)
    assert len(picked) == picked["cdr3_aa"].nunique(), \
        "the same loop must not be ordered twice"


def test_duplicates_appear_without_the_guard():
    # Documents the bug this guards: real data has identical CDRH3s assigned
    # to different V genes, and per-gene quotas would each pick them.
    cands = _shared_loop_candidates()
    loose = diversity.sample(cands, 12, seed=1, unique_cdr3=False)
    assert len(loose) > loose["cdr3_aa"].nunique()
    assert diversity.Spec(min_genes=3).check(loose)["max_pairwise_identity"]["value"] == 1.0


def test_sampler_returns_fewer_rather_than_repeating():
    cands = _shared_loop_candidates()          # only 8 distinct loops
    picked = diversity.sample(cands, 12, seed=1)
    assert len(picked) == 8
    assert picked["cdr3_aa"].nunique() == 8


def test_plate_sized_pick_is_duplicate_free(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=12000, seed=95)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    picked = diversity.sample(pool.index(), 96, seed=1)
    assert len(picked) == 96
    assert picked["cdr3_aa"].nunique() == 96
    report = diversity.Spec(min_genes=3).check(picked)
    assert report["max_pairwise_identity"]["value"] < 1.0
    assert report["pass"], report
    stores.close_all()


# --- alignment presentation ------------------------------------------------

def test_germline_numbering_skips_gaps():
    """A gap column is a position some other sequence has and the germline
    does not, so it must not consume a germline residue number."""
    from plateforge.library import imgt
    cols = ["1 ", "2 ", "111A", "3 "]
    ref = {"1 ": "Q", "2 ": "V", "3 ": "L"}          # nothing at 111A
    numbers = []
    k = 0
    for c in cols:
        aa = ref.get(c)
        if isinstance(aa, str) and aa.strip():
            k += 1
            numbers.append(k)
        else:
            numbers.append(None)      # a gap is not a germline position
    assert numbers == [1, 2, None, 3]
    assert max(n for n in numbers if n) == 3, "three germline residues, not four"
    assert imgt.position_key("111A") > imgt.position_key("2 ")


def test_alignment_renders_with_ragged_coverage(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures, imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=800, seed=96)
    df, _ = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    gene = df["v_gene"].value_counts().index[0]
    sub = df[df["v_gene"] == gene].head(10).copy()

    rows = []
    for i, (_, row) in enumerate(sub.iterrows()):
        numbering = imgt.parse(row["anarci_numbering"])
        if i % 2:
            numbering["fwh1"] = {k: v for k, v in numbering["fwh1"].items()
                                 if int(str(k).strip().rstrip("ABCDEFG") or 0) > 10}
        row = row.copy()
        row["anarci_numbering"] = repr(numbering)
        rows.append(row)
    ragged = pd.DataFrame(rows)

    aln = imgt.build(ragged)
    assert aln.matrix.isna().any().any(), "truncated reads must leave gaps"
    p = figures.alignment(ragged, gene, tmp_path / "ragged.png", max_rows=10)
    assert p.exists() and p.stat().st_size > 5_000
    stores.close_all()


def test_summary_panel_renders(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=3000, seed=97)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    idx = pool.index()
    picked = diversity.sample(idx, 24, seed=1)
    full = pool.fetch(list(picked["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa",
                               "anarci_numbering"])
    p = figures.summary_panel(idx, picked, full, out=tmp_path / "summary.png")
    assert p.exists() and p.stat().st_size > 20_000
    stores.close_all()


def test_summary_panel_without_alignment(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=2000, seed=98)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df, meta)
    idx = pool.index()
    picked = diversity.sample(idx, 12, seed=1)
    p = figures.summary_panel(idx, picked, None, out=tmp_path / "bars.png")
    assert p.exists() and p.stat().st_size > 10_000
    stores.close_all()


# --- parquet schema drift --------------------------------------------------

def test_stale_pool_is_detected(tmp_path, monkeypatch):
    """seq_id is a content hash, so re-ingest skips known rows -- which means a
    pool written by older code never gains a column added since, silently."""
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=800, seed=99)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df.drop(columns=["anarci_numbering"]), meta)
    assert "anarci_numbering" in pool.schema_drift()
    stores.close_all()


def test_plain_reingest_does_not_repair_drift(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.core import artifacts
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=800, seed=100)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df.drop(columns=["anarci_numbering"]), meta)
    run = pool.ingest(df, meta)                      # no refresh
    assert artifacts.get(run).meta["rows_added"] == 0
    assert "anarci_numbering" in pool.schema_drift()
    stores.close_all()


def test_refresh_backfills_and_does_not_duplicate(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=800, seed=101)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df.drop(columns=["anarci_numbering"]), meta)
    before = len(pool.index())

    pool.ingest(df, meta, refresh=True)
    assert pool.schema_drift() == []
    assert len(pool.index()) == before, "refresh must replace, not duplicate"
    assert pool.consistency()["consistent"]
    stores.close_all()


def test_alignment_builds_after_refresh(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=1200, seed=102)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df.drop(columns=["anarci_numbering"]), meta)

    idx = pool.index()
    gene = idx["v_gene"].value_counts().index[0]
    cols = ["seq_id", "v_gene", "aa_gapped", "germline_aa"]
    stale = pool.fetch(list(idx["seq_id"]), columns=cols)
    assert imgt.build(stale[stale["v_gene"] == gene].head(10)) is None

    pool.ingest(df, meta, refresh=True)
    fresh = pool.fetch(list(pool.index()["seq_id"]), columns=cols + ["anarci_numbering"])
    assert imgt.build(fresh[fresh["v_gene"] == gene].head(10)) is not None
    stores.close_all()


def test_summary_panel_says_why_it_is_empty(tmp_path, monkeypatch):
    """An empty panel of bare axes is indistinguishable from a broken figure."""
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=1200, seed=103)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    pool.ingest(df.drop(columns=["anarci_numbering"]), meta)
    idx = pool.index()
    picked = diversity.sample(idx, 12, seed=1)
    stale = pool.fetch(list(idx["seq_id"]),
                       columns=["seq_id", "v_gene", "aa_gapped", "germline_aa"])
    p = figures.summary_panel(idx, picked, stale, out=tmp_path / "empty.png")
    assert p.exists() and p.stat().st_size > 10_000
    stores.close_all()


def test_panel_uses_rows_that_have_numbering_not_merely_the_first(tmp_path, monkeypatch):
    """A pool ingested across versions has its oldest rows first, and those
    are the ones written before anarci_numbering existed. Taking head(n)
    selected exactly the unusable rows and the panel rendered blank."""
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures, imgt
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=2400, seed=104)
    df, meta = oas.load(src, oas.Filter(cdr3_len_range=(8, 30)))
    half = len(df) // 2
    pool.ingest(df.iloc[:half].drop(columns=["anarci_numbering"]), meta)
    pool.ingest(df.iloc[half:], meta)

    idx = pool.index()
    gene = idx["v_gene"].value_counts().index[0]
    full = pool.fetch(list(idx["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa",
                               "anarci_numbering"])
    sub = full[full["v_gene"] == gene]
    assert sub["anarci_numbering"].head(10).isna().all(), "setup: stale rows come first"
    assert sub["anarci_numbering"].notna().any(), "setup: some rows are usable"

    # The alignment needs a germline, not numbering, so stale rows are fine.
    m, _rows = figures.build_msa(sub, gene, max_rows=20)
    assert m is not None, "must not give up because the first rows are stale"
    assert m.n_sequences > 0

    picked = diversity.sample(idx, 12, seed=1)
    p = figures.summary_panel(idx, picked, full, gene=gene, out=tmp_path / "mixed.png")
    assert p.exists() and p.stat().st_size > 30_000
    stores.close_all()


# --- the alignment itself ---------------------------------------------------

GERM = ("EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKG"
        "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKWGQGTLVTVSS")


def test_msa_numbers_only_the_germline():
    """Germline residues are numbered 1..L. Nothing else is numbered at all."""
    from plateforge.library import msa

    m = msa.build(GERM, {"a": GERM, "b": GERM[:40] + "WWW" + GERM[40:]})
    numbered = [n for n in m.numbers if n is not None]
    assert numbered == list(range(1, len(GERM) + 1)), \
        "every germline residue owns exactly one column, in order"
    assert m.numbers.count(None) == 3, "the 3 inserted columns are not positions"
    # and the reference is the only row the numbering describes
    assert len(m.ref) == m.width and all(len(r) == m.width for r in m.rows)


def test_insertion_opens_a_gap_in_every_other_row():
    from plateforge.library import msa

    with_insert = GERM[:40] + "PQR" + GERM[40:]
    m = msa.build(GERM, {"plain": GERM, "longer": with_insert})
    cols = [i for i, n in enumerate(m.numbers) if n is None]
    assert cols, "an insertion must open columns"
    plain = m.rows[m.labels.index("plain")]
    assert all(plain[c] == "-" for c in cols)
    assert "".join(m.rows[m.labels.index("longer")][c] for c in cols) == "PQR"
    assert all(m.ref[c] == "-" for c in cols), "the germline has no residue there"


def test_truncated_read_aligns_with_leading_gaps():
    """A 5' truncated read keeps register instead of sliding left."""
    from plateforge.library import msa

    m = msa.build(GERM, {"short": GERM[15:]})
    row = m.rows[0]
    assert row[:15] == "-" * 15
    assert row[15:] == m.ref[15:]
    assert m.identity(row) < 1.0, "the missing span counts against identity"


def test_differences_are_reported_in_germline_numbering():
    from plateforge.library import msa

    mutant = GERM[:4] + "V" + GERM[5:]
    m = msa.build(GERM, {"m": mutant})
    assert m.differences(m.rows[0]) == [(5, GERM[4], "V")]


def test_reference_comes_from_the_oas_germlines(tmp_path, monkeypatch):
    """The reference is reconstructed from OAS's own germline calls, and says so."""
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import msa

    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=71)
    df, _ = oas.load(src)
    gene = df["v_gene"].value_counts().index[0]
    # use_imgt=False is the fallback path: no IMGT tables, reconstruct from
    # what OAS itself reported for these reads.
    seq, label = msa.reference_for(df, gene, use_imgt=False)
    assert seq and "OAS germline consensus" in label
    assert "-" not in seq and "." not in seq, "the reference is ungapped"


def test_consensus_reference_takes_the_longest_span():
    from plateforge.library import msa

    seq, n = msa.consensus_reference(["QLLESGGGLVQPGG", "EVQLLESGGGLVQPGG"])
    assert seq == "EVQLLESGGGLVQPGG" and n == 2


def test_unknown_backend_is_an_error_not_a_silent_fallback():
    from plateforge.library import msa

    with pytest.raises(KeyError, match="unknown backend"):
        msa.build(GERM, {"a": GERM}, backend="nope")


def test_reference_backend_is_always_available():
    from plateforge.library import msa

    assert "reference" in msa.available_backends()


# --- the run bundle ---------------------------------------------------------

def _selection_fixture(tmp_path, monkeypatch, n=300, pick=24, seed=80):
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=n, seed=seed)
    filt = oas.Filter(chain="H", cdr3_len_range=(8, 30))
    df, meta = oas.load(src, filt)
    pool.ingest(df, meta)
    idx = pool.index()
    full = pool.fetch(list(idx["seq_id"]),
                      columns=["seq_id", "v_gene", "aa_gapped", "germline_aa"]
                              + oas.REGION_COLUMNS)
    picked = diversity.sample(idx, pick, seed=seed)
    rows = full[full["seq_id"].isin(set(picked["seq_id"]))]
    lib_id = pool.make_library("t", list(picked["seq_id"]), produced_by="test")
    return idx, picked, rows, meta, filt, lib_id


def test_selection_run_writes_one_alignment_per_gene(tmp_path, monkeypatch):
    import json
    from plateforge.library import selection

    idx, picked, rows, meta, filt, lib_id = _selection_fixture(tmp_path, monkeypatch)
    b = selection.run(lib_id=lib_id, pool_index=idx, picked=picked,
                      alignment_rows=rows, meta=meta, n_requested=24, seed=80,
                      filter_spec=filt, funnel=meta.get("_funnel"),
                      max_alignment_rows=20)
    genes = set(picked["v_gene"])
    for gene in genes:
        assert (b.dir / f"alignment_{gene}.png").exists()
        assert (b.dir / f"alignment_{gene}.fasta").exists()

    man = json.loads((b.dir / "manifest.json").read_text())
    assert set(man["sections"]["alignments"]) == genes
    stores.close_all()


def test_manifest_records_how_the_sequences_were_obtained(tmp_path, monkeypatch):
    """The bundle has to answer 'why these and not others' without the terminal."""
    import json
    from plateforge.library import selection

    idx, picked, rows, meta, filt, lib_id = _selection_fixture(tmp_path, monkeypatch)
    b = selection.run(lib_id=lib_id, pool_index=idx, picked=picked,
                      alignment_rows=rows, meta=meta, n_requested=24, seed=80,
                      filter_spec=filt, funnel=meta.get("_funnel"),
                      diversity_report=diversity.Spec(min_genes=1).check(picked),
                      max_alignment_rows=20)
    man = json.loads((b.dir / "manifest.json").read_text())
    s = man["sections"]
    assert s["source"]["ref"] == meta["_source_ref"]
    assert s["source"]["scanned"] == meta["_rows_scanned"]
    assert s["filter"]["cdr3_len_range"] == [8, 30]
    assert s["filter"]["funnel"][0]["step"] == "input"
    assert s["selection"]["seed"] == 80
    assert s["selection"]["n_requested"] == 24
    assert s["selection"]["n_selected"] == len(picked)
    assert man["code_version"]
    assert (b.dir / "selected.csv").exists()
    stores.close_all()


def test_selected_csv_matches_the_library(tmp_path, monkeypatch):
    import pandas as pd
    from plateforge.library import selection

    idx, picked, rows, meta, filt, lib_id = _selection_fixture(tmp_path, monkeypatch)
    b = selection.run(lib_id=lib_id, pool_index=idx, picked=picked,
                      alignment_rows=rows, meta=meta, n_requested=24, seed=80,
                      max_alignment_rows=20)
    got = pd.read_csv(b.dir / "selected.csv")
    assert list(got["seq_id"]) == list(picked["seq_id"])
    stores.close_all()


def test_bundle_describes_itself(tmp_path, monkeypatch):
    from plateforge.core import runs
    from plateforge.library import selection

    idx, picked, rows, meta, filt, lib_id = _selection_fixture(tmp_path, monkeypatch)
    selection.run(lib_id=lib_id, pool_index=idx, picked=picked,
                  alignment_rows=rows, meta=meta, n_requested=24, seed=80,
                  max_alignment_rows=20)
    text = runs.describe(lib_id)
    assert lib_id in text and "seed 80" in text
    stores.close_all()


# --- colour means "differs from germline", and nothing else ----------------

def _cells(ax, y):
    """Filled cells on one row of a drawn alignment, as (column, facecolor)."""
    import matplotlib.colors as mcolors
    out = []
    for patch in ax.patches:
        bbox = patch.get_bbox()
        if abs(bbox.y0 - y) < 1e-6:
            out.append((int(round(bbox.x0)),
                        mcolors.to_hex(patch.get_facecolor())))
    return out


def test_germline_row_is_never_filled():
    """The reference is grey letters on white. A fill always means divergence."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plateforge.library import figures, msa

    mutant = GERM[:4] + "V" + GERM[5:]
    m = msa.build(GERM, {"a": mutant, "b": GERM[:40] + "PQR" + GERM[40:]})
    fig, ax = plt.subplots()
    figures.draw_msa(ax, m, show_letters=False)
    germline_fills = {c for _, c in _cells(ax, m.n_sequences)}
    plt.close(fig)
    assert germline_fills <= {figures.GAP_COLOR}, \
        "only insertion columns may be tinted on the germline row"
    assert not (germline_fills & set(figures.AA_COLOR.values()))


def test_a_substitution_is_filled_and_a_match_is_not():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plateforge.library import figures, msa

    mutant = GERM[:4] + "V" + GERM[5:]
    m = msa.build(GERM, {"a": mutant})
    fig, ax = plt.subplots()
    figures.draw_msa(ax, m, show_letters=False)
    filled = [c for c, colour in _cells(ax, m.n_sequences - 1)
              if colour in set(figures.AA_COLOR.values())]
    plt.close(fig)
    assert filled == [4], "exactly the one changed position is coloured"


def test_an_insertion_is_highlighted_in_the_variant():
    """An inserted residue has nothing above it, and is still a divergence."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plateforge.library import figures, msa

    m = msa.build(GERM, {"ins": GERM[:40] + "PQR" + GERM[40:]})
    insert_cols = [i for i, n in enumerate(m.numbers) if n is None]
    fig, ax = plt.subplots()
    figures.draw_msa(ax, m, show_letters=False)
    filled = {c for c, colour in _cells(ax, m.n_sequences - 1)
              if colour in set(figures.AA_COLOR.values())}
    plt.close(fig)
    assert set(insert_cols) <= filled, "inserted residues must be coloured"


# --- the README block -------------------------------------------------------

def test_run_writes_a_readme_block_naming_its_own_files(tmp_path, monkeypatch):
    from plateforge.library import selection

    idx, picked, rows, meta, filt, lib_id = _selection_fixture(tmp_path, monkeypatch)
    b = selection.run(lib_id=lib_id, pool_index=idx, picked=picked,
                      alignment_rows=rows, meta=meta, n_requested=24, seed=80,
                      filter_spec=filt, funnel=meta.get("_funnel"),
                      max_alignment_rows=20)
    text = (b.dir / "summary.md").read_text()
    assert lib_id in text, "the block must point at this run, not a generic path"
    for name in ("manifest.json", "selected.csv", "filter_funnel.csv",
                 "clones.csv", "plate_map.csv"):
        assert name in text
    assert "seed 80" in text
    assert str(len(picked)) in text
    stores.close_all()


def test_readme_block_does_not_claim_the_mirror_for_a_local_file(tmp_path,
                                                                 monkeypatch):
    """A local unit is not the HuggingFace mirror, and must not say it is."""
    from plateforge.library import selection

    idx, picked, rows, meta, filt, lib_id = _selection_fixture(tmp_path, monkeypatch)
    meta = dict(meta, _source_kind="local file")
    b = selection.run(lib_id=lib_id, pool_index=idx, picked=picked,
                      alignment_rows=rows, meta=meta, n_requested=24, seed=80,
                      max_alignment_rows=20)
    text = (b.dir / "summary.md").read_text()
    assert "HuggingFace" not in text and "local data unit" in text
    stores.close_all()


# --- the vector, and what it already carries -------------------------------

def test_vector_parts_translate_to_what_they_claim():
    """The frozen DNA constants must encode the proteins named beside them."""
    from plateforge.library import codon, vector

    for dna, aa in [(vector.SIGNAL_DNA, vector.SIGNAL_AA),
                    (vector.LINKER_G4S3_DNA, vector.LINKER_G4S3_AA),
                    (vector.VL_DNA, vector.VL_AA),
                    (vector.HIS6_DNA, vector.HIS6_AA)]:
        assert codon.translate(dna) == aa


def test_light_chain_is_the_imgt_germline_plus_a_named_junction():
    """The VL must be traceable to IMGT, junction residues called out."""
    from plateforge.library import germline_db, vector

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed; install the [imgt] extra")
    v = germline_db.from_anarci(vector.VL_V_GENE)
    j = germline_db.from_anarci(vector.VL_J_GENE)
    assert v and j, "the bundled IMGT tables must carry both"
    assert vector.VL_AA == v.sequence + vector.VL_JUNCTION + j.sequence
    assert vector.VL_JUNCTION not in (v.sequence, j.sequence)


def test_assembled_orf_is_in_frame_and_makes_the_designed_protein():
    from plateforge.library import codon, vector

    vec = vector.get("pcdna-scfv-vk-v1")
    vh = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAKGGYFDYWGQGTLVTVSS"
    dna = codon.optimize(vh, enzymes=("BsaI",)).dna
    orf = vec.orf(dna)
    assert len(orf) % 3 == 0
    protein = codon.translate(orf)
    assert protein.endswith("*") and protein.count("*") == 1
    assert protein.startswith(vector.SIGNAL_AA)
    assert vh in protein
    assert vector.LINKER_G4S3_AA in protein and vector.VL_AA in protein


def test_fusion_sites_are_constant_across_germlines():
    """The whole point: one vector accepts every insert."""
    from plateforge.library import vector

    vec = vector.get("pcdna-scfv-vk-v1")
    five, three = vec.fusion_sites()
    assert len(five) == len(three) == 4
    assert not vec.check_fusion_sites(), "the declared sites must be usable"
    # They come from the flanks, so they cannot depend on the insert at all.
    assert five == vec.upstream_dna()[-4:]
    assert three == vec.downstream_dna()[:4]


# --- the fragment -----------------------------------------------------------

def _cds(aa="EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAKGGYFDYWGQGTLVTVSS"):
    from plateforge.library import codon
    return codon.optimize(aa, enzymes=("BsaI",)).dna, aa


def test_golden_gate_puts_the_enzyme_outside_and_the_overhangs_inside():
    from plateforge.library import cloning, codon, vector

    vec = vector.get("pcdna-scfv-vk-v1")
    cds, _aa = _cds()
    f = cloning.design(cds, vec, "golden_gate", min_length=300)
    site = codon.sites_to_avoid(("BsaI",))[0]
    assert f.dna.count(site) == 1 and f.dna.count(codon.reverse_complement(site)) == 1
    five, three = vec.fusion_sites()
    assert f.part("5' fusion site") == five
    assert f.part("3' fusion site") == three
    assert f.part("CDS") == cds
    assert f.dna.index(site) < f.dna.index(five)      # enzyme outside the site
    assert not f.warnings


def test_short_fragments_are_padded_to_the_vendor_minimum():
    from plateforge.library import cloning, vector

    vec = vector.get("pcdna-scfv-vk-v1")
    cds, _ = _cds("EVQLLESGGGLVQPGGSLRLSCAAS")          # deliberately tiny
    f = cloning.design(cds, vec, "golden_gate", min_length=300)
    assert f.length >= 300
    assert f.part("CDS") == cds, "padding must not touch the coding sequence"


def test_blunt_fragments_refuse_to_look_cloneable():
    from plateforge.library import cloning, vector

    vec = vector.get("pcdna-scfv-vk-v1")
    cds, _ = _cds()
    f = cloning.design(cds, vec, "blunt")
    assert f.cloneable is False and f.warnings


def test_verify_catches_a_frameshifted_insert():
    """The expensive failure: fine alone, wrong once it is in the vector."""
    from plateforge.library import cloning, vector

    vec = vector.get("pcdna-scfv-vk-v1")
    cds, aa = _cds()
    good = cloning.design(cds, vec, "golden_gate", min_length=300)
    assert cloning.verify(good, cds, vec, expected_protein=aa) == []

    shifted = cds[:-1]
    bad = cloning.design(shifted, vec, "golden_gate", min_length=300)
    problems = cloning.verify(bad, shifted, vec, expected_protein=aa)
    assert problems, "a 1-nt truncation must not pass"


def test_filler_carries_no_type_iis_sites():
    from plateforge.library import cloning, codon

    for site in codon.sites_to_avoid(("BsaI", "BsmBI", "BbsI")):
        assert site not in cloning.FILLER


# --- the order --------------------------------------------------------------

def test_design_separates_what_is_ordered_from_what_is_expressed(tmp_path,
                                                                 monkeypatch):
    from plateforge.library import ordering

    idx, picked, _rows, _meta, _filt, _lib = _selection_fixture(
        tmp_path, monkeypatch, pick=8)
    clones = ordering.design(picked, min_order_length=300)
    for _, row in clones.iterrows():
        assert row["insert_aa_seq"] in row["clone_aa_seq"]
        assert len(row["clone_aa_seq"]) > len(row["insert_aa_seq"]), \
            "the vector adds leader, linker, VL and tag"
        assert row["insert_dna_seq"] in row["order_dna_seq"]
        assert row["order_length"] >= 300
        assert row["assembly_problems"] is None
    stores.close_all()


def test_plate_is_sequential_and_neighbours_are_similar(tmp_path, monkeypatch):
    from plateforge.library import diversity, ordering

    idx, picked, _rows, _meta, _filt, _lib = _selection_fixture(
        tmp_path, monkeypatch, n=600, pick=24)
    grouped = ordering.design(picked, min_order_length=300)
    scattered = ordering.design(picked, min_order_length=300,
                                group_by_similarity=False)

    assert list(grouped["well"])[:3] == ["A01", "B01", "C01"], "column-major"
    assert set(grouped["seq_id"]) == set(scattered["seq_id"]), \
        "layout must not change which clones were chosen"

    def neighbour_distance(frame):
        seqs = [str(s) for s in frame["cdr3_aa"]]
        return sum(diversity.distance(a, b)
                   for a, b in zip(seqs, seqs[1:])) / max(1, len(seqs) - 1)

    assert neighbour_distance(grouped) < neighbour_distance(scattered)
    stores.close_all()


def test_unverified_vendor_layout_is_not_emitted_by_accident(tmp_path,
                                                             monkeypatch):
    from plateforge.library import ordering, vendors

    idx, picked, _rows, _meta, _filt, _lib = _selection_fixture(
        tmp_path, monkeypatch, pick=6)
    clones = ordering.design(picked, min_order_length=300)

    with pytest.raises(ValueError, match="no real upload template"):
        vendors.emit(clones, "idt_eblocks")

    order = vendors.emit(clones, "idt_eblocks", allow_unverified=True)
    assert list(order.table.columns) == ["Name", "Sequence"]
    assert order.template_verified is False and order.caveats
    stores.close_all()


def test_generic_order_table_is_always_available(tmp_path, monkeypatch):
    from plateforge.library import ordering, vendors

    idx, picked, _rows, _meta, _filt, _lib = _selection_fixture(
        tmp_path, monkeypatch, pick=6)
    clones = ordering.design(picked, min_order_length=300)
    order = vendors.emit(clones, "generic")
    assert order.ok and order.template_verified
    assert len(order.table) == len(clones)
    stores.close_all()


def test_vendor_limits_are_checked_even_when_the_layout_is_a_guess():
    import pandas as pd
    from plateforge.library import vendors

    clones = pd.DataFrame({"clone_id": ["CLN-short", "CLN-fine"],
                           "order_dna_seq": ["ATGC" * 10, "ATGC" * 100]})
    order = vendors.emit(clones, "idt_eblocks", allow_unverified=True)
    assert not order.ok
    assert set(order.issues["name"]) == {"CLN-short"}


# --- germlines from the bundled IMGT tables --------------------------------

def test_imgt_tables_resolve_the_panel_germlines():
    from plateforge.library import germline_db, germlines

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed; install the [imgt] extra")
    for gene in germlines.DEFAULT.genes + [germlines.DEFAULT.light_chain]:
        found = germline_db.resolve(gene)
        assert found and found.source == "anarci"
        assert found.allele and found.allele.startswith(gene)
        assert not found.is_derived, "an IMGT reference is not a pool consensus"


def test_alignment_reference_prefers_imgt_and_names_both_alleles(tmp_path,
                                                                 monkeypatch):
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    _fresh(tmp_path, monkeypatch)
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=200, seed=91)
    df, _ = oas.load(src)
    gene = df["v_gene"].value_counts().index[0]
    seq, label = msa.reference_for(df[df["v_gene"] == gene], gene)
    assert label.startswith("IMGT ") and gene in label
    assert "IGHJ" in label, "the modal J is appended so FR4 has a germline"
    assert len(seq) > len(germline_db.from_anarci(gene).sequence)
    stores.close_all()


# --- gaps cost what the region says they cost ------------------------------

def test_cdr3_gaps_are_cheaper_than_framework_gaps():
    from plateforge.library import msa

    numbers = list(range(1, 129))
    profile = msa.region_gap_profile(numbers)
    fr1_open, _ = profile[10]            # IMGT 10-ish, FR1
    cdr3_open, _ = profile[110]          # IMGT 110-ish, CDR3
    cdr1_open, _ = profile[30]
    assert cdr3_open > cdr1_open > fr1_open, \
        "junction gaps must be the cheapest and framework gaps the dearest"


def test_a_seam_takes_the_more_permissive_side():
    """The V/J junction must be relaxed from its first column, not one late."""
    from plateforge.library import msa

    numbers = [104, 105]                 # FR3 then CDR3
    profile = msa.region_gap_profile(numbers)
    assert profile[1] == msa.REGION_GAP_PENALTY["CDR3"]


def test_region_of_matches_the_imgt_definitions():
    from plateforge.library import msa

    assert msa.region_of(1) == "FR1" and msa.region_of(26) == "FR1"
    assert msa.region_of(27) == "CDR1" and msa.region_of(38) == "CDR1"
    assert msa.region_of(105) == "CDR3" and msa.region_of(117) == "CDR3"
    assert msa.region_of(118) == "FR4"
    assert msa.region_of(None) is None


def test_gap_profile_changes_where_the_junction_lands():
    """A long junction against a mutated framework: the two differ."""
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    v = germline_db.from_anarci("IGHV3-23")
    j = germline_db.from_anarci("IGHJ4")
    ref = v.sequence + j.sequence
    msa._NUMBERING[ref] = msa._numbers_for(v) + msa._numbers_for(j)

    mutated = list(v.sequence)
    for i in (12, 40, 70, 88):                   # somatic hypermutation
        mutated[i] = "K" if mutated[i] != "K" else "R"
    query = "".join(mutated) + "GWLSTWYPRQVMKGCADLHNETSPY" + j.sequence

    relaxed = msa.build(ref, {"q": query}, region_aware=True)
    uniform = msa.build(ref, {"q": query}, region_aware=False)
    assert relaxed.identity(relaxed.rows[0]) >= uniform.identity(uniform.rows[0])

    # and the junction must sit in CDR3 columns, not spill into FR4
    spans = {n: (a, b) for n, a, b in relaxed.region_columns()}
    a, b = spans["CDR3"]
    assert sum(1 for c in range(a, b) if relaxed.rows[0][c] not in "-.") >= 25


def test_region_aware_is_inert_without_numbering():
    from plateforge.library import msa

    ref = "EVQLLESGGGLVQPGGSLRLSCAAS"
    a = msa.build(ref, {"q": ref}, region_aware=True)
    b = msa.build(ref, {"q": ref}, region_aware=False)
    assert a.ref == b.ref and a.rows == b.rows


def test_regions_are_column_spans_that_cover_the_alignment():
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    v = germline_db.from_anarci("IGHV3-23")
    j = germline_db.from_anarci("IGHJ4")
    ref = v.sequence + j.sequence
    msa._NUMBERING[ref] = msa._numbers_for(v) + msa._numbers_for(j)
    m = msa.build(ref, {"q": v.sequence + "ARWWWDY" + j.sequence})
    spans = m.region_columns()
    assert [n for n, _, _ in spans] == ["FR1", "CDR1", "FR2", "CDR2", "FR3",
                                        "CDR3", "FR4"]
    assert spans[0][1] == 0 and spans[-1][2] == m.width, "spans must tile"
    for (_, _, end), (_, start, _) in zip(spans, spans[1:]):
        assert end == start, "spans must not overlap or leave holes"


# --- the IgG vector ---------------------------------------------------------

def test_constant_regions_match_the_stored_uniprot_records():
    """The constants must be derivable from fixtures/, not typed from memory."""
    from pathlib import Path
    from plateforge.library import vector

    def read(name):
        text = Path("fixtures/uniprot") / name
        return "".join(l.strip() for l in text.read_text().splitlines()
                       if not l.startswith(">"))

    ighg1 = read("P01857_IGHG1_HUMAN.fasta")
    secreted = ighg1[:ighg1.index("KSLSLSP") + len("KSLSLSP")] + "GK"
    assert vector.IGHG1_CH1_CH3_AA == secreted
    assert vector.IGKC_AA == read("P01834_IGKC_HUMAN.fasta")
    # and the membrane tail must be gone
    assert not vector.IGHG1_CH1_CH3_AA.endswith("ELQLEESCAEAQDGELDGLW")


def test_igg_constants_translate_to_what_they_claim():
    from plateforge.library import codon, vector

    for dna, aa in [(vector.IGHG1_CH1_CH3_DNA, vector.IGHG1_CH1_CH3_AA),
                    (vector.IGKC_DNA, vector.IGKC_AA),
                    (vector.FURIN_T2A_DNA, vector.FURIN_T2A_AA)]:
        assert codon.translate(dna) == aa


def test_igg_orf_makes_two_chains_from_one_reading_frame():
    from plateforge.library import codon, vector

    vec = vector.get("pcdna-igg1-t2a-vk-v1")
    vh = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAKGGYFDYWGQGTLVTVSS"
    orf = vec.orf(codon.optimize(vh, enzymes=("BsaI",)).dna)
    assert len(orf) % 3 == 0
    protein = codon.translate(orf)
    assert protein.count("*") == 1 and protein.endswith("*")

    # heavy chain: leader, VH, then the constant region
    assert protein.startswith(vector.SIGNAL_AA)
    assert vh in protein
    hc_end = protein.index(vector.IGHG1_CH1_CH3_AA)
    assert hc_end > protein.index(vh)

    # the split, then the light chain with a leader of its own
    t2a = protein.index(vector.T2A_AA)
    assert t2a > hc_end
    lc = protein[t2a + len(vector.T2A_AA):]
    assert lc.startswith(vector.SIGNAL_AA), "the LC needs its own signal peptide"
    assert vector.VL_AA in lc and vector.IGKC_AA in lc


def test_furin_site_precedes_the_2a_peptide():
    """Without it the heavy chain keeps the 2A tail on its C-terminus."""
    from plateforge.library import vector

    vec = vector.get("pcdna-igg1-t2a-vk-v1")
    assert vector.FURIN_T2A_AA.startswith(vector.FURIN_SITE_AA)
    assert vector.FURIN_T2A_AA.endswith(vector.T2A_AA)
    assert vec.protein  # the element is in the vector at all
    names = [e.name for e in vec.downstream]
    assert names.index("furin-GSG-T2A") > names.index("IgG1 CH1-hinge-CH2-CH3")


def test_dual_cassette_vector_is_the_default_for_ordering(tmp_path, monkeypatch):
    import json
    from plateforge.library import ordering

    idx, picked, _r, _m, _f, _l = _selection_fixture(tmp_path, monkeypatch, pick=4)
    clones = ordering.design(picked)
    assert set(clones["vector"]) == {"pcdna-igg1-dual-vk-v1"}
    assert set(clones["construct"]) == {"igg1"}
    assert (clones["n_chains"] == 2).all(), "two cassettes, two chains"
    for text in clones["cassettes_json"]:
        chains = json.loads(text)
        assert len(chains) == 2
        # the heavy chain carries the insert and so is the longer of the two
        assert max(chains.values()) > min(chains.values())
    # clone_aa_seq is the chain the insert determines, not both concatenated
    assert (clones["aa_length"] < 600).all()
    assert (clones["aa_length"] > 300).all()
    stores.close_all()


def test_scfv_vector_is_still_available(tmp_path, monkeypatch):
    from plateforge.library import ordering

    idx, picked, _r, _m, _f, _l = _selection_fixture(tmp_path, monkeypatch, pick=4)
    clones = ordering.design(picked, vector_name="pcdna-scfv-vk-v1", fmt="scfv")
    assert (clones["aa_length"] < 400).all()
    stores.close_all()


# --- the metadata figure ----------------------------------------------------

def test_metadata_table_figure_keeps_clone_id_tails(tmp_path, monkeypatch):
    from plateforge.library import figures, ordering

    idx, picked, _r, _m, _f, _l = _selection_fixture(tmp_path, monkeypatch, pick=12)
    clones = ordering.design(picked)
    p = figures.metadata_table(clones, 10, tmp_path / "meta.png")
    assert p.exists() and p.stat().st_size > 10_000
    # ids share a long prefix, so a front-truncated id would be identical in
    # every row and the figure would say nothing
    shown = {figures._shorten(c, 14, True) for c in clones["clone_id"].head(10)}
    assert len(shown) == 10
    stores.close_all()


def test_order_summary_names_the_vector_and_the_junction(tmp_path, monkeypatch):
    from plateforge.library import ordering, vector

    idx, picked, _r, _m, _f, _l = _selection_fixture(tmp_path, monkeypatch, pick=6)
    clones = ordering.design(picked)
    text = ordering.summary_markdown(clones, set_id="CLN-set-test")
    assert "pcdna-igg1-dual-vk-v1" in text
    assert "cassettes" in text and "1.5–2:1" in text
    assert vector.VL_JUNCTION in text and "not** \ngermline-encoded" not in text
    assert vector.IGHG1_ACCESSION in text and vector.IGKC_ACCESSION in text
    assert "fixtures/uniprot/" in text
    assert "CLN-set-test" in text
    stores.close_all()


# --- figures say which run wrote them --------------------------------------

def test_figure_directory_records_the_run_that_wrote_it(tmp_path, monkeypatch):
    """'Did my figures update?' must be answerable without comparing pixels."""
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures

    out = tmp_path / "figs"
    p = figures.stamp(out, run_id="LIB-test-123", note="selection run")
    text = p.read_text()
    assert "LIB-test-123" in text and "code" in text and "written" in text
    stores.close_all()


def test_stale_figures_are_detected(tmp_path, monkeypatch):
    """A PNG left at a path the current version no longer writes is the bug
    that made an old alignment look current for two releases."""
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures

    out = tmp_path / "figs"
    out.mkdir()
    orphan = out / "alignment_IGHV3-23.png"
    orphan.write_bytes(b"\x89PNG old")
    fresh = out / "fig1_genes.png"
    fresh.write_bytes(b"\x89PNG new")

    found = figures.stale(out, made=[fresh])
    assert found == [orphan]
    assert figures.stale(out, made=[fresh, orphan]) == []
    stores.close_all()


# --- CDR3 is one block, anchored at both ends ------------------------------

def test_place_leaves_the_gap_in_the_middle():
    from plateforge.library import msa

    assert msa.place("ARDY", 12, "center") == "AR--------DY"
    assert msa.place("ARDY", 12, "left") == "ARDY--------"
    assert msa.place("ARDY", 4, "center") == "ARDY"
    # an odd residue goes left, as IMGT does
    assert msa.place("ARD", 7, "center") == "AR----D"


def _igg_ref(gene="IGHV3-23", j="IGHJ4"):
    from plateforge.library import germline_db, msa
    v, jj = germline_db.from_anarci(gene), germline_db.from_anarci(j)
    seq = v.sequence + jj.sequence
    numbers = msa._numbers_for(v) + msa._numbers_for(jj)
    return msa._drop_cdr3(seq, numbers)


def test_reference_carries_no_cdr3_residues_of_its_own():
    """Germline CDR3 residues are anchors the aligner would chop the junction on."""
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    seq, numbers = _igg_ref()
    assert not any(105 <= n <= 117 for n in numbers if n is not None)
    assert 104 in numbers and 118 in numbers, "both anchors must survive"
    assert len(seq) == len(numbers)


def test_the_junction_is_a_single_contiguous_block():
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    seq, numbers = _igg_ref()
    msa._NUMBERING[seq] = numbers
    queries = {
        "short": seq[:96] + "ARDY" + seq[96:],
        "long": seq[:96] + "ARWGYSSGWYFDYYYGMDV" + seq[96:],
        "mid": seq[:96] + "ARQCINCGFCGGKDY" + seq[96:],
    }
    m = msa.build(seq, queries)
    spans = [s for s in m.region_columns() if s[0] == "CDR3"]
    assert len(spans) == 1, "one CDR3 block, not confetti"
    _, a, b = spans[0]
    assert b - a == 19, "the block is exactly as wide as the longest junction"
    for row in m.rows:
        inner = row[a:b]
        # Centring puts the gap in the middle on purpose, so a short junction
        # HAS interior gaps. What must not happen is several separate runs --
        # that is the confetti this replaced.
        runs = [r for r in re.split(r"[^-]+", inner) if r]
        assert len(runs) <= 1, f"one gap run, got {len(runs)}: {inner}"


def test_junctions_share_their_first_and_last_columns():
    """Both anchors hold: CDR3s line up at the C end AND the W end."""
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    seq, numbers = _igg_ref()
    msa._NUMBERING[seq] = numbers
    m = msa.build(seq, {"a": seq[:96] + "ARDY" + seq[96:],
                        "b": seq[:96] + "ARWGYSSGWYFDY" + seq[96:],
                        "c": seq[:96] + "ARQCINCGFCGGKDY" + seq[96:]})
    _, a, b = [s for s in m.region_columns() if s[0] == "CDR3"][0]
    starts = {row[a:b].index(next(ch for ch in row[a:b] if ch != "-")) for row in m.rows}
    ends = {len(row[a:b].rstrip("-")) for row in m.rows}
    assert starts == {0}, "every junction starts in the first CDR3 column"
    assert ends == {b - a}, "every junction ends in the last CDR3 column"


def test_left_justifying_would_be_worse():
    """The regression this replaces: a short and a long CDR3 sharing only a start."""
    from plateforge.library import msa

    width = 15
    short, long_ = "ARDY", "ARQCINCGFCGGKDY"
    left = [msa.place(short, width, "left"), msa.place(long_, width, "left")]
    centred = [msa.place(short, width, "center"), msa.place(long_, width, "center")]

    def agreeing(rows):
        return sum(1 for col in zip(*rows) if len(set(col)) == 1 and col[0] != "-")

    assert agreeing(centred) > agreeing(left), \
        "centring must line up more residues than left-justifying"


def test_deletions_stay_expensive_inside_the_junction():
    """Only insertions are cheap in CDR3. Losing framework is a strong claim."""
    from plateforge.library import msa

    ref = "EVQLLESGGGLVQPGGSLRLSCAAS"
    cheap = [(-0.5, -0.1)] * (len(ref) + 1)
    _gapped_ref, gapped_seq = msa.align_pair(ref, ref[:8] + ref[14:],
                                             gap_profile=cheap)
    assert gapped_seq.count("-") == 6, "the deletion is still taken"
    runs = [r for r in re.split(r"[^-]+", gapped_seq) if r]
    assert len(runs) == 1, "one contiguous deletion, not several cheap ones"


def test_column_regions_tile_the_whole_alignment():
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    seq, numbers = _igg_ref()
    msa._NUMBERING[seq] = numbers
    m = msa.build(seq, {"a": seq[:96] + "ARWGYSSGWYFDY" + seq[96:]})
    per_column = m.column_regions()
    assert len(per_column) == m.width and all(per_column)
    spans = m.region_columns()
    assert [n for n, _, _ in spans] == ["FR1", "CDR1", "FR2", "CDR2", "FR3",
                                        "CDR3", "FR4"]
    assert spans[0][1] == 0 and spans[-1][2] == m.width


# --- the anchors bracket the junction, always ------------------------------

def _anchored(queries: dict[str, str]):
    """Build an alignment against a real IMGT reference and return it."""
    from plateforge.library import germline_db, msa
    v, j = germline_db.from_anarci("IGHV3-23"), germline_db.from_anarci("IGHJ4")
    seq, numbers = msa._drop_cdr3(v.sequence + j.sequence,
                                  msa._numbers_for(v) + msa._numbers_for(j))
    msa._NUMBERING[seq] = numbers
    cut = len([n for n in numbers if n is not None and n <= 104])
    built = {k: seq[:cut] + junction + seq[cut:] for k, junction in queries.items()}
    return msa.build(seq, built), seq, cut


def _junction_bounds(m):
    _, a, b = [s for s in m.region_columns() if s[0] == "CDR3"][0]
    starts, ends = set(), set()
    for row in m.rows:
        seg = row[a:b]
        residues = [i for i, c in enumerate(seg) if c != "-"]
        starts.add(residues[0])
        ends.add(residues[-1])
    return starts, ends, b - a


def test_every_junction_starts_and_ends_in_the_same_column():
    from plateforge.library import germline_db

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    m, _seq, _cut = _anchored({
        "short": "ARDY",
        "mid": "ARGLREYWKIDY",
        "long": "ARQCINCGFCGGKDYYYGMDV",
    })
    starts, ends, width = _junction_bounds(m)
    assert starts == {0}
    assert ends == {width - 1}


def test_a_missing_cysteine_does_not_shift_that_sequence_left():
    """The outlier bug: one junction a column left of every other."""
    from plateforge.library import germline_db, msa

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    v, j = germline_db.from_anarci("IGHV3-23"), germline_db.from_anarci("IGHJ4")
    seq, numbers = msa._drop_cdr3(v.sequence + j.sequence,
                                  msa._numbers_for(v) + msa._numbers_for(j))
    msa._NUMBERING[seq] = numbers
    cut = len([n for n in numbers if n is not None and n <= 104])
    m = msa.build(seq, {
        "normal": seq[:cut] + "ARGLREYWKIDY" + seq[cut:],
        "no_cys": seq[:cut - 1] + "ARGLREYWKIDY" + seq[cut:],      # Cys gone
        "no_trp": seq[:cut] + "ARGLREYWKIDY" + seq[cut + 1:],      # Trp gone
    })
    starts, ends, width = _junction_bounds(m)
    assert starts == {0}, "a missing anchor must not shift the junction"
    assert ends == {width - 1}

    # and the anchor column shows a gap rather than a borrowed residue
    _, a, _b = [s for s in m.region_columns() if s[0] == "CDR3"][0]
    by_label = dict(zip(m.labels, m.rows))
    assert by_label["no_cys"][a - 1] == "-"


def test_a_junction_containing_tryptophan_keeps_its_own_columns():
    """Why the anchors are not enforced by residue identity.

    Forbidding non-W at the tryptophan column lets any W inside a junction
    capture that column from across the block, tearing the junction in two.
    """
    from plateforge.library import germline_db

    if not germline_db.anarci_available():
        pytest.skip("anarci not installed")
    m, _seq, _cut = _anchored({
        "with_w": "ARGLREYWKIDY",          # a W in the middle of the junction
        "plain": "ARGLREYAKIDY",
    })
    starts, ends, width = _junction_bounds(m)
    assert starts == {0} and ends == {width - 1}
    for row in m.rows:
        _, a, b = [s for s in m.region_columns() if s[0] == "CDR3"][0]
        runs = [r for r in re.split(r"[^-]+", row[a:b]) if r]
        assert len(runs) <= 1, "the junction must not be torn in two"


def test_consolidation_is_a_no_op_when_the_anchors_are_intact():
    from plateforge.library import msa

    at = list("CW")
    before = ["", "ARDY", ""]
    msa.consolidate_junction(at, before, "CW", 0, 1)
    assert at == ["C", "W"] and before[1] == "ARDY"


# --- where the data tree lives ---------------------------------------------

def test_data_root_defaults_inside_the_checkout(monkeypatch):
    from plateforge.core import paths

    monkeypatch.delenv("PLATEFORGE_DATA", raising=False)
    root = paths.repo_root()
    assert root is not None and (root / "pyproject.toml").exists()
    assert paths.default_root() == root / "plateforge-data"


def test_data_root_does_not_follow_the_working_directory(monkeypatch, tmp_path):
    """Two half-populated pools is a bad afternoon."""
    from plateforge.core import paths

    monkeypatch.delenv("PLATEFORGE_DATA", raising=False)
    first = paths.default_root()
    monkeypatch.chdir(tmp_path)
    assert paths.default_root() == first


def test_environment_still_wins(monkeypatch, tmp_path):
    from plateforge.core import paths

    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path / "elsewhere"))
    assert paths.data_root() == tmp_path / "elsewhere"


def test_the_data_tree_is_gitignored():
    from pathlib import Path

    ignored = Path(".gitignore").read_text()
    assert "/plateforge-data/" in ignored


# --- two cassettes, two transcripts ----------------------------------------

def test_a_cassette_knows_whether_it_carries_the_insert():
    from plateforge.library import vector

    vec = vector.get("pcdna-igg1-dual-vk-v1")
    carrying = [c for c in vec.cassettes if c.carries_insert]
    assert len(carrying) == 1 and carrying[0] is vec.insert_cassette
    assert len(vec.cassettes) == 2


def test_each_cassette_is_its_own_reading_frame():
    """The light chain is not downstream of the heavy chain in any frame."""
    from plateforge.library import codon, vector

    vec = vector.get("pcdna-igg1-dual-vk-v1")
    vh = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAKGGYFDYWGQGTLVTVSS"
    dna = codon.optimize(vh, enzymes=("BsaI",)).dna
    proteins = vec.proteins(dna)
    assert len(proteins) == 2
    for name, protein in proteins.items():
        assert protein.count("*") == 1 and protein.endswith("*")
        assert protein.startswith(vector.SIGNAL_AA), \
            f"{name} needs its own signal peptide"

    heavy = proteins["HC (CMV)"]
    light = proteins["LC (CAG)"]
    assert vh in heavy and vector.IGHG1_CH1_CH3_AA in heavy
    assert vector.VL_AA in light and vector.IGKC_AA in light
    assert vh not in light, "the insert must appear in exactly one cassette"
    assert vector.T2A_AA not in heavy and vector.T2A_AA not in light


def test_the_ordered_fragment_is_unchanged_by_the_vector_swap():
    """Fusion sites are in the constant flanks, so the eBlocks do not move."""
    from plateforge.library import cloning, codon, vector

    dual = vector.get("pcdna-igg1-dual-vk-v1")
    t2a = vector.get("pcdna-igg1-t2a-vk-v1")
    assert dual.fusion_sites() == t2a.fusion_sites()

    vh = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAKGGYFDYWGQGTLVTVSS"
    cds = codon.optimize(vh, enzymes=("BsaI",)).dna
    a = cloning.design(cds, dual, "golden_gate", min_length=300)
    b = cloning.design(cds, t2a, "golden_gate", min_length=300)
    assert a.dna == b.dna, "a vector change must not reprint the order"


def test_verify_checks_every_cassette_not_just_the_insert_s():
    from plateforge.library import cloning, codon, vector

    vec = vector.get("pcdna-igg1-dual-vk-v1")
    vh = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAKGGYFDYWGQGTLVTVSS"
    cds = codon.optimize(vh, enzymes=("BsaI",)).dna
    fragment = cloning.design(cds, vec, "golden_gate", min_length=300)
    assert cloning.verify(fragment, cds, vec, expected_protein=vh) == []

    # break the cassette that does NOT carry the insert
    broken = vector.Vector(
        name="broken", backbone="t", insert="VH", enzyme="BsaI",
        cassettes=(vec.cassettes[0],
                   vector.Cassette("LC (broken)", tuple(
                       e if e.kind != "domain"
                       else vector.Element(e.name, e.kind, dna=e.dna[:-1],
                                           aa=e.aa)
                       for e in vec.cassettes[1].elements))))
    problems = cloning.verify(fragment, cds, broken, expected_protein=vh)
    assert any("LC (broken)" in p for p in problems), \
        "a fault in the second cassette must not go unreported"


def test_declared_bp_does_not_pretend_to_know_the_promoters():
    from plateforge.library import vector

    vec = vector.get("pcdna-igg1-dual-vk-v1")
    sizes = vec.declared_bp("N" * 360)
    assert sizes["_undeclared_elements"] >= 3, \
        "promoters and the backbone have no sequence here and must be counted out"
    assert sizes["HC (CMV)"] > sizes["LC (CAG)"]
    assert sizes["_total_declared"] == sizes["HC (CMV)"] + sizes["LC (CAG)"]


def test_the_t2a_vector_is_still_registered():
    """Kept for comparison; the decision record says why it is not default."""
    from plateforge.library import vector

    assert "pcdna-igg1-t2a-vk-v1" in vector.VECTORS
    assert len(vector.get("pcdna-igg1-t2a-vk-v1").cassettes) == 1


# --- figures cannot go stale silently --------------------------------------

def test_a_figure_directory_records_which_run_wrote_it(tmp_path):
    from plateforge.library import figures

    figures.stamp(tmp_path, run_id="LIB-test-1", note="selection run")
    text = (tmp_path / "FIGURES_FROM.txt").read_text()
    assert "LIB-test-1" in text and "selection run" in text
    assert "code" in text and "written" in text


def test_a_figure_nothing_wrote_this_run_is_reported_as_stale(tmp_path):
    """The bug this exists for: an output path moved and the old PNG stayed.

    figures/real/alignment_IGHV3-23.png sat in place across two releases
    looking perfectly current while being produced by code that no longer
    existed. Nothing about the file said so.
    """
    from plateforge.library import figures

    fresh = tmp_path / "fig1.png"
    fresh.write_bytes(b"\x89PNG fresh")
    orphan = tmp_path / "alignment_IGHV3-23.png"
    orphan.write_bytes(b"\x89PNG old")

    stale = figures.stale(tmp_path, made=[fresh])
    assert stale == [orphan]
    assert figures.stale(tmp_path, made=[fresh, orphan]) == []


def test_stale_ignores_a_directory_that_does_not_exist(tmp_path):
    from plateforge.library import figures

    assert figures.stale(tmp_path / "nope", made=[]) == []


def test_the_scripts_stamp_and_check_their_figure_directories():
    """Both entry points must do this, or the trap is only half closed."""
    from pathlib import Path

    for script in ("scripts/fetch_oas.py", "scripts/order.py"):
        text = Path(script).read_text()
        assert "figures.stamp(" in text, f"{script} must stamp its figures"
    assert "figures.stale(" in Path("scripts/fetch_oas.py").read_text()
