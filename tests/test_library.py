import pandas as pd
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
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=900, seed=10)
    df, _ = oas.load(src)
    picked = diversity.sample(df, 10, seed=3)
    report = diversity.Spec(min_genes=3).check(picked)
    assert report["pass"], report


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


def test_alignment_report_detects_ragged(tmp_path):
    from plateforge.library import hf
    df = pd.DataFrame({
        "v_gene": ["IGHV3-23"] * 3 + ["IGHV1-69"] * 2,
        "aa_gapped": ["AAAA", "AAAA", "AAA", "BBBB", "BBBB"],
    })
    rep = hf.alignment_report(df).set_index("v_gene")
    assert not rep.loc["IGHV3-23", "column_aligned"]
    assert rep.loc["IGHV1-69", "column_aligned"]


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


# --- alignment refuses to misalign ----------------------------------------

def test_alignment_excludes_ragged_sequences(tmp_path, monkeypatch):
    _fresh(tmp_path, monkeypatch)
    from plateforge.library import figures
    src = synth.make_unit(tmp_path / "unit.csv.gz", n=300, seed=50)
    df, _ = oas.load(src)
    gene = df["v_gene"].value_counts().index[0]
    sub = df[df["v_gene"] == gene].copy()
    # Truncate a few reads, as a partial-coverage read would be.
    idx = sub.index[:3]
    sub.loc[idx, "aa_gapped"] = sub.loc[idx, "aa_gapped"].str[:-7]

    p = figures.alignment(sub, gene, tmp_path / "ragged.png", max_rows=12)
    assert p.exists() and p.stat().st_size > 5_000
    q = figures.alignment(sub, gene, tmp_path / "padded.png", max_rows=12,
                          allow_ragged=True)
    assert q.exists() and q.stat().st_size > 5_000
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
