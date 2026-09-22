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
