import pytest

from plateforge.core import artifacts, ids, registry, stores, wells


def test_well_normalization():
    assert wells.normalize("a1") == "A01"
    assert wells.normalize("A1") == "A01"
    assert wells.normalize(" H 12 ") == "H12"
    assert wells.normalize("1A") == "A01"
    assert wells.normalize("A01") == "A01"


def test_well_bounds():
    with pytest.raises(ValueError):
        wells.normalize("A13", 96)
    assert wells.normalize("A13", 384) == "A13"
    assert wells.normalize("P24", 384) == "P24"


def test_well_ordering():
    col = wells.all_wells(96, "column")
    assert col[:3] == ["A01", "B01", "C01"]
    assert len(col) == 96
    row = wells.all_wells(96, "row")
    assert row[:3] == ["A01", "A02", "A03"]


def test_edges():
    assert wells.is_edge("A01")
    assert wells.is_edge("H12")
    assert not wells.is_edge("D06")


def test_id_roundtrip():
    plate = ids.mint("PLT", "20260912-001")
    assert ids.parse(plate) == ("PLT", "20260912-001")
    assert ids.is_a(plate, "PLT")
    assert not ids.is_a(plate, "CLN")


def test_unregistered_prefix_fails_loudly():
    with pytest.raises(KeyError):
        ids.mint("ZZZ")
    with pytest.raises(KeyError):
        ids.parse("ZZZ-123")


def test_registry():
    reg = registry.Registry("widget")

    @reg.register("alpha", ext=".txt")
    def alpha():
        return 1

    assert reg.get("alpha")() == 1
    assert reg.meta("alpha")["ext"] == ".txt"
    assert "alpha" in reg
    with pytest.raises(KeyError):
        reg.get("beta")


def test_artifact_lineage(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()

    lib = artifacts.register(ids.mint("LIB", "demo"), "sequence_set", "library.build",
                             meta={"n": 10})
    plate = artifacts.register(ids.mint("PLT", "p1"), "plate", "assay.design",
                               parents={lib: "sampled_from"})
    xfr = artifacts.register(ids.mint("XFR", "x1"), "transfer_list", "emit.worklist",
                             parents={plate: "rearray_of"})

    assert artifacts.get(lib).meta["n"] == 10
    assert artifacts.parents(plate) == [(lib, "sampled_from")]
    assert lib in artifacts.lineage(xfr, "up")
    assert xfr in artifacts.lineage(lib, "down")
    assert len(artifacts.find(obj_type="plate")) == 1
    stores.close_all()


def test_meta_escape_hatch(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()
    oid = artifacts.register(ids.mint("RUN", "r1"), "analysis", "assay.analyze")
    artifacts.set_meta(oid, "unexpected_field", "recorded anyway")
    assert artifacts.get(oid).meta["unexpected_field"] == "recorded anyway"
    assert artifacts.find(unexpected_field="recorded anyway")
    stores.close_all()


def test_paths_layout(tmp_path, monkeypatch):
    from plateforge.core import paths
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    assert paths.store_path("assay").parent.name == "stores"
    assert paths.bulk_path("sequences").suffix == ".parquet"
    assert paths.raw_dir("oas").is_dir()
    assert paths.output_dir("PLT-x1").is_dir()


def test_bulk_roundtrip(tmp_path, monkeypatch):
    import pandas as pd
    from plateforge.core import bulk
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    df = pd.DataFrame({"seq_id": ["SEQ-a", "SEQ-b"], "chain": ["H", "H"], "aa": ["QVQ", "EVQ"]})
    bulk.write("sequences", df)
    assert bulk.exists("sequences")
    back = bulk.read("sequences", columns=["seq_id", "chain"])
    assert list(back["seq_id"]) == ["SEQ-a", "SEQ-b"]
    assert bulk.info("sequences")["rows"] == 2


def test_code_version_recorded(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()
    oid = artifacts.register(ids.mint("RUN", "cv"), "analysis", "assay.analyze")
    assert artifacts.get(oid).code_version
    stores.close_all()


def test_supersession(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()
    first = artifacts.register(ids.mint("RUN", "v1"), "analysis", "assay.analyze")
    second = artifacts.register(ids.mint("RUN", "v2"), "analysis", "assay.analyze")
    artifacts.supersede(first, second)
    assert artifacts.is_superseded(first)
    assert not artifacts.is_superseded(second)
    assert [a.obj_id for a in artifacts.current("analysis")] == [second]
    stores.close_all()


def test_stamped_ids_do_not_collide():
    made = {ids.mint_stamped("RUN", "ingest") for _ in range(200)}
    assert len(made) == 200
    assert all(ids.is_a(m, "RUN") for m in made)
