"""The plate lineage layer: containers, steps, and what happened.

The property nearly every test here is really checking is the one from
decision 0019 -- the plate is not the sample. A clone is tracked; plates are
things it passes through.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from plateforge.assay import plates, protocol, steps
from plateforge.core import artifacts, stores


@pytest.fixture
def fresh(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    stores.close_all()
    yield tmp_path
    stores.close_all()


def _dna_plate(n: int = 96, barcode: str = "PF-DNA-001") -> plates.Plate:
    plate = plates.create("dna", 96, barcode=barcode, label="order")
    wells = [f"{r}{c:02d}" for c in range(1, 13) for r in "ABCDEFGH"][:n]
    rows = pd.DataFrame({"well": wells,
                         "clone_id": [f"CLN-{i:03d}" for i in range(n)]})
    plates.fill(plate, rows, volume_ul=50.0)
    plates.save(plate)
    return plate


def _run_chain(start: plates.Plate, through: list[str],
               concentrations: dict | None = None,
               barcodes: dict | None = None) -> plates.Plate:
    current = start
    for name in through:
        plan = steps.get(name).plan(current, barcode=(barcodes or {}).get(name))
        steps.commit(plan)
        steps.record(plan, operator="test")
        if plan.outputs:
            current = plan.outputs[0]
        if name == "quantify_bli" and concentrations is not None:
            frame = current.wells
            for i in frame.index:
                if frame.loc[i, "role"] != "empty":
                    frame.loc[i, "concentration"] = concentrations.get(
                        frame.loc[i, "well"], 20.0)
                    frame.loc[i, "conc_units"] = "ug/mL"
            plates.save(current)
    return current


# --- plates -----------------------------------------------------------------

def test_a_new_plate_has_every_well_and_they_are_all_empty(fresh):
    plate = plates.create("dna", 96)
    assert len(plate.wells) == 96
    assert set(plate.wells["role"]) == {"empty"}
    assert list(plate.wells["well"])[:3] == ["A01", "B01", "C01"]
    assert plate.occupied.empty


def test_wells_are_normalised_on_the_way_in(fresh):
    """Rule 3: A01, never A1. The caller may not know that."""
    plate = plates.create("dna", 96)
    plates.fill(plate, pd.DataFrame({"well": ["a1", "B2"],
                                     "clone_id": ["CLN-1", "CLN-2"]}))
    got = plate.wells[plate.wells["role"] == "sample"]
    assert set(got["well"]) == {"A01", "B02"}


def test_a_plate_round_trips_through_the_store(fresh):
    plate = _dna_plate(12)
    back = plates.load(plate.plate_id)
    assert back.barcode == plate.barcode
    assert back.content_kind == "dna"
    assert back.clones == plate.clones
    assert back.map_of("clone_id") == plate.map_of("clone_id")


def test_a_plate_is_a_registered_artifact(fresh):
    plate = _dna_plate(4)
    record = artifacts.get(plate.plate_id)
    assert record.obj_type == "plate"
    assert record.meta["barcode"] == "PF-DNA-001"


# --- the chain --------------------------------------------------------------

def test_the_whole_chain_runs_and_keeps_the_layout(fresh):
    """The property that matters: the IgG prep has the ordering plate's map."""
    dna = _dna_plate()
    final = _run_chain(dna, protocol.DEFAULT_CHAIN, concentrations={})
    assert final.content_kind == "dilution"
    assert plates.layout_matches(dna, final), \
        "every clone must still be in the well it was ordered into"


def test_every_step_declares_what_it_moves_between(fresh):
    table = steps.available().set_index("step")
    for name in protocol.DEFAULT_CHAIN:
        assert name in table.index
        assert table.loc[name, "produces"], f"{name} must declare an output kind"


def test_a_step_refuses_the_wrong_content(fresh):
    """A step that expects supernatant must say so before a robot moves anything."""
    dna = _dna_plate(4)
    with pytest.raises(ValueError, match="consumes 'supernatant'"):
        steps.get("add_beads").plan(dna)


def test_steps_chain_by_content_kind(fresh):
    """Each step's output is the next step's input. That is the contract."""
    table = steps.available().set_index("step")
    chain = protocol.DEFAULT_CHAIN
    for earlier, later in zip(chain, chain[1:]):
        assert table.loc[earlier, "produces"] == table.loc[later, "consumes"], \
            f"{earlier} -> {later} does not connect"


# --- lineage ----------------------------------------------------------------

def test_a_well_traces_back_to_the_plate_it_was_ordered_on(fresh):
    dna = _dna_plate()
    final = _run_chain(dna, protocol.DEFAULT_CHAIN, concentrations={})
    chain = plates.trace(final.plate_id, "D07")

    assert chain[0]["plate_id"] == final.plate_id
    assert chain[-1]["plate_id"] == dna.plate_id
    assert {hop["clone_id"] for hop in chain} == {dna.well("D07")["clone_id"]}, \
        "one clone the whole way; a broken link would show as a second id"
    kinds = [hop["content_kind"] for hop in chain]
    assert kinds[0] == "dilution" and kinds[-1] == "dna"


def test_lineage_is_per_well_so_a_rearray_is_not_a_special_case(fresh):
    """The reason source_plate_id lives on the well and not the plate."""
    source = _dna_plate(8)
    dest = plates.create("dna", 96, label="rearray")
    # deliberately scramble: A01 -> H12, B01 -> G12, ...
    moved = [("A01", "H12"), ("B01", "G12"), ("C01", "F12")]
    frame = dest.wells.copy()
    index = {w: i for i, w in enumerate(frame["well"])}
    for src, dst in moved:
        row = source.well(src)
        i = index[dst]
        frame.loc[i, ["clone_id", "role", "source_plate_id", "source_well"]] = [
            row["clone_id"], "sample", source.plate_id, src]
    dest.wells = frame
    plates.save(dest)

    for src, dst in moved:
        chain = plates.trace(dest.plate_id, dst)
        assert len(chain) == 2
        assert chain[1]["plate_id"] == source.plate_id
        assert chain[1]["well"] == src
        assert chain[0]["clone_id"] == chain[1]["clone_id"]


def test_where_is_finds_a_clone_in_every_container_it_passed_through(fresh):
    dna = _dna_plate()
    _run_chain(dna, protocol.DEFAULT_CHAIN, concentrations={})
    clone = dna.well("A01")["clone_id"]
    found = plates.where_is(clone)
    assert set(found["content_kind"]) >= {"dna", "supernatant", "igg_prep",
                                          "dilution"}
    assert set(found["well"]) == {"A01"}


def test_one_barcode_can_hold_several_states(fresh):
    """Adding beads to a plate does not make it a different piece of plastic.

    Each state is its own PLT because artifacts are never modified (rule 10).
    They share a barcode because the container is the same. Making the column
    UNIQUE meant INSERT OR REPLACE deleted the earlier state and broke the
    trace halfway along -- which is how this was found.
    """
    dna = _dna_plate()
    _run_chain(dna, ["transfect", "harvest_supernatant", "add_beads",
                     "immunoprecipitate", "elute_neutralise"],
               barcodes={"harvest_supernatant": "PF-SUP-001"})
    states = plates.states_of("PF-SUP-001")
    assert len(states) == 4, "supernatant, beads, washed beads, eluate"
    assert list(states["content_kind"]) == ["supernatant", "beads", "beads",
                                            "eluate"]
    assert states["plate_id"].nunique() == len(states), "distinct artifacts"


# --- steps as records -------------------------------------------------------

def test_plan_has_no_side_effects_until_commit(fresh):
    dna = _dna_plate(4)
    plan = steps.get("transfect").plan(dna)
    assert steps.history().empty
    with pytest.raises(KeyError):
        plates.load(plan.outputs[0].plate_id)
    steps.commit(plan)
    assert len(steps.history()) == 1
    assert plates.load(plan.outputs[0].plate_id).content_kind == "cells"


def test_a_step_is_planned_before_it_is_done(fresh):
    dna = _dna_plate(4)
    plan = steps.get("transfect").plan(dna)
    steps.commit(plan)
    assert steps.history().iloc[0]["status"] == "planned"
    steps.record(plan, operator="shivan", deviations="column 7 short")
    row = steps.history().iloc[0]
    assert row["status"] == "done"
    assert row["operator"] == "shivan"
    assert row["deviations"] == "column 7 short"
    assert row["finished_at"]


def test_a_human_step_records_exactly_like_a_robot_step(fresh):
    """The protocol must not fork on who did the work."""
    dna = _dna_plate(8)
    current = _run_chain(dna, ["transfect", "harvest_supernatant", "add_beads"])
    manual = steps.get("immunoprecipitate")
    assert manual.robot is False

    plan = manual.plan(current)
    steps.commit(plan)
    run = steps.record(plan, operator="shivan", lots=["LOT-beads-42"],
                       deviations="B03 lost beads")
    assert run.status == "done" and run.lots == ["LOT-beads-42"]
    stored = steps.history(plan.outputs[0].plate_id).iloc[0]
    assert json.loads(stored["inputs_json"]) == [current.plate_id]


def test_history_filters_to_one_plate(fresh):
    dna = _dna_plate(4)
    current = _run_chain(dna, ["transfect", "harvest_supernatant"])
    touching = steps.history(current.plate_id)
    assert len(touching) == 1
    assert touching.iloc[0]["step_type"] == "harvest_supernatant"


def test_reagent_lots_are_recorded_against_the_step(fresh):
    dna = _dna_plate(4)
    current = _run_chain(dna, ["transfect", "harvest_supernatant"])
    plan = steps.get("add_beads").plan(current)
    steps.commit(plan)
    steps.record(plan, lots=["LOT-beads-42", "LOT-buffer-7"])
    conn = stores.connect("assay")
    lots = conn.execute("SELECT lots_json FROM steps WHERE run_id = ?",
                        (plan.run_id,)).fetchone()[0]
    assert json.loads(lots) == ["LOT-beads-42", "LOT-buffer-7"]


# --- transfers --------------------------------------------------------------

def test_a_one_to_one_transfer_is_ninety_six_lines_at_one_volume(fresh):
    dna = _dna_plate()
    cells = _run_chain(dna, ["transfect"])
    plan = steps.get("harvest_supernatant").plan(cells, volume_ul=500.0)
    assert len(plan.transfers) == 96
    assert {t.volume_ul for t in plan.transfers} == {500.0}
    assert all(t.source_well == t.dest_well for t in plan.transfers)


def test_a_reagent_addition_moves_nothing_between_plates(fresh):
    """Beads come from a trough, so there is no source plate to transfer from."""
    dna = _dna_plate(8)
    current = _run_chain(dna, ["transfect", "harvest_supernatant"])
    plan = steps.get("add_beads").plan(current, bead_slurry_ul=25.0)
    assert plan.transfers == []
    assert any("trough" in line for line in plan.instructions)


def test_transfers_are_persisted_with_the_step(fresh):
    dna = _dna_plate(8)
    cells = _run_chain(dna, ["transfect"])
    plan = steps.get("harvest_supernatant").plan(cells)
    steps.commit(plan)
    conn = stores.connect("assay")
    n = conn.execute("SELECT COUNT(*) FROM transfers WHERE run_id = ?",
                     (plan.run_id,)).fetchone()[0]
    assert n == len(plan.transfers) == 8


# --- normalisation ----------------------------------------------------------

def test_normalisation_computes_a_different_volume_per_well(fresh):
    """The one step whose volumes differ per well. That is the hard capability."""
    dna = _dna_plate(8)
    wells = list(dna.wells["well"])[:8]
    concentrations = {w: c for w, c in zip(wells, [10, 20, 40, 50, 15, 30, 25, 60])}
    igg = _run_chain(dna, protocol.DEFAULT_CHAIN[:-1],
                     concentrations=concentrations)

    plan = steps.get("normalise").plan(igg, target_conc=10.0,
                                       target_volume_ul=100.0)
    volumes = {t.source_well: t.volume_ul for t in plan.transfers}
    assert len(set(volumes.values())) > 1, "uniform volumes would defeat the point"
    # c1v1 = c2v2, so a well at twice the target needs half the volume
    assert volumes[wells[0]] == pytest.approx(100.0)     # 10 ug/mL -> all of it
    assert volumes[wells[1]] == pytest.approx(50.0)      # 20 ug/mL -> half
    assert volumes[wells[2]] == pytest.approx(25.0)      # 40 ug/mL -> a quarter


def test_a_well_that_cannot_reach_the_target_is_flagged_not_silently_wrong(fresh):
    dna = _dna_plate(4)
    wells = list(dna.wells["well"])[:4]
    # 2 ug/mL cannot make 100 uL at 10 ug/mL however much you transfer
    igg = _run_chain(dna, protocol.DEFAULT_CHAIN[:-1],
                     concentrations={wells[0]: 2.0, wells[1]: 40.0,
                                     wells[2]: 40.0, wells[3]: 40.0})
    plan = steps.get("normalise").plan(igg, target_conc=10.0,
                                       target_volume_ul=100.0)
    out = plan.outputs[0]
    flagged = out.wells[out.wells["flags"].notna()]
    assert len(flagged) == 1 and flagged.iloc[0]["well"] == wells[0]
    assert "too dilute" in flagged.iloc[0]["flags"]
    assert plan.params["wells_flagged"] == 1


def test_a_well_with_no_measurement_is_skipped_rather_than_guessed(fresh):
    dna = _dna_plate(4)
    igg = _run_chain(dna, protocol.DEFAULT_CHAIN[:-1])   # no concentrations set
    plan = steps.get("normalise").plan(igg)
    assert plan.transfers == []
    assert plan.params["wells_flagged"] == 4


def test_bli_reads_the_plate_without_moving_anything(fresh):
    dna = _dna_plate(4)
    igg = _run_chain(dna, protocol.DEFAULT_CHAIN[:6])
    plan = steps.get("quantify_bli").plan(igg)
    assert plan.outputs == [] and plan.transfers == []
    assert plan.instrument == "bli"


# --- protocol defaults are placeholders, and say so -------------------------

def test_protocol_volumes_are_parameters_not_constants(fresh):
    """Q19: volumes and incubations are wet-lab decisions, not software ones."""
    dna = _dna_plate(4)
    cells = _run_chain(dna, ["transfect"])
    plan = steps.get("harvest_supernatant").plan(cells, volume_ul=250.0)
    assert {t.volume_ul for t in plan.transfers} == {250.0}
    assert plan.params["volume_ul"] == 250.0


def test_every_step_declares_its_assumed_defaults():
    for name in protocol.DEFAULT_CHAIN:
        step = steps.get(name)
        assert getattr(step, "assumed", None), \
            f"{name} must declare which numbers are placeholders"


# --- simulated ELISA --------------------------------------------------------

CLONES = [f"CLN-SIM-{i:03d}" for i in range(96)]


def test_a_pair_is_the_same_clones_in_the_same_wells():
    """The pairing is the whole point: a hit beats its OWN control well."""
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES, seed=3, scenario="nominal")
    assert list(pair.target["well"]) == list(pair.control["well"])
    assert list(pair.target["clone_id"]) == list(pair.control["clone_id"])
    assert len(pair.target) == len(CLONES)


def test_the_same_seed_gives_the_same_plates():
    from plateforge.assay import elisa

    a = elisa.simulate_pair(CLONES, seed=11, scenario="nominal")
    b = elisa.simulate_pair(CLONES, seed=11, scenario="nominal")
    assert a.target["od450"].tolist() == b.target["od450"].tolist()
    assert a.control["od450"].tolist() == b.control["od450"].tolist()


def test_binders_beat_their_own_control_well():
    from plateforge.assay import elisa

    truth = pd.concat([elisa.simulate_pair(CLONES, seed=s, scenario="nominal").truth()
                       for s in range(12)], ignore_index=True)
    binders = truth[truth["is_binder"] & ~truth["is_control_binder"]]
    others = truth[~truth["is_binder"] & ~truth["is_control_binder"]]
    assert binders["ratio"].median() > 3.0
    assert others["ratio"].median() < 2.0
    assert binders["od450_target"].median() > others["od450_target"].median() * 3


def test_some_wells_read_lower_than_their_control():
    """Non-binders are a coin flip between the two plates, and should be."""
    from plateforge.assay import elisa

    truth = elisa.simulate_pair(CLONES, seed=5, scenario="nominal").truth()
    flipped = truth[truth["delta"] < 0]
    assert 5 < len(flipped) < len(truth), "some, not none and not all"


def test_a_well_that_never_expressed_is_blank_on_both_plates():
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES, seed=9, scenario="nominal",
                               expression_failure_rate=0.5)
    truth = pair.truth()
    dead = truth[~truth["expressed"]]
    assert len(dead) > 10
    assert dead["od450_target"].max() < 0.3
    assert dead["od450_control"].max() < 0.3


def test_overdeveloped_tmb_saturates_instead_of_brightening():
    """The failure that destroys contrast: no headroom, so hits flatten."""
    from plateforge.assay import elisa

    clean = elisa.simulate_pair(CLONES, seed=4, scenario="nominal")
    burnt = elisa.simulate_pair(CLONES, seed=4, scenario="tmb_overdeveloped")
    assert burnt.target["od450"].median() > clean.target["od450"].median()
    assert burnt.control["od450"].median() > clean.control["od450"].median()

    # contrast between hits and background collapses
    def contrast(pair):
        t = pair.truth()
        hits = t[t["is_binder"]]["od450_target"]
        rest = t[~t["is_binder"]]["od450_target"]
        return hits.median() / max(rest.median(), 1e-6)

    if (clean.truth()["is_binder"]).sum() >= 3:
        assert contrast(burnt) < contrast(clean)


def test_readings_stay_inside_the_readers_range():
    from plateforge.assay import elisa

    for scenario in elisa.SCENARIOS:
        pair = elisa.simulate_pair(CLONES, seed=2, scenario=scenario)
        for frame in (pair.target, pair.control):
            assert frame["od450"].min() >= 0.0
            assert frame["od450"].max() <= 4.0 + 1e-9


def test_a_badly_washed_control_looks_like_a_target_plate():
    """The confusable case, generated on purpose."""
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES, seed=6, scenario="poor_wash_control")
    clean = elisa.simulate_pair(CLONES, seed=6, scenario="nominal")
    assert pair.control["od450"].median() > clean.control["od450"].median() * 1.5


def test_some_panels_yield_nothing_at_all():
    """A campaign with no hits is common, and a dataset without them is easy."""
    from plateforge.assay import elisa

    pairs = elisa.simulate_many(CLONES, 120, seed=2)
    counts = [sum(b.is_binder for b in p.behaviours.values()) for p in pairs]
    assert min(counts) == 0, "no zero-hit panel means the task is too easy"
    assert max(counts) > 15


def test_the_control_antigen_has_real_binders_too():
    """Without them, only the target plate can ever have a right tail."""
    from plateforge.assay import elisa

    pairs = elisa.simulate_many(CLONES, 60, seed=3)
    control_binders = sum(b.is_control_binder
                          for p in pairs for b in p.behaviours.values())
    assert control_binders > 0
    anti_tag = [b for p in pairs for b in p.behaviours.values()
                if b.is_binder and b.is_control_binder]
    assert anti_tag, "anti-tag clones light up both plates and must exist"


def test_plate_features_carry_no_label():
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    feats = elisa.plate_features(pair.target["od450"])
    assert feats and all(isinstance(v, float) for v in feats.values())
    for leak in ("is_target", "plate_kind", "scenario", "clone_id"):
        assert leak not in feats


def test_the_feature_table_has_one_row_per_plate():
    from plateforge.assay import elisa

    pairs = elisa.simulate_many(CLONES, 10, seed=4)
    table = elisa.feature_table(pairs)
    assert len(table) == 20
    assert table["is_target"].sum() == 10
    assert table.groupby("pair_id").size().eq(2).all()


def test_an_unknown_scenario_is_an_error_not_a_default():
    from plateforge.assay import elisa

    with pytest.raises(KeyError, match="unknown scenario"):
        elisa.simulate_pair(CLONES, seed=1, scenario="nope")


def test_the_heatmap_renders(tmp_path):
    from plateforge.assay import elisa

    pairs = elisa.simulate_many(CLONES, 4, seed=8)
    out = elisa.heatmap(pairs, tmp_path / "plates.png", n=4)
    assert out.exists() and out.stat().st_size > 10_000


# --- readings in the store --------------------------------------------------

def test_readings_round_trip_against_a_plate(fresh):
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    plate = plates.create("assay", 96, label="elisa", is_virtual=True,
                          notes="SIMULATED")
    plates.fill(plate, pair.target, volume_ul=100.0)
    plates.save(plate)
    n = plates.save_reads(plate.plate_id, pair.target, channel="od450",
                          units="OD450")
    assert n == 96

    back = plates.read_values(plate.plate_id, channel="od450")
    assert len(back) == 96
    assert set(back["clone_id"]) == set(CLONES), "readings join to the clones"
    assert back["value"].max() <= 4.0


# --- partial plates ---------------------------------------------------------

def test_a_plate_comes_back_with_every_well_even_when_half_used():
    """Because that is what a plate reader hands you."""
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES[:40], seed=4, scenario="nominal")
    assert len(pair.target) == 96
    assert (pair.target["role"] == "sample").sum() == 40
    assert (pair.target["role"] == "empty").sum() == 56
    assert pair.target[pair.target["role"] == "empty"]["clone_id"].isna().all()


def test_empty_wells_read_at_background_not_zero():
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES[:40], seed=4, scenario="nominal")
    blanks = pair.target[pair.target["role"] == "empty"]["od450"]
    assert blanks.median() > 0.0
    assert blanks.median() < 0.2
    assert blanks.std() > 0, "a blank is measured, not assumed"


def test_features_ignore_empty_wells():
    """The bug this prevents: blanks take over the median and every ratio.

    Measured on the same plate with 96 and 40 samples, before the fix:
    p99_over_median went 6.46 -> 15.04 for identical biology.
    """
    from plateforge.assay import elisa

    full = elisa.simulate_pair(CLONES, seed=4, scenario="nominal")
    part = elisa.simulate_pair(CLONES[:40], seed=4, scenario="nominal")

    a = elisa.plate_features(full.target["od450"], full.target["role"])
    b = elisa.plate_features(part.target["od450"], part.target["role"])
    naive = elisa.plate_features(part.target["od450"])

    assert abs(b["median"] - a["median"]) < 0.5 * a["median"]
    assert naive["median"] < 0.5 * a["median"], "the naive version IS broken"
    assert b["p99_over_median"] < naive["p99_over_median"] * 0.7


def test_features_report_how_full_the_plate_was():
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES[:24], seed=4, scenario="nominal")
    feats = elisa.plate_features(pair.target["od450"], pair.target["role"])
    assert feats["n_samples"] == 24
    assert feats["fill_fraction"] == pytest.approx(24 / 96)
    assert feats["blank_median"] > 0
    assert feats["signal_over_blank"] > 1


def test_a_full_plate_has_no_blank_features_and_says_so():
    """NaN by construction, not by accident: there were no blank wells."""
    import math
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES, seed=4, scenario="nominal")
    feats = elisa.plate_features(pair.target["od450"], pair.target["role"])
    assert feats["fill_fraction"] == 1.0
    assert math.isnan(feats["blank_median"])


def test_occupancy_is_exact_when_the_count_is_known():
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES[:40], seed=4, scenario="nominal")
    guess = elisa.infer_roles(pair.target["od450"], n_samples=40)
    agreement = (guess.values == pair.target["role"].values).mean()
    assert agreement > 0.8
    assert (guess == "sample").sum() == 40


def test_occupancy_declines_to_guess_without_a_count():
    """You cannot tell an empty well from a failed sample well by OD alone."""
    from plateforge.assay import elisa

    pair = elisa.simulate_pair(CLONES[:40], seed=4, scenario="nominal")
    guess = elisa.infer_roles(pair.target["od450"])
    assert set(guess) == {"sample"}, "no silent threshold guessing"


def test_fill_fraction_varies_across_a_simulated_campaign():
    from plateforge.assay import elisa

    pairs = elisa.simulate_many(CLONES, 60, seed=5)
    fills = {round(
        (p.target["role"] == "sample").mean(), 2) for p in pairs}
    assert len(fills) > 2, "a campaign of only full plates teaches the wrong thing"


# --- reading an export ------------------------------------------------------

def test_a_grid_is_found_by_shape_not_by_vendor_layout(tmp_path):
    from plateforge.assay import elisa, readers

    pairs = elisa.simulate_many(CLONES, 2, seed=9)
    sheets = []
    for i, pair in enumerate(pairs):
        sheets += [(f"Plate {2 * i + 1}", pair.target),
                   (f"Plate {2 * i + 2}", pair.control)]
    path = readers.write_gen5_like(tmp_path / "scan.xlsx", sheets)

    grids = readers.read(path, "biotek_gen5")
    assert len(grids) == 4
    assert all(g.plate_format == 96 for g in grids)
    assert all(g.n_values == 96 for g in grids)


def test_values_round_trip_through_a_workbook(tmp_path):
    from plateforge.assay import elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=3, scenario="nominal")
    path = readers.write_gen5_like(tmp_path / "one.xlsx",
                                   [("Plate 1", pair.target)])
    grid = readers.read(path, "biotek_gen5")[0]
    original = dict(zip(pair.target["well"], pair.target["od450"]))
    assert grid.values == pytest.approx(original)


def test_the_preamble_is_kept_and_is_not_full_of_nan(tmp_path):
    from plateforge.assay import elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=3, scenario="nominal")
    path = readers.write_gen5_like(tmp_path / "one.xlsx",
                                   [("Plate 7", pair.target)])
    grid = readers.read(path, "biotek_gen5")[0]
    assert any("Plate 7" in line for line in grid.preamble)
    assert not any("nan" in line.lower() for line in grid.preamble)


def test_sheet_order_is_preserved(tmp_path):
    """Scan order is the strongest structural clue for pairing."""
    from plateforge.assay import elisa, readers

    pairs = elisa.simulate_many(CLONES, 3, seed=2)
    sheets = []
    for i, pair in enumerate(pairs):
        sheets += [(f"T{i}", pair.target), (f"C{i}", pair.control)]
    path = readers.write_gen5_like(tmp_path / "s.xlsx", sheets)
    grids = readers.read(path, "biotek_gen5")
    assert [g.sheet for g in grids] == ["T0", "C0", "T1", "C1", "T2", "C2"]


def test_a_grid_that_does_not_start_at_the_corner_is_still_found(tmp_path):
    from plateforge.assay import elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=3, scenario="nominal")
    path = readers.write_gen5_like(tmp_path / "off.xlsx",
                                   [("P", pair.target)], preamble=True)
    grids = readers.read(path, "biotek_gen5")
    assert grids and grids[0].header_row > 0


def test_readers_declare_their_confidence():
    from plateforge.assay import readers

    table = readers.available()
    assert set(table["confidence"]) <= {"verified", "heuristic", "sketch"}
    ours = table[table["reader"] == "long_csv"].iloc[0]
    assert ours["confidence"] == "verified", "our own format is not a guess"
    gen5 = table[table["reader"] == "biotek_gen5"].iloc[0]
    assert gen5["confidence"] == "heuristic"


# --- assigning grids to plates ---------------------------------------------

def _trained_model():
    from plateforge.assay import assign, elisa
    return assign.train(elisa.simulate_many(CLONES, 80, seed=1))


def test_an_explicit_choice_beats_everything_else(fresh):
    from plateforge.assay import assign, elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    path = readers.write_gen5_like(fresh / "s.xlsx",
                                   [("PF-001", pair.target),
                                    ("PF-002", pair.control)])
    grids = readers.read(path)
    proposals = assign.propose(grids, explicit={0: "target"},
                               barcodes={"PF-001": "PLT-x"})
    assert proposals[0].source == "explicit"
    assert proposals[0].confidence == 1.0
    assert not proposals[0].needs_confirmation


def test_a_barcode_in_the_sheet_is_matched(fresh):
    from plateforge.assay import assign, elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    path = readers.write_gen5_like(fresh / "s.xlsx",
                                   [("PF-SUP-0042", pair.target)])
    grids = readers.read(path)
    proposals = assign.propose(grids, barcodes={"PF-SUP-0042": "PLT-abc"})
    assert proposals[0].plate_id == "PLT-abc"
    assert proposals[0].source == "barcode"
    assert "PF-SUP-0042" in proposals[0].reason


def test_two_barcodes_in_one_sheet_is_unresolved_not_a_coin_flip(fresh):
    from plateforge.assay import assign, elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    path = readers.write_gen5_like(fresh / "s.xlsx",
                                   [("PF-A1 and PF-B2", pair.target)])
    grids = readers.read(path)
    proposals = assign.propose(grids, barcodes={"PF-A1": "PLT-a",
                                                "PF-B2": "PLT-b"})
    assert proposals[0].source == "unresolved"
    assert proposals[0].needs_confirmation
    assert set(proposals[0].alternatives) == {"PLT-a", "PLT-b"}


def test_the_sheet_often_just_says_which_it_is(fresh):
    from plateforge.assay import assign, elisa, readers

    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    path = readers.write_gen5_like(
        fresh / "s.xlsx", [("target antigen", pair.target),
                           ("HLA control", pair.control)])
    grids = readers.read(path)
    proposals = assign.propose(grids)
    assert proposals[0].plate_kind == "target"
    assert proposals[1].plate_kind == "control"


def test_a_model_guess_always_needs_confirmation(fresh):
    from plateforge.assay import assign, elisa, readers

    model = _trained_model()
    pair = elisa.simulate_pair(CLONES, seed=31, scenario="nominal")
    path = readers.write_gen5_like(fresh / "s.xlsx",
                                   [("sheet one", pair.target),
                                    ("sheet two", pair.control)])
    grids = readers.read(path)
    proposals = assign.propose(grids, model=model)
    for proposal in proposals:
        assert proposal.source == "model"
        assert proposal.needs_confirmation, "a guess is never a decision"
        assert "p(target)" in proposal.reason
        assert 0.5 <= proposal.confidence <= 1.0


def test_a_guess_says_when_it_did_not_know_the_plate_map(fresh):
    from plateforge.assay import assign, elisa, readers

    model = _trained_model()
    pair = elisa.simulate_pair(CLONES[:40], seed=31, scenario="nominal")
    path = readers.write_gen5_like(fresh / "s.xlsx", [("s", pair.target)])
    grids = readers.read(path)
    proposals = assign.propose(grids, model=model)
    assert "PLATE MAP UNKNOWN" in proposals[0].reason

    told = assign.propose(grids, model=model, n_samples=40)
    assert "n_samples" in told[0].reason


def test_a_model_round_trips_through_disk(fresh):
    from plateforge.assay import assign

    model = _trained_model()
    directory = assign.save(model, label="test-model")
    assert (directory / "model.joblib").exists()
    meta = json.loads((directory / "model.json").read_text())
    assert meta["kind"] == "plate_target_vs_control"
    assert "SIMULATED" in json.dumps(meta["trained_on"])
    assert "re-fit on real data" in meta["warning"]

    back = assign.load(directory.name)
    assert back.features == model.features
    assert back.accuracy == pytest.approx(model.accuracy)
    assert assign.latest() == directory.name


def test_the_report_names_what_still_needs_a_human(fresh):
    from plateforge.assay import assign, elisa, readers

    model = _trained_model()
    pair = elisa.simulate_pair(CLONES, seed=1, scenario="nominal")
    path = readers.write_gen5_like(fresh / "s.xlsx",
                                   [("PF-KNOWN", pair.target), ("?", pair.control)])
    grids = readers.read(path)
    proposals = assign.propose(grids, barcodes={"PF-KNOWN": "PLT-a"},
                               model=model)
    table = assign.report(proposals)
    assert list(table.columns)[:4] == ["grid_index", "sheet", "plate_id",
                                       "plate_kind"]
    assert len(assign.unresolved(proposals)) >= 1
