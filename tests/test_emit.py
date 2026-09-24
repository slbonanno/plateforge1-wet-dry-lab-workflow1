"""Worklists, and the lab that decides which one to write.

Rule 6 is the thing under test here as much as any format detail: an emitter
that has never seen a real file from its instrument must say so, loudly,
every time.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from plateforge.core import lab, stores
from plateforge.emit import worklists


@pytest.fixture
def fresh(tmp_path, monkeypatch):
    monkeypatch.setenv("PLATEFORGE_DATA", str(tmp_path))
    monkeypatch.delenv("PLATEFORGE_LAB", raising=False)
    stores.close_all()
    yield tmp_path
    stores.close_all()


@pytest.fixture
def transfers():
    return pd.DataFrame({
        "source_plate_id": ["P1"] * 4,
        "source_well": ["A01", "B01", "A02", "H12"],
        "dest_plate_id": ["P2"] * 4,
        "dest_well": ["A01", "B01", "A02", "H12"],
        "volume_ul": [500.0, 250.5, 33.3, 12.0],
    })


# --- rule 6 -----------------------------------------------------------------

def test_an_unverified_format_is_never_emitted_by_accident(transfers):
    with pytest.raises(ValueError, match="confidence is 'documented'"):
        worklists.emit(transfers, "tecan_evo")
    got = worklists.emit(transfers, "tecan_evo", allow_unverified=True)
    assert got.confidence == "documented" and got.caveats


def test_our_own_format_needs_no_permission(transfers):
    """It cannot be wrong about anyone else's file."""
    got = worklists.emit(transfers, "generic")
    assert got.verified and got.caveats == []


def test_nothing_is_verified_yet():
    """The day a fixture lands, this test changes. Until then it is the truth."""
    table = worklists.available()
    verified = set(table[table["confidence"] == "verified"]["instrument"])
    assert verified == {"generic"}, \
        "an instrument format is 'verified' only once a real file is in fixtures/"


def test_every_emitter_declares_its_confidence():
    table = worklists.available()
    assert set(table["confidence"]) <= {"verified", "documented", "sketch"}
    assert table["software"].notna().all()


def test_an_unverified_file_is_written_with_a_caveat_beside_it(tmp_path, transfers):
    got = worklists.emit(transfers, "tecan_evo", allow_unverified=True)
    path = got.write(tmp_path, stem="harvest")
    caveat = path.with_suffix(".CAVEAT.txt")
    assert path.exists() and caveat.exists()
    assert "documented" in caveat.read_text()
    assert "robotools" in caveat.read_text()


def test_a_verified_file_gets_no_caveat(tmp_path, transfers):
    path = worklists.emit(transfers, "generic").write(tmp_path)
    assert not path.with_suffix(".CAVEAT.txt").exists()


# --- Tecan ------------------------------------------------------------------

def test_tecan_positions_run_down_each_column():
    """A01=1, B01=2, A02=9 on a 96-well plate. Getting this wrong transposes."""
    assert worklists.tecan_position("A01") == 1
    assert worklists.tecan_position("B01") == 2
    assert worklists.tecan_position("H01") == 8
    assert worklists.tecan_position("A02") == 9
    assert worklists.tecan_position("H12") == 96
    assert worklists.tecan_position("a1") == 1, "normalised on the way in"


def test_tecan_emits_aspirate_dispense_wash_per_transfer(transfers):
    got = worklists.emit(transfers, "tecan_evo", allow_unverified=True)
    lines = got.text.strip().splitlines()
    assert len(lines) == len(transfers) * 3
    assert lines[0].startswith("A;") and lines[1].startswith("D;")
    assert lines[2] == "W1;"
    # ten fields after the letter
    assert lines[0].count(";") == 10


def test_tecan_volumes_are_two_decimal_places(transfers):
    got = worklists.emit(transfers, "tecan_evo", allow_unverified=True)
    assert ";500.00;" in got.text and ";250.50;" in got.text


def test_tecan_flags_volumes_the_instrument_cannot_do(transfers):
    got = worklists.emit(transfers, "tecan_evo", allow_unverified=True,
                         min_ul=50.0, max_ul=400.0)
    assert any("below" in i for i in got.issues)
    assert any("above" in i for i in got.issues)


def test_wash_can_be_turned_off(transfers):
    got = worklists.emit(transfers, "tecan_evo", allow_unverified=True,
                         wash_between=False)
    assert "W1;" not in got.text


# --- Opentrons --------------------------------------------------------------

def test_opentrons_emits_a_runnable_shaped_protocol(transfers):
    got = worklists.emit(transfers, "opentrons", allow_unverified=True)
    assert got.suffix == ".py"
    assert "def run(protocol: protocol_api.ProtocolContext):" in got.text
    assert "apiLevel" in got.text
    assert got.text.count("pipette.transfer(") == len(transfers)
    compile(got.text, "protocol.py", "exec")       # it must at least parse


def test_opentrons_says_its_deck_is_a_guess(transfers):
    got = worklists.emit(transfers, "opentrons", allow_unverified=True)
    assert any("deck" in c.lower() for c in got.caveats)


# --- the ones that refuse ---------------------------------------------------

def test_integra_refuses_and_says_what_would_unblock_it(transfers):
    with pytest.raises(NotImplementedError, match="fixtures/integra/"):
        worklists.emit(transfers, "integra_vialab", allow_unverified=True)


@pytest.mark.parametrize("name", ["hamilton_venus", "beckman_biomek",
                                  "agilent_bravo"])
def test_a_closed_format_refuses_rather_than_guessing(transfers, name):
    with pytest.raises(NotImplementedError):
        worklists.emit(transfers, name, allow_unverified=True)


def test_an_empty_worklist_is_reported_not_silently_written():
    empty = pd.DataFrame(columns=["source_well", "dest_well", "volume_ul"])
    got = worklists.emit(empty, "tecan_evo", allow_unverified=True)
    assert got.lines == 0
    assert "no transfers" in " ".join(got.issues)


# --- the lab ----------------------------------------------------------------

def test_an_unconfigured_lab_still_runs(fresh):
    here = lab.load()
    assert not here.configured
    assert here.liquid_handler == "generic"
    assert "not configured" in here.describe()


def test_a_written_lab_config_is_read_back(fresh):
    lab.write({"name": "Bonanno lab", "liquid_handler": "integra_vialab",
               "liquid_handler_model": "ASSIST PLUS",
               "has_single_channel_module": True, "bli": "octet"})
    here = lab.load()
    assert here.configured
    assert here.name == "Bonanno lab"
    assert here.liquid_handler == "integra_vialab"
    assert here.has_single_channel_module is True
    assert "ASSIST PLUS" in here.describe()


def test_the_config_is_json_with_a_comment_key(fresh):
    """Rule 8: JSON has no comments, so a _comment key carries the intent."""
    written = lab.write({"name": "x", "liquid_handler": "generic"})
    data = json.loads(written.read_text())
    assert "_comment" in data
    assert data["name"] == "x"


def test_a_run_may_override_the_lab_and_the_override_is_visible(fresh):
    lab.write({"name": "Bonanno lab", "liquid_handler": "integra_vialab"})
    borrowed = lab.load({"liquid_handler": "opentrons"})
    assert borrowed.liquid_handler == "opentrons"
    assert "overrides" in borrowed.source, "an exception must be visible as one"


def test_instrument_for_prefers_an_explicit_choice(fresh):
    lab.write({"name": "x", "liquid_handler": "tecan_evo"})
    assert lab.instrument_for("liquid_handler") == "tecan_evo"
    assert lab.instrument_for("liquid_handler", override="opentrons") == "opentrons"


def test_asking_for_an_instrument_the_lab_lacks_is_an_error(fresh):
    lab.write({"name": "x", "liquid_handler": "generic"})
    with pytest.raises(ValueError, match="no bli configured"):
        lab.instrument_for("bli")


def test_unknown_config_keys_are_kept_not_dropped(fresh):
    """Rule 7: a lab will have something this schema did not anticipate."""
    lab.write({"name": "x", "liquid_handler": "generic",
               "centrifuge": "Eppendorf 5810R"})
    here = lab.load()
    assert here.extra["centrifuge"] == "Eppendorf 5810R"
    assert here.as_dict()["centrifuge"] == "Eppendorf 5810R"


def test_a_broken_config_says_so_rather_than_falling_back(fresh):
    lab.path().write_text("{not json")
    with pytest.raises(ValueError, match="not readable JSON"):
        lab.load()


# --- steps to worklists -----------------------------------------------------

def test_a_planned_step_can_be_emitted_for_an_instrument(fresh):
    """The join: assay plans the transfers, emit writes the file."""
    from plateforge.assay import plates, steps

    plate = plates.create("cells", 96, label="culture")
    rows = pd.DataFrame({"well": ["A01", "B01", "A02"],
                         "clone_id": ["CLN-1", "CLN-2", "CLN-3"]})
    plates.fill(plate, rows, volume_ul=1000.0)
    plates.save(plate)

    plan = steps.get("harvest_supernatant").plan(plate, volume_ul=500.0)
    got = worklists.emit(plan.transfer_frame, "tecan_evo", allow_unverified=True)
    assert got.lines == len(plan.transfers) * 3
    assert ";500.00;" in got.text
