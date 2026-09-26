"""Sequencing: ordering it, reading what comes back, and deciding what is in
each well.

Everything about the *formats* here is checked against real files in
`fixtures/` (rule 6): three Azenta/GENEWIZ `.ab1` traces, the pcDNA3.1(+)
SnapGene map, and the three vendors' own order-form templates. Everything
about the *logic* is checked against generated traces, because covering nine
verdicts needs defects on demand and a real plate cannot supply them.
"""
from __future__ import annotations

import random
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from plateforge.core import wells
from plateforge.library import abif, codon, dna, sanger, snapgene
from plateforge.library import vector as vectors

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
TRACES = FIXTURES / "sanger"
BACKBONE = FIXTURES / "backbone" / "pcDNA3.1.dna"
FORMS = FIXTURES / "sequencing"

pytestmark = pytest.mark.skipif(not TRACES.exists(),
                                reason="real fixture files are not present")


# --- reading a real trace ----------------------------------------------------

def test_a_real_trace_says_where_it_came_from():
    """The whole rerun design rests on this: the instrument writes the plate,
    the well and the run time into the file, so none of it has to be guessed
    from a filename."""
    trace = abif.read(TRACES / "194-M13F.ab1")
    assert trace.sample == "194-M13F"
    assert trace.well == "F5"                      # vendors do not zero-pad
    assert trace.plate_id == "04A000129809"
    assert trace.instrument == "3730"
    assert trace.vendor == "Azenta/GENEWIZ"
    assert trace.run_started == datetime(2024, 2, 8, 0, 57, 21)


def test_quality_is_read_as_scores_not_as_text():
    """PCON is a char array whose bytes are Phred scores. Decoding it as a
    string gives every base a plausible-looking quality in the 30s-60s and is
    completely meaningless."""
    trace = abif.read(TRACES / "194-M13F.ab1")
    assert len(trace.quality) == len(trace.sequence) == 929
    assert max(trace.quality) <= 62                # Phred, not ASCII
    assert trace.bases_at_least(20) == 828


def test_the_three_fixtures_are_one_good_pair_and_one_failure():
    good = abif.read(TRACES / "194-M13F.ab1")
    failed = abif.read(TRACES / "194-3xP3.ab1")
    assert good.mean_quality > 40
    assert failed.mean_quality < 15
    assert failed.bases_at_least(20) < 100


def test_a_failed_read_trims_to_nothing_rather_than_to_garbage():
    failed = abif.read(TRACES / "194-3xP3.ab1")
    assert abif.mott_trim(failed.quality) == (0, 0)


def test_trimming_removes_the_called_ns_at_the_start():
    trace = abif.read(TRACES / "194-M13F.ab1")
    start, end = abif.mott_trim(trace.quality)
    assert start > 0 and end > start
    assert "N" not in trace.sequence[start:start + 50]
    assert trace.sequence[:start].count("N") > 5


def test_a_file_that_is_not_abif_is_refused_by_name():
    with pytest.raises(abif.NotABIF):
        abif.read(BACKBONE)


# --- the strongest check available: two real reads of one molecule ----------

def test_the_forward_and_reverse_reads_of_one_clone_agree_exactly():
    """194-M13F and 194-M13R are the same insert read from both ends. Their
    overlap has to be a perfect reverse complement, and nothing short of a
    correct parser, quality decode, trim, reverse complement and aligner will
    produce that. This one test covers the whole chain on real data."""
    forward, reverse = (abif.read(TRACES / f"194-{p}.ab1") for p in ("M13F", "M13R"))
    f_start, f_end = abif.mott_trim(forward.quality)
    r_start, r_end = abif.mott_trim(reverse.quality)

    aligned = dna.align_to(dna.reverse_complement(reverse.sequence[r_start:r_end]),
                           forward.sequence[f_start:f_end])
    assert aligned is not None
    assert aligned.length > 200
    assert aligned.identity == 1.0
    assert aligned.mismatches == 0 and aligned.gaps == 0


# --- writing a trace, so the logic can be tested ----------------------------

def test_a_written_trace_reads_back_through_the_same_parser(tmp_path):
    written = abif.write(tmp_path / "t.ab1", "ACGT" * 50, [38] * 200,
                         sample="CLN-042", well="B07", plate_id="PF-1",
                         run_started=datetime(2026, 3, 4, 9, 30, 15))
    back = abif.read(written)
    assert back.sequence == "ACGT" * 50
    assert back.quality == [38] * 200
    assert back.sample == "CLN-042" and back.well == "B07"
    assert back.plate_id == "PF-1"
    assert back.run_started == datetime(2026, 3, 4, 9, 30, 15)


# --- the backbone ------------------------------------------------------------

def test_the_backbone_is_the_plasmid_it_claims_to_be():
    plasmid = snapgene.read(BACKBONE)
    assert len(plasmid) == 5428
    assert plasmid.circular
    assert snapgene.check(plasmid, 5428,
                          ("CMV promoter", "MCS", "bGH poly(A) signal",
                           "AmpR", "NeoR/KanR", "f1 ori")) == []


def test_the_mcs_is_where_the_map_says_it_is():
    plasmid = snapgene.read(BACKBONE)
    mcs = plasmid.feature("MCS")
    sequence = mcs.sequence(plasmid)
    assert sequence.startswith("GCTAGC")            # NheI
    assert "AAGCTT" in sequence                     # HindIII
    assert "GCGGCCGC" in sequence                   # NotI
    assert plasmid.feature("CMV promoter").end < mcs.start


def test_a_slice_wraps_the_origin_on_a_circular_plasmid():
    """A primer near the end of the numbering reads straight through base 1."""
    plasmid = snapgene.read(BACKBONE)
    wrapped = plasmid.slice(len(plasmid) - 10, len(plasmid) + 10)
    assert len(wrapped) == 20
    assert wrapped[:10] == plasmid.sequence[-10:]
    assert wrapped[10:] == plasmid.sequence[:10]


def test_a_backbone_of_the_wrong_length_is_reported_not_used():
    plasmid = snapgene.read(BACKBONE)
    problems = snapgene.check(plasmid, 5429, ("Nonexistent feature",))
    assert any("5429" in p for p in problems)
    assert any("Nonexistent" in p for p in problems)


# --- the aligner -------------------------------------------------------------

@pytest.fixture(scope="module")
def reference() -> str:
    return snapgene.read(BACKBONE).sequence[:1500]


def _mutate(seq, subs=0, delete=None, insert=None, seed=1):
    rng = random.Random(seed)
    chars = list(seq)
    for _ in range(subs):
        at = rng.randrange(len(chars))
        chars[at] = rng.choice([b for b in "ACGT" if b != chars[at]])
    out = "".join(chars)
    if delete:
        at, n = delete
        out = out[:at] + out[at + n:]
    if insert:
        at, text = insert
        out = out[:at] + text + out[at:]
    return out


def test_an_exact_read_aligns_with_no_differences(reference):
    aligned = dna.align_to(reference, reference[300:900])
    assert aligned.identity == 1.0
    assert dna.differences(aligned) == []
    assert (aligned.ref_start, aligned.ref_end) == (300, 900)


def test_junk_at_the_start_of_a_read_overhangs_rather_than_aligning(reference):
    """A basecaller emits noise before the template starts. Scoring it as gaps
    turns 40 bases of nothing into a page of invented indels."""
    junk = "N" * 13 + "TTGACCATAGGCTAGCCATAGGCTA"
    aligned = dna.align_to(reference, junk + reference[300:900])
    # The junk overhangs rather than being aligned. A couple of its last bases
    # may match the reference by chance -- that is the aligner being right,
    # not wrong -- so what is asserted is the property: no gaps were invented,
    # the read ends where the template does, and nearly all the junk is
    # outside the alignment.
    assert aligned.gaps == 0
    assert aligned.ref_end == 900
    assert aligned.query_start >= len(junk) - 5
    assert aligned.ref_start >= 300 - 5
    assert not [d for d in dna.differences(aligned) if d["kind"] != "substitution"]


@pytest.mark.parametrize("label,kwargs,kind,length", [
    ("one base deleted", {"delete": (200, 1)}, "deletion", 1),
    ("three deleted", {"delete": (200, 3)}, "deletion", 3),
    ("six deleted", {"delete": (200, 6)}, "deletion", 6),
    ("twelve deleted", {"delete": (300, 12)}, "deletion", 12),
    ("one inserted", {"insert": (150, "G")}, "insertion", 1),
    ("nine inserted", {"insert": (150, "GGCATTACG")}, "insertion", 9),
])
def test_an_indel_is_reported_as_one_event(reference, label, kwargs, kind, length):
    """Affine gaps with a state-aware traceback. Without either, one 6 bp
    deletion comes back as two or three smaller ones -- same score, wrong
    biology, and it changes an in-frame deletion into an apparent frameshift."""
    inner = reference[300:900]
    aligned = dna.align_to(reference, _mutate(inner, **kwargs))
    found = dna.differences(aligned)
    assert len(found) == 1, f"{label}: {found}"
    assert found[0]["kind"] == kind
    assert found[0]["length"] == length


def test_substitutions_are_found_at_the_right_positions(reference):
    inner = reference[300:900]
    changed = list(inner)
    for at in (57, 86, 173):
        changed[at] = "A" if changed[at] != "A" else "C"
    aligned = dna.align_to(reference, "".join(changed))
    found = dna.differences(aligned)
    assert [d["kind"] for d in found] == ["substitution"] * 3
    assert [d["ref_start"] for d in found] == [357, 386, 473]


def test_the_screen_finds_the_right_reference_among_decoys(reference):
    rng = random.Random(3)
    panel = {"correct": reference}
    for i in range(95):
        panel[f"decoy{i}"] = "".join(rng.choice("ACGT") for _ in range(1200))
    ranked = dna.screen(reference[300:900], panel)
    assert ranked[0].name == "correct"
    assert ranked[0].score > 20 * ranked[1].score


def test_the_screen_notices_a_read_on_the_reverse_strand(reference):
    panel = {"correct": reference}
    ranked = dna.screen(dna.reverse_complement(reference[300:900]), panel)
    assert ranked[0].name == "correct"
    assert ranked[0].forward is False


# --- verdicts ----------------------------------------------------------------

def _vh(seed: int) -> str:
    rng = random.Random(seed)
    protein = ("EVQLVESGGGLVQPGGSLRLSCAAS"
               + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(95)))
    return codon.optimize(protein).dna


@pytest.fixture(scope="module")
def panel() -> dict:
    vector = vectors.get("pcdna-igg1-dual-vk-v1")
    frame = pd.DataFrame({
        "clone_id": [f"CLN-{i:03d}" for i in range(12)],
        "well": [wells.from_rc(r, c) for c in range(1, 3) for r in range(1, 9)][:12],
        "cds_dna": [_vh(i) for i in range(12)],
    })
    return sanger.expectations_from(frame, vector)


def _trace(tmp_path, expectation, sequence, name, *, quality=45, junk=35):
    """A read of `sequence` with a low-quality 5' run, as a real trace has."""
    rng = random.Random(len(name))
    noise = "".join(rng.choice("ACGT") for _ in range(junk))
    return abif.write(tmp_path / f"{name}.ab1", noise + sequence,
                      [6] * len(noise) + [quality] * len(sequence),
                      sample=expectation.clone_id, well=expectation.well,
                      plate_id="PF-SEQ-1",
                      run_started=datetime(2026, 2, 1, 9, 0, 0))


def _window(expectation, sequence=None, pad=30) -> str:
    text = sequence if sequence is not None else expectation.context
    start = max(expectation.insert_at - pad, 0)
    return text[start:expectation.insert_at + len(expectation.insert) + pad]


def _call(tmp_path, panel, expectation, sequence, name, **kwargs):
    path = _trace(tmp_path, expectation, sequence, name, **kwargs)
    others = {e.clone_id: e.insert for e in panel.values()
              if e.clone_id != expectation.clone_id}
    return sanger.call_one(abif.read(path), expectation, others)


def test_a_perfect_read_is_exact(tmp_path, panel):
    expectation = panel["A01"]
    call = _call(tmp_path, panel, expectation, _window(expectation), "exact")
    assert call.verdict == sanger.EXACT
    assert call.usable
    assert call.identity == 1.0
    assert call.differences == []


def test_a_silent_change_is_usable_and_says_so(tmp_path, panel):
    expectation = panel["B01"]
    context = list(expectation.context)
    at = expectation.insert_at + 5                  # third base of a codon
    for alternative in "ACGT":
        if alternative == context[at]:
            continue
        context[at] = alternative
        insert = "".join(context)[expectation.insert_at:
                                  expectation.insert_at + len(expectation.insert)]
        if codon.translate(insert).rstrip("*") == expectation.protein:
            break
    call = _call(tmp_path, panel, expectation,
                 _window(expectation, "".join(context)), "silent")
    assert call.verdict == sanger.SILENT
    assert call.usable
    assert len(call.differences) == 1


def test_a_missense_change_names_the_residue(tmp_path, panel):
    expectation = panel["C01"]
    context = list(expectation.context)
    at = expectation.insert_at + 30
    context[at] = "A" if context[at] != "A" else "C"
    call = _call(tmp_path, panel, expectation,
                 _window(expectation, "".join(context)), "missense")
    assert call.verdict == sanger.MISSENSE
    assert not call.usable
    assert call.protein_changes
    assert call.protein_changes[0][0].isalpha()     # e.g. "L11M"


def test_a_one_base_deletion_is_a_frameshift(tmp_path, panel):
    expectation = panel["D01"]
    at = expectation.insert_at + 50
    broken = expectation.context[:at] + expectation.context[at + 1:]
    call = _call(tmp_path, panel, expectation,
                 _window(expectation, broken), "frameshift")
    assert call.verdict == sanger.FRAMESHIFT
    assert "frame" in " ".join(call.notes)


def test_a_three_base_deletion_is_an_indel_not_a_frameshift(tmp_path, panel):
    expectation = panel["E01"]
    at = expectation.insert_at + 60
    shortened = expectation.context[:at] + expectation.context[at + 3:]
    call = _call(tmp_path, panel, expectation,
                 _window(expectation, shortened), "inframe")
    assert call.verdict == sanger.INDEL
    assert "in frame" in " ".join(call.notes)


def test_empty_backbone_is_recognised_as_such(tmp_path, panel):
    expectation = panel["F01"]
    without = (expectation.context[:expectation.insert_at]
               + expectation.context[expectation.insert_at + len(expectation.insert):])
    call = _call(tmp_path, panel, expectation, without, "empty")
    assert call.verdict == sanger.NO_INSERT
    # and it does not claim to have matched some other clone
    assert call.matched == ""
    assert call.confidence == 0.0


def test_a_swapped_well_names_the_clone_it_actually_holds(tmp_path, panel):
    """The expensive failure. Every later step would carry the wrong antibody
    under the right name."""
    expectation, actual = panel["G01"], panel["B02"]
    call = _call(tmp_path, panel, expectation, _window(actual), "swap")
    assert call.verdict == sanger.WRONG_CLONE
    assert call.matched == actual.clone_id
    assert call.identity > 0.95                     # it IS that clone, cleanly
    assert actual.clone_id in " ".join(call.notes)


def test_two_colonies_in_one_well_are_called_mixed(tmp_path, panel):
    expectation = panel["H01"]
    window = list(_window(expectation))
    for at in range(20, len(window), 25):
        window[at] = "R"                            # A or G, both present
    call = _call(tmp_path, panel, expectation, "".join(window), "mixed")
    assert call.verdict == sanger.MIXED
    assert "two templates" in " ".join(call.notes)


def test_a_dead_trace_is_low_quality_and_nothing_is_read_into_it(tmp_path, panel):
    expectation = panel["A02"]
    call = _call(tmp_path, panel, expectation, _window(expectation), "dead",
                 quality=6)
    assert call.verdict == sanger.LOW_QUALITY
    assert call.differences == []
    assert "Re-run" in " ".join(call.notes)


def test_a_partial_read_is_not_called(tmp_path, panel):
    """Half an insert cannot say the other half is correct."""
    expectation = panel["B02"]
    half = expectation.context[expectation.insert_at:
                               expectation.insert_at + len(expectation.insert) // 2]
    call = _call(tmp_path, panel, expectation, half, "partial")
    assert call.verdict == sanger.LOW_QUALITY
    assert call.coverage < 0.8


def test_a_whole_plate_is_summarised_worst_news_first(tmp_path, panel):
    traces = []
    for index, (well, expectation) in enumerate(sorted(panel.items())):
        sequence = _window(expectation)
        if index == 3:                              # one broken well
            at = expectation.insert_at + 40
            sequence = _window(expectation,
                               expectation.context[:at] + expectation.context[at + 1:])
        traces.append(abif.read(_trace(tmp_path, expectation, sequence, f"w{index}")))

    calls, report = sanger.verify(traces, panel)
    assert len(calls) == len(panel)
    assert calls.iloc[0]["verdict"] == sanger.EXACT        # sorted best first
    assert sanger.FRAMESHIFT in set(calls["verdict"])
    assert len(sanger.usable_wells(calls)) == len(panel) - 1
    assert report is not None and not report.superseded

    counts = sanger.summarise(calls)
    assert counts["wells"].sum() == len(panel)
    assert bool(counts[counts["verdict"] == sanger.EXACT]["usable"].iloc[0])


# --- reruns ------------------------------------------------------------------

@pytest.mark.parametrize("name,stem,marked", [
    ("194-3xP3_RR", "194-3xP3", True),
    ("194-M13R-rerun", "194-M13R", True),
    ("CLN-012_redo", "CLN-012", True),
    ("well-A01-run2", "well-A01", True),
    ("sample (2)", "sample", True),
    ("CLN-004", "CLN-004", False),
    ("PLATE_2", "PLATE_2", False),        # a plate, not a second attempt
])
def test_rerun_markers_are_stripped_conservatively(name, stem, marked):
    assert sanger.stem_of(name) == (stem, marked)


def test_a_rerun_supersedes_its_original_and_nothing_is_deleted(tmp_path):
    original = abif.write(tmp_path / "CLN-007.ab1", "ACGT" * 150, [40] * 600,
                          sample="CLN-007", well="C02", plate_id="PF-1",
                          run_started=datetime(2026, 2, 1, 9, 0, 0))
    repeat = abif.write(tmp_path / "CLN-007_RR.ab1", "ACGT" * 200, [44] * 800,
                        sample="CLN-007_RR", well="H11", plate_id="PF-2",
                        run_started=datetime(2026, 2, 8, 14, 0, 0))
    report = sanger.resolve_reruns([abif.read(original), abif.read(repeat)])

    assert [t.path.name for t in report.keep] == ["CLN-007_RR.ab1"]
    assert len(report.superseded) == 1
    dropped, winner = report.superseded[0]
    assert dropped.path.name == "CLN-007.ab1"
    assert winner.path.name == "CLN-007_RR.ab1"
    assert len(report.table) == 2


def test_a_rerun_lands_on_a_different_plate_and_well_and_is_still_grouped(tmp_path):
    """Which is why the sample name is the key, not the plate and well: a
    rerun is cherry-picked into a fresh plate."""
    a = abif.write(tmp_path / "x.ab1", "ACGT" * 100, [40] * 400, sample="CLN-009",
                   well="A01", plate_id="PLATE-A",
                   run_started=datetime(2026, 2, 1, 9, 0, 0))
    b = abif.write(tmp_path / "y.ab1", "ACGT" * 100, [40] * 400,
                   sample="CLN-009 (2)", well="G11", plate_id="PLATE-B",
                   run_started=datetime(2026, 2, 3, 9, 0, 0))
    report = sanger.resolve_reruns([abif.read(a), abif.read(b)])
    assert len(report.keep) == 1 and len(report.superseded) == 1
    assert report.keep[0].path.name == "y.ab1"


def test_a_worse_rerun_is_flagged_rather_than_trusted_silently(tmp_path):
    """A rerun can fail too. Taking the later run is a default, not a fact."""
    good = abif.write(tmp_path / "g.ab1", "ACGT" * 150, [44] * 600,
                      sample="CLN-003", well="A01", plate_id="P1",
                      run_started=datetime(2026, 2, 1, 9, 0, 0))
    worse = abif.write(tmp_path / "w.ab1", "ACGT" * 150, [8] * 600,
                       sample="CLN-003_RR", well="A01", plate_id="P2",
                       run_started=datetime(2026, 2, 5, 9, 0, 0))
    report = sanger.resolve_reruns([abif.read(good), abif.read(worse)])
    assert report.keep[0].path.name == "w.ab1"
    assert "WARNING" in " ".join(report.table["why"])

    by_quality = sanger.resolve_reruns([abif.read(good), abif.read(worse)],
                                       prefer="best")
    assert by_quality.keep[0].path.name == "g.ab1"


def test_runs_that_cannot_be_ordered_go_to_a_human(tmp_path):
    same = datetime(2026, 2, 1, 9, 0, 0)
    a = abif.write(tmp_path / "a.ab1", "ACGT" * 100, [40] * 400,
                   sample="CLN-011", well="A01", plate_id="P", run_started=same)
    b = abif.write(tmp_path / "b.ab1", "ACGT" * 100, [40] * 400,
                   sample="CLN-011_RR", well="A02", plate_id="P", run_started=same)
    report = sanger.resolve_reruns([abif.read(a), abif.read(b)])
    assert not report.keep and len(report.unresolved) == 2
    assert all(report.table["outcome"] == "unresolved")


# --- order forms -------------------------------------------------------------

@pytest.fixture
def plate() -> pd.DataFrame:
    return pd.DataFrame({
        "well": [wells.from_rc(r, c) for c in range(1, 13) for r in range(1, 9)][:95],
        "clone_id": [f"CLN-{i:03d}" for i in range(95)],
        "length_bp": 6800,
        "concentration_ng_ul": 60,
    })


def test_every_form_is_verified_against_a_real_template():
    from plateforge.emit import sequencing

    table = sequencing.available().set_index("form")
    assert set(table.index) >= {"genewiz", "elim", "ucberkeley", "generic"}
    assert (table["confidence"] == "verified").all()
    assert bool(table.loc["genewiz", "default"])
    for name in ("genewiz", "elim", "ucberkeley"):
        assert (FIXTURES.parent / table.loc[name, "template"]).exists()


def test_the_genewiz_well_columns_match_the_vendors_own_pairing(plate):
    """Azenta's form gives `Well (H)` and `Well (V)` as two orderings of the
    same 96 positions, paired row by row. Getting them the wrong way round
    transposes the entire plate, and nothing downstream could detect it."""
    openpyxl = pytest.importorskip("openpyxl")
    from plateforge.emit import sequencing

    form = sequencing.order(plate, "genewiz", primer="CMV_F")
    sheet = openpyxl.load_workbook(FORMS / "azenta_sanger_form_v2.xlsx",
                                   data_only=True)["Azenta Sanger Sequencing"]
    theirs = [(sheet.cell(r, 2).value, sheet.cell(r, 3).value)
              for r in range(2, 2 + len(plate))]
    ours = list(zip(form.frame["Well (H)"], form.frame["Well (V)"]))
    assert ours == theirs


def test_the_genewiz_headers_are_the_templates_headers_exactly(plate):
    openpyxl = pytest.importorskip("openpyxl")
    from plateforge.emit import sequencing

    sheet = openpyxl.load_workbook(FORMS / "azenta_sanger_form_v2.xlsx",
                                   data_only=True)["Azenta Sanger Sequencing"]
    theirs = [sheet.cell(1, c).value for c in range(2, 12)]
    form = sequencing.order(plate, "genewiz")
    assert list(form.frame.columns) == theirs


def test_the_elim_form_is_row_major_and_uses_one_digit_wells(plate):
    from plateforge.emit import sequencing

    form = sequencing.order(plate, "elim", primer="CMV_F")
    assert list(form.frame["Well"][:13]) == [f"A{c}" for c in range(1, 13)] + ["B1"]
    assert form.frame["Template Type"].eq("Plasmid").all()


def test_the_elim_form_rejects_a_name_it_would_reject_at_upload(plate):
    from plateforge.emit import sequencing

    bad = plate.copy()
    bad.loc[0, "clone_id"] = "CLN 000/α"
    with pytest.raises(ValueError, match="characters outside"):
        sequencing.order(bad, "elim")


def test_the_berkeley_form_is_column_major_in_two_blocks(plate):
    from plateforge.emit import sequencing

    form = sequencing.order(plate, "ucberkeley")
    assert list(form.frame["well #"][:9]) == ["A1", "B1", "C1", "D1", "E1",
                                              "F1", "G1", "H1", "A2"]
    assert list(form.frame["well # (7-12)"][:3]) == ["A7", "B7", "C7"]


def test_the_berkeley_form_refuses_a_completely_full_plate(plate):
    """Their template says to leave a well empty; it is how they confirm the
    plate's orientation, so a full plate is an error, not a preference."""
    from plateforge.emit import sequencing

    full = pd.DataFrame({
        "well": [wells.from_rc(r, c) for c in range(1, 13) for r in range(1, 9)],
        "clone_id": [f"CLN-{i:03d}" for i in range(96)]})
    with pytest.raises(ValueError, match="at least one well"):
        sequencing.order(full, "ucberkeley")
    assert sequencing.order(full, "ucberkeley", allow_issues=True).issues


def test_wells_are_normalised_on_the_way_in_whatever_the_caller_used():
    from plateforge.emit import sequencing

    ragged = pd.DataFrame({"well": ["a1", "B2", "c03"],
                           "clone_id": ["X", "Y", "Z"]})
    form = sequencing.order(ragged, "generic")
    assert list(form.frame["well"]) == ["A01", "B02", "C03"]


def test_a_form_writes_a_file_and_its_notes(tmp_path, plate):
    pytest.importorskip("openpyxl")
    from plateforge.emit import sequencing

    form = sequencing.order(plate, "genewiz", primer="CMV_F")
    path = form.write(tmp_path)
    assert path.exists() and path.suffix == ".xlsx"
    assert path.with_suffix(".NOTES.txt").exists()
    assert "column-major" in path.with_suffix(".NOTES.txt").read_text()
