"""Verifying clones from their Sanger traces: what is in each well, really.

The question a lab member asks after sequencing is not "what is the sequence".
It is **"which wells can I use"**, and the honest answers are more varied than
pass/fail:

    exact         the insert is what we ordered
    silent        bases differ, the protein does not -- usable
    missense      the protein differs -- usable only if you meant it to
    frameshift    an indel that is not a multiple of three -- dead
    indel         an in-frame insertion or deletion -- alive, not what we asked
    no_insert     the read is vector all the way through -- empty backbone
    wrong_clone   the read matches a different well's insert -- a plate swap
    mixed         two colonies in one well; the trace has double peaks
    low_quality   the trace never cleared the quality floor
    unreadable    no usable alignment at all

`wrong_clone` is the one worth building deliberately for. It is cheap to
detect -- align each read against every expected insert, not just its own --
and it is the failure that silently poisons everything downstream, because
every later step will faithfully carry the wrong antibody under the right
name.

## Screening on the insert, not the contig

Every well's expected sequence shares the same vector flanks, so a read is
*supposed* to match all 96 references across the leader and the constant
region. Screening on the whole contig therefore ranks everything equally and
finds no swaps at all. Measured on the real fixture traces, an M13F read from
an unrelated construct still shares 302 12-mers with pcDNA3.1, purely through
shared backbone. So the screen runs on the **insert** alone, where the 96
references genuinely differ.

## Reruns

See `resolve_reruns`. Short version: group by the sample name the submitter
typed, order by the run timestamp the instrument wrote, and never delete --
supersede (rule 10), and surface what could not be resolved instead of
guessing.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime

import pandas as pd

from ..core import wells as wellmod
from . import abif, codon, dna

EXACT = "exact"
SILENT = "silent"
MISSENSE = "missense"
FRAMESHIFT = "frameshift"
INDEL = "indel"
NO_INSERT = "no_insert"
WRONG_CLONE = "wrong_clone"
MIXED = "mixed"
LOW_QUALITY = "low_quality"
UNREADABLE = "unreadable"

USABLE = {EXACT, SILENT}
ORDER = [EXACT, SILENT, MISSENSE, INDEL, FRAMESHIFT, MIXED, NO_INSERT,
         WRONG_CLONE, LOW_QUALITY, UNREADABLE]

DEFAULTS = {
    "trim_cutoff": 20,          # Phred, for Mott trimming
    "min_trimmed_bp": 60,       # below this the read carries no information
    "min_insert_coverage": 0.80,  # of the expected insert, to call it at all
    "min_identity": 0.90,       # below this the read is not this sequence
    "swap_margin": 1.25,        # a rival must beat the expected by this much
    "min_screen_kmers": 20,     # fewer shared k-mers than this is not a match
}

# Tokens vendors append when they re-run a sample. Deliberately conservative,
# and deliberately not `_2`-style bare numbers: `PLATE_2` is a plate, not a
# second attempt, and wrongly grouping two distinct samples silently drops one
# of them as superseded. Over-grouping is the expensive mistake here;
# under-grouping just means a human is asked. Note that `\b` does not help --
# `_` is a word character, so there is no boundary between `P3` and `_RR`,
# which is how the first version of this quietly matched nothing.
RERUN_TOKENS = re.compile(
    r"(?:[-_ .]+(?:rr|rerun|re-?run|repeat|redo|run[2-9]|v[2-9]|r[2-9])"
    r"|\s*\(\d+\))\s*$", re.IGNORECASE)


def stem_of(name: str) -> tuple[str, bool]:
    """A sample name with any rerun marker stripped, and whether one was there."""
    cleaned = str(name).strip()
    marked = False
    while True:
        shorter = RERUN_TOKENS.sub("", cleaned).strip(" -_")
        if shorter == cleaned or not shorter:
            break
        cleaned, marked = shorter, True
    return cleaned, marked


@dataclass
class Expectation:
    """What a well is supposed to contain.

    `insert` is the ordered CDS -- the part that differs between wells, and
    the only part that can identify a clone. `context` is the whole expected
    contig, insert in its vector flanks, which is what a read is aligned
    against so that the flanks are matched rather than treated as mismatches.
    """
    clone_id: str
    well: str
    insert: str
    context: str = ""
    protein: str = ""

    def __post_init__(self):
        self.well = wellmod.normalize(self.well)
        self.insert = self.insert.upper()
        self.context = (self.context or self.insert).upper()
        if not self.protein:
            self.protein = codon.translate(self.insert).rstrip("*")

    @property
    def insert_at(self) -> int:
        """Where the insert starts inside the context."""
        return max(self.context.find(self.insert), 0)


@dataclass
class Call:
    """One trace, judged."""
    file: str
    verdict: str
    clone_id: str = ""
    well: str = ""
    sample: str = ""
    identity: float = 0.0
    coverage: float = 0.0
    confidence: float = 0.0
    trimmed_bp: int = 0
    mean_quality: float = 0.0
    matched: str = ""              # which expectation the read actually fits
    runner_up: str = ""
    differences: list = field(default_factory=list)
    protein_changes: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    run_started: datetime | None = None

    @property
    def usable(self) -> bool:
        return self.verdict in USABLE

    def row(self) -> dict:
        return {"file": self.file, "well": self.well, "clone_id": self.clone_id,
                "sample": self.sample, "verdict": self.verdict,
                "identity": round(self.identity, 4),
                "coverage": round(self.coverage, 4),
                "confidence": round(self.confidence, 3),
                "trimmed_bp": self.trimmed_bp,
                "mean_quality": round(self.mean_quality, 1),
                "matched": self.matched, "runner_up": self.runner_up,
                "n_differences": len(self.differences),
                "protein_changes": ";".join(self.protein_changes),
                "notes": "; ".join(self.notes),
                "run_started": self.run_started}


# --- mixed traces ------------------------------------------------------------

def mixed_fraction(trace: abif.Trace, start: int, end: int) -> float:
    """How much of the good window the basecaller called ambiguous.

    Two colonies in one well give a trace with two bases at every position
    where the clones differ, and the basecaller emits IUPAC ambiguity codes
    (R, Y, S, W, K, M) rather than N. A handful is noise; a few percent spread
    through an otherwise high-quality read is two templates.
    """
    window = trace.sequence[start:end].upper()
    if not window:
        return 0.0
    return sum(1 for base in window if base in "RYSWKMBDHV") / len(window)


# --- calling one read --------------------------------------------------------

def _protein_changes(expected: str, observed: str) -> list[str]:
    changes = []
    for index, (before, after) in enumerate(zip(expected, observed), start=1):
        if before != after:
            changes.append(f"{before}{index}{after}")
    if len(observed) != len(expected):
        changes.append(f"length {len(expected)}->{len(observed)}")
    return changes


def call_one(trace: abif.Trace, expectation: Expectation | None,
             others: dict[str, str] | None = None, **params) -> Call:
    """Judge one trace against what its well was supposed to contain.

    `others` maps clone_id to insert for every *other* well, so a read that
    belongs somewhere else can be named rather than merely rejected.
    """
    p = DEFAULTS | params
    call = Call(file=trace.path.name, verdict=UNREADABLE, sample=trace.sample,
                well=wellmod.normalize(trace.well) if trace.well else "",
                mean_quality=trace.mean_quality, run_started=trace.run_started)
    if expectation is not None:
        call.clone_id, call.well = expectation.clone_id, expectation.well

    start, end = abif.mott_trim(trace.quality, p["trim_cutoff"])
    call.trimmed_bp = end - start
    if call.trimmed_bp < p["min_trimmed_bp"]:
        call.verdict = LOW_QUALITY
        call.notes.append(
            f"only {call.trimmed_bp} bp survived quality trimming at Q"
            f"{p['trim_cutoff']}; mean quality {trace.mean_quality:.0f}. "
            "Re-run this well rather than reading anything into it.")
        return call

    read = trace.sequence[start:end]
    ambiguous = mixed_fraction(trace, start, end)

    if expectation is None:
        call.notes.append("no expected sequence for this well; read not judged")
        return call

    # Which insert does this read actually carry? Screened on inserts alone,
    # because the vector flanks are identical across every well.
    library = dict(others or {})
    library[expectation.clone_id] = expectation.insert
    ranked = dna.screen(read, library, top=3)
    best = ranked[0] if ranked else None
    second = ranked[1] if len(ranked) > 1 else None
    mine = next((c for c in ranked if c.name == expectation.clone_id), None)

    if best is not None:
        call.matched = best.name
        call.runner_up = second.name if second else ""
        rival = second.score if second else 0
        call.confidence = round(1.0 - (rival / best.score), 3) if best.score else 0.0

    if best is None or best.score < p["min_screen_kmers"]:
        # Nothing in the panel matched. Whatever the screen's best row was, it
        # is noise, and leaving it in `matched` makes an empty well look like
        # it holds a named clone.
        call.matched = call.runner_up = ""
        call.confidence = 0.0
        # Either the well is empty backbone or the read is junk; the alignment
        # against the context tells which.
        against_context = dna.align_to(expectation.context, read)
        if against_context and against_context.identity >= p["min_identity"]:
            call.verdict = NO_INSERT
            call.identity = against_context.identity
            call.notes.append(
                "the read matches the vector but carries none of the expected "
                "insert: empty backbone, or the insert dropped out.")
        else:
            call.verdict = UNREADABLE
            call.notes.append(
                "the read matches neither the expected insert nor the vector; "
                "wrong template, wrong primer, or a failed reaction.")
        return call

    if (mine is None or best.name != expectation.clone_id) and \
            best.score >= (mine.score if mine else 0) * p["swap_margin"]:
        call.verdict = WRONG_CLONE
        # How well it matches the clone it actually is, which is the number
        # that says whether this is a swap or merely a bad read.
        rival_insert = library.get(best.name, "")
        against_rival = dna.align_to(
            rival_insert, read if best.forward else dna.reverse_complement(read))
        if against_rival is not None:
            call.identity = against_rival.identity
        call.notes.append(
            f"this read matches {best.name} ({best.score} shared 12-mers"
            + (f", {call.identity:.1%} identity" if call.identity else "")
            + f") better than {expectation.clone_id} "
            f"({mine.score if mine else 0}). Check the pick list before "
            "trusting any later result for this well.")
        return call

    # It is the right clone. Now: is it the right *sequence*?
    if not best.forward:
        read = dna.reverse_complement(read)
        call.notes.append("read is on the reverse strand; reverse-complemented")

    aligned = dna.align_to(expectation.context, read)
    if aligned is None:
        call.verdict = UNREADABLE
        call.notes.append("no alignment against the expected contig")
        return call

    call.identity = aligned.identity
    insert_start = expectation.insert_at
    insert_end = insert_start + len(expectation.insert)
    covered = max(0, min(aligned.ref_end, insert_end) - max(aligned.ref_start, insert_start))
    call.coverage = covered / len(expectation.insert) if expectation.insert else 0.0

    # Only differences that land inside the insert are this clone's problem.
    # A mismatch out in the constant region is the vector's, and reporting it
    # per well is how a real defect gets lost in 96 copies of the same noise.
    found = [d for d in dna.differences(aligned)
             if insert_start <= d["ref_start"] < insert_end]
    call.differences = found

    if call.coverage < p["min_insert_coverage"]:
        call.verdict = LOW_QUALITY
        call.notes.append(
            f"the read covers only {call.coverage:.0%} of the insert; too "
            "little to call. A longer read, or sequence from the other end.")
        return call

    if ambiguous >= 0.02:
        call.verdict = MIXED
        call.notes.append(
            f"{ambiguous:.1%} of the good window is an ambiguity code, not a "
            "single base: two templates in one well. Re-streak and re-pick.")
        return call

    if call.identity < p["min_identity"]:
        call.verdict = UNREADABLE
        call.notes.append(f"identity {call.identity:.1%} against the expected "
                          "contig is too low to interpret")
        return call

    indels = [d for d in found if d["kind"] in ("insertion", "deletion")]
    if indels:
        shift = sum(d["length"] * (1 if d["kind"] == "insertion" else -1)
                    for d in indels)
        call.verdict = FRAMESHIFT if shift % 3 else INDEL
        call.notes.append(
            f"{len(indels)} indel(s), net {shift:+d} bp"
            + ("; the reading frame is broken from the first one onward."
               if shift % 3 else "; in frame, but not the sequence ordered."))
        return call

    if not found:
        call.verdict = EXACT
        return call

    observed = list(expectation.context)
    for difference in found:
        observed[difference["ref_start"]] = difference["bases"]
    observed_insert = "".join(observed)[insert_start:insert_end]
    protein = codon.translate(observed_insert).rstrip("*")
    if protein == expectation.protein:
        call.verdict = SILENT
        call.notes.append(f"{len(found)} silent base change(s); protein is "
                          "the sequence ordered.")
    else:
        call.verdict = MISSENSE
        call.protein_changes = _protein_changes(expectation.protein, protein)
        call.notes.append(", ".join(call.protein_changes[:6]))
    return call


# --- reruns ------------------------------------------------------------------

@dataclass
class RerunReport:
    """Which traces to use, which were replaced, and which we could not tell."""
    keep: list
    superseded: list
    unresolved: list
    table: pd.DataFrame


def resolve_reruns(traces: list[abif.Trace], prefer: str = "later") -> RerunReport:
    """Work out which of several traces of the same sample to believe.

    The grouping key is the **sample name the submitter typed**, with any
    recognised rerun marker stripped -- not the filename, and not the plate
    and well. A rerun is usually cherry-picked into a fresh plate, so it
    arrives with a different barcode and a different well; the sample name is
    the only thing that follows it.

    Ordering is by the run timestamp the sequencer wrote into the trace, never
    by the file's modification time. Downloading rewrites mtime, and the rerun
    is often downloaded first, so mtime ordering is not merely unreliable, it
    is backwards in exactly the case that matters.

    `prefer="later"` takes the most recent run. `prefer="best"` takes the one
    with the most Q20 bases, which is the more honest rule -- **a rerun can
    fail too, and nothing here assumes the second attempt was better** -- but
    it is not the default, because a rerun is usually done for a reason the
    trace cannot see.

    Nothing is deleted. Superseded traces are returned so the caller can
    record the replacement (rule 10), and anything ambiguous lands in
    `unresolved` for a human rather than being resolved by guesswork.
    """
    groups: dict[str, list] = {}
    labels: dict[str, set] = {}
    for trace in traces:
        label = trace.sample or trace.path.stem
        stem, _marked = stem_of(label)
        key = stem or label
        groups.setdefault(key, []).append(trace)
        labels.setdefault(key, set()).add(label)

    keep, superseded, unresolved, rows = [], [], [], []
    for stem, members in sorted(groups.items()):
        if len(members) == 1:
            keep.append(members[0])
            rows.append({"stem": stem, "file": members[0].path.name,
                         "outcome": "only run", "why": ""})
            continue

        timed = [t for t in members if t.run_started is not None]
        if len(timed) < len(members):
            unresolved.extend(members)
            for trace in members:
                rows.append({"stem": stem, "file": trace.path.name,
                             "outcome": "unresolved",
                             "why": "at least one trace has no run timestamp, "
                                    "so the attempts cannot be ordered"})
            continue

        if prefer == "best":
            ranked = sorted(members, key=lambda t: t.bases_at_least(20), reverse=True)
            reason = "most Q20 bases"
        else:
            ranked = sorted(members, key=lambda t: t.run_started, reverse=True)
            reason = "latest run"

        # Two runs at the same instant, or the later run being clearly worse,
        # are both cases where the rule stops being obviously right.
        winner, rest = ranked[0], ranked[1:]
        tied = [t for t in rest if t.run_started == winner.run_started]
        if tied:
            unresolved.extend(members)
            for trace in members:
                rows.append({"stem": stem, "file": trace.path.name,
                             "outcome": "unresolved",
                             "why": "two runs share a timestamp"})
            continue

        keep.append(winner)
        rows.append({"stem": stem, "file": winner.path.name, "outcome": "kept",
                     "why": f"{reason}; grouped from sample names "
                            + ", ".join(sorted(labels.get(stem, set())))})
        for trace in rest:
            superseded.append((trace, winner))
            worse = winner.bases_at_least(20) < trace.bases_at_least(20)
            rows.append({
                "stem": stem, "file": trace.path.name, "outcome": "superseded",
                "why": f"replaced by {winner.path.name}" + (
                    f" -- WARNING: the replacement has fewer Q20 bases "
                    f"({winner.bases_at_least(20)} against "
                    f"{trace.bases_at_least(20)}); check before relying on it"
                    if worse else "")})

    return RerunReport(keep=keep, superseded=superseded, unresolved=unresolved,
                       table=pd.DataFrame(rows))


# --- calling a plate ---------------------------------------------------------

def expectations_from(frame: pd.DataFrame, vector=None, *,
                      clone_column: str = "clone_id",
                      well_column: str = "well",
                      dna_column: str = "cds_dna") -> dict[str, Expectation]:
    """Build one expectation per well from a clone table.

    `vector` is anything with `orf(insert)` -- a `library.vector.Vector` --
    used to wrap each insert in its flanks. Without one the insert is its own
    context, which still works and simply cannot recognise empty backbone.
    """
    out = {}
    for _, row in frame.iterrows():
        insert = str(row[dna_column]).upper()
        context = vector.orf(insert) if vector is not None else insert
        expectation = Expectation(clone_id=str(row[clone_column]),
                                  well=str(row[well_column]),
                                  insert=insert, context=context)
        out[expectation.well] = expectation
    return out


def match_traces(traces: list[abif.Trace],
                 expectations: dict[str, Expectation],
                 by: str = "well") -> list[tuple]:
    """Pair each trace with the expectation for its well.

    `by="well"` uses the well the instrument recorded in the trace, which is
    the only identifier that did not pass through a human. `by="sample"`
    matches the sample name against the clone id instead, for submissions
    where the plate was rearranged between ordering and sequencing.
    """
    pairs = []
    for trace in traces:
        if by == "sample":
            stem, _ = stem_of(trace.sample)
            found = next((e for e in expectations.values()
                          if e.clone_id == stem or stem.startswith(e.clone_id)), None)
        else:
            well = wellmod.normalize(trace.well) if trace.well else ""
            found = expectations.get(well)
        pairs.append((trace, found))
    return pairs


def verify(traces: list[abif.Trace], expectations: dict[str, Expectation],
           *, by: str = "well", resolve: bool = True,
           prefer: str = "later", **params) -> tuple[pd.DataFrame, RerunReport | None]:
    """Judge a whole plate of traces. Returns the calls and the rerun report."""
    report = resolve_reruns(traces, prefer=prefer) if resolve else None
    usable = report.keep + report.unresolved if report else list(traces)

    inserts = {e.clone_id: e.insert for e in expectations.values()}
    calls = []
    for trace, expectation in match_traces(usable, expectations, by=by):
        others = {k: v for k, v in inserts.items()
                  if expectation is None or k != expectation.clone_id}
        calls.append(call_one(trace, expectation, others, **params).row())

    frame = pd.DataFrame(calls)
    if not frame.empty:
        frame["verdict"] = pd.Categorical(frame["verdict"], categories=ORDER,
                                          ordered=True)
        frame = frame.sort_values(["verdict", "well"]).reset_index(drop=True)
    return frame, report


def summarise(calls: pd.DataFrame) -> pd.DataFrame:
    """How many wells landed in each verdict, worst-news-first."""
    if calls.empty:
        return pd.DataFrame(columns=["verdict", "wells", "usable"])
    counts = (calls.groupby("verdict", observed=True).size()
              .rename("wells").reset_index())
    counts["usable"] = counts["verdict"].isin(USABLE)
    return counts


def usable_wells(calls: pd.DataFrame) -> list[str]:
    """The wells worth carrying forward, in plate order."""
    if calls.empty:
        return []
    good = calls[calls["verdict"].isin(USABLE)]
    return sorted(good["well"].dropna().unique())
