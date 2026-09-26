"""DNA alignment, for comparing a sequencing read against what we expect.

`library.msa` aligns protein against a germline with BLOSUM62 and cares about
IMGT regions. None of that applies here: this is DNA against DNA, the two are
supposed to be identical, and the interesting output is *where they are not*.

Two pieces, and the second exists because of the first:

**A k-mer screen.** The expensive question is not "does this read match its
own well" but "does it match some *other* well better", because that is a
plate swap and it is the error that quietly poisons a campaign. Answering it
means comparing one read against all 96 expected inserts. Full dynamic
programming 96 times per read is minutes; counting shared 12-mers is
milliseconds, and it is a reliable screen because a real match shares hundreds
of exact 12-mers and an unrelated antibody V region shares a handful.

**A banded alignment.** Once the screen has picked candidates, the shared
k-mers also say *where* on the reference the read sits. Aligning inside a band
around that diagonal is a few hundred thousand cells instead of a few million,
and the band is derived from the data rather than assumed -- if the k-mer
offsets disagree, the band widens to cover them, which is exactly the case
where an indel is present.

End gaps are free at both ends of both sequences (an overlap alignment). The
read starts somewhere upstream in the promoter and runs out partway through
the constant region; neither overhang is a mutation and neither should be
scored as one.
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                           "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

MATCH, MISMATCH = 2.0, -4.0
# Opening a gap is deliberately expensive relative to extending one. With a
# cheap open, a single 6 bp deletion is reported as three small deletions
# scattered over eight bases, because the pieces pick up enough coincidental
# matches between them to pay for the extra opens. That is arithmetically
# optimal and biologically nonsense: it is one event, and an in-frame one.
GAP_OPEN, GAP_EXTEND = -14.0, -1.0
K = 12


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# --- the screen --------------------------------------------------------------

def kmers(seq: str, k: int = K) -> dict[str, list[int]]:
    """Every k-mer and where it starts. Ambiguous bases are skipped.

    A basecaller emits N generously at the start of a trace, and a k-mer
    containing one matches nothing, so indexing them wastes work and dilutes
    the score.
    """
    index: dict[str, list[int]] = {}
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        word = seq[i:i + k]
        if "N" in word:
            continue
        index.setdefault(word, []).append(i)
    return index


def shared(query_index: dict, reference: str, k: int = K) -> int:
    """How many of the query's k-mers occur in the reference."""
    reference = reference.upper()
    present = {reference[i:i + k] for i in range(len(reference) - k + 1)}
    return sum(1 for word in query_index if word in present)


@dataclass
class Candidate:
    """One reference the read might belong to, and how strongly."""
    name: str
    score: int                 # shared k-mers
    forward: bool
    fraction: float            # of the read's k-mers, so reads compare fairly


def screen(read: str, references: dict[str, str], k: int = K,
           top: int = 3) -> list[Candidate]:
    """Rank references by shared k-mers, both strands, best first.

    Returns at most `top`, but always enough to see a runner-up: the gap
    between first and second is what says whether a match is unambiguous, and
    a single result cannot show it.
    """
    forward_index = kmers(read, k)
    reverse_index = kmers(reverse_complement(read), k)
    total = max(len(forward_index), 1)

    found = []
    for name, reference in references.items():
        ahead = shared(forward_index, reference, k)
        behind = shared(reverse_index, reference, k)
        if behind > ahead:
            found.append(Candidate(name, behind, False, behind / max(len(reverse_index), 1)))
        else:
            found.append(Candidate(name, ahead, True, ahead / total))
    found.sort(key=lambda c: c.score, reverse=True)
    return found[:max(top, 2)]


def diagonal(read: str, reference: str, k: int = K) -> tuple[int, int] | None:
    """Where the read sits on the reference, from the shared k-mers.

    Returns `(offset, spread)` where reference position is about
    `query position + offset`. `spread` is how much the individual k-mer
    offsets disagree, which is the band the alignment needs -- an indel of 6 bp
    shows up here as a spread of 6 before any alignment has been run.
    """
    index = kmers(reference, k)
    offsets = []
    read = read.upper()
    for i in range(len(read) - k + 1):
        word = read[i:i + k]
        hits = index.get(word)
        if hits and len(hits) == 1:      # unique anchors only; repeats mislead
            offsets.append(hits[0] - i)
    if not offsets:
        return None
    middle = int(statistics.median(offsets))
    spread = max(abs(o - middle) for o in offsets)
    return middle, spread


# --- the alignment -----------------------------------------------------------

@dataclass
class Alignment:
    """Two gapped strings, and where on the reference they start and stop."""
    reference: str              # gapped
    query: str                  # gapped
    score: float
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int

    @property
    def length(self) -> int:
        return len(self.reference)

    @property
    def matches(self) -> int:
        return sum(1 for a, b in zip(self.reference, self.query) if a == b and a != "-")

    @property
    def mismatches(self) -> int:
        return sum(1 for a, b in zip(self.reference, self.query)
                   if a != b and a != "-" and b != "-")

    @property
    def gaps(self) -> int:
        return sum(1 for a, b in zip(self.reference, self.query) if a == "-" or b == "-")

    @property
    def identity(self) -> float:
        covered = self.matches + self.mismatches + self.gaps
        return self.matches / covered if covered else 0.0


def align(reference: str, query: str, *, band: int | None = None,
          offset: int | None = None, match: float = MATCH,
          mismatch: float = MISMATCH, gap_open: float = GAP_OPEN,
          gap_extend: float = GAP_EXTEND) -> Alignment | None:
    """Align `query` to `reference`, inside a band when one is given.

    Local at both ends: the alignment may begin and end anywhere, so the
    junk a basecaller puts at the start of a trace, and the point where the
    read runs out partway through the reference, are both free. Scoring them
    as gaps is what turns 40 bases of noise into a page of invented indels.

    Affine gaps, with an expensive open: one 6 bp deletion should cost far
    less than six scattered single-base deletions.

    Returns None when nothing aligns at all, rather than a best-effort
    alignment of two unrelated sequences.
    """
    reference, query = reference.upper(), query.upper()
    n, m = len(query), len(reference)
    if not n or not m:
        return None

    if band is None or offset is None:
        low_of = lambda i: 0                                  # noqa: E731
        high_of = lambda i: m                                 # noqa: E731
    else:
        def low_of(i, _o=offset, _b=band):
            return max(0, i + _o - _b)

        def high_of(i, _o=offset, _b=band):
            return min(m, i + _o + _b + 1)

    negative = float("-inf")
    # best[j]: the best score reaching (i, j) in any state
    # up[j]:   best score reaching (i, j) with a gap in the reference
    # left:    best score reaching (i, j) with a gap in the query
    best = [0.0] * (m + 1)
    up = [negative] * (m + 1)
    trace: list[bytes] = []
    end_i = end_j = 0
    top = 0.0

    # One byte per cell, packed: bits 0-1 say which state won (DIAG, UP, LEFT,
    # or STOP for a local restart); bit 2 is set when the UP gap was opened
    # here rather than extended, bit 3 the same for LEFT. Both bits are needed.
    # Traceback through affine gaps is a state machine: once inside a gap you
    # must keep following that gap's chain until it opens. Re-reading each
    # cell's winning state instead lets a gap close and immediately reopen,
    # which costs nothing in the score -- the DP is still right -- but turns a
    # single 6 bp deletion into two smaller ones in the reported alignment.
    DIAG, UP, LEFT, STOP = 0, 1, 2, 3
    UP_OPENED, LEFT_OPENED = 4, 8

    for i in range(1, n + 1):
        low, high = low_of(i - 1), high_of(i - 1)
        new_best = [0.0] * (m + 1)
        new_up = [negative] * (m + 1)
        left = negative
        left_opened = True
        row = bytearray(m + 1)
        for j in range(max(low, 1), high + 1):
            if j > m:
                break
            diagonal_score = best[j - 1] + (match if query[i - 1] == reference[j - 1]
                                            else mismatch)
            open_up, extend_up = best[j] + gap_open, up[j] + gap_extend
            new_up[j] = max(open_up, extend_up)
            open_left, extend_left = new_best[j - 1] + gap_open, left + gap_extend
            left = max(open_left, extend_left)
            left_opened = open_left >= extend_left

            take = max(diagonal_score, new_up[j], left)
            if take <= 0:
                # Starting fresh here beats carrying a negative score in, so
                # everything before this point is a free overhang.
                new_best[j], row[j] = 0.0, STOP
                continue
            new_best[j] = take
            state = (DIAG if take == diagonal_score
                     else (UP if take == new_up[j] else LEFT))
            flags = (UP_OPENED if open_up >= extend_up else 0) | \
                    (LEFT_OPENED if left_opened else 0)
            row[j] = state | flags
            if take > top:
                top, end_i, end_j = take, i, j
        trace.append(bytes(row))
        best, up = new_best, new_up

    if top <= 0:
        return None

    ref_out, query_out = [], []
    i, j = end_i, end_j
    state = DIAG
    while i > 0 and j > 0:
        cell = trace[i - 1][j]
        move = state if state != DIAG else (cell & 3)
        if move == STOP:
            break
        if move == DIAG:
            ref_out.append(reference[j - 1])
            query_out.append(query[i - 1])
            i, j = i - 1, j - 1
            state = DIAG
        elif move == UP:                      # gap in the reference
            ref_out.append("-")
            query_out.append(query[i - 1])
            state = DIAG if cell & UP_OPENED else UP
            i -= 1
        else:                                 # gap in the query
            ref_out.append(reference[j - 1])
            query_out.append("-")
            state = DIAG if cell & LEFT_OPENED else LEFT
            j -= 1

    return Alignment(reference="".join(reversed(ref_out)),
                     query="".join(reversed(query_out)),
                     score=top, ref_start=j, ref_end=end_j,
                     query_start=i, query_end=end_i)


def align_to(reference: str, query: str, *, k: int = K,
             minimum_band: int = 24, **kwargs) -> Alignment | None:
    """Align, choosing the band from the data.

    Falls back to a full alignment when the k-mer screen finds no anchor at
    all, because "no shared 12-mer" is itself informative -- it means the read
    is not this sequence -- and the full alignment is what proves it.
    """
    placed = diagonal(query, reference, k)
    if placed is None:
        return align(reference, query, **kwargs)
    offset, spread = placed
    band = max(minimum_band, spread + minimum_band)
    return align(reference, query, band=band, offset=offset, **kwargs)


def differences(alignment: Alignment) -> list[dict]:
    """Every place the read disagrees with the reference.

    Runs of gaps are reported as one event, not one per base: a 3 bp deletion
    is a single in-frame deletion and calling it three deletions makes an
    honest read look catastrophic.
    """
    out: list[dict] = []
    ref_at = alignment.ref_start
    run: dict | None = None

    for ref_base, query_base in zip(alignment.reference, alignment.query):
        if ref_base == "-" or query_base == "-":
            kind = "insertion" if ref_base == "-" else "deletion"
            if run and run["kind"] == kind and run["ref_end"] == ref_at:
                run["length"] += 1
                run["bases"] += query_base if kind == "insertion" else ref_base
                if kind == "deletion":
                    run["ref_end"] = ref_at + 1
            else:
                if run:
                    out.append(run)
                run = {"kind": kind, "ref_start": ref_at,
                       "ref_end": ref_at + (0 if kind == "insertion" else 1),
                       "length": 1,
                       "bases": query_base if kind == "insertion" else ref_base}
        else:
            if run:
                out.append(run)
                run = None
            if ref_base != query_base:
                out.append({"kind": "substitution", "ref_start": ref_at,
                            "ref_end": ref_at + 1, "length": 1,
                            "bases": query_base, "was": ref_base})
        if ref_base != "-":
            ref_at += 1
    if run:
        out.append(run)
    return out
