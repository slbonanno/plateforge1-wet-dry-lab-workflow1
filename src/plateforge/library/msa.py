"""Multiple sequence alignment against a germline reference.

What this is for: given a set of VH sequences that all call the same V gene,
lay them out in shared columns so a figure can show where each one departs
from its germline.

Why not IMGT numbering. Numbering is a *classification* -- it tells you which
structural position a residue occupies. It is the right tool for comparing
across genes, and the wrong tool for showing an alignment of one family
against its own germline, because it needs every sequence to carry ANARCI
output and it silently drops any sequence that does not. An alignment does not
need numbering. It needs gaps.

How it works. Every sequence is aligned to the germline pairwise
(Needleman-Wunsch, BLOSUM62, affine gaps), and the pairwise results are then
merged into one column space anchored on the germline:

    germline   E V Q L L E S G G G . . . L V Q P
    seq_1      E V Q L V E S G G G . . . L V Q P
    seq_2      E V Q L - E S G G G L Y R L V Q P
                                   ^^^^^ insertion columns

Germline residues own their columns, so germline column 5 is germline residue
5 in every row. Where a sequence carries residues the germline does not, an
insertion column opens and every other row shows a gap there. Insertion
columns are not germline positions and so are never numbered -- which is what
makes the numbering on the reference row honest.

This is a reference-anchored progressive alignment, not a full-blown
tree-guided MSA. For a set of sequences that are all one germline's children
the two agree, and this one has no external dependency, is deterministic, and
puts the germline in column space by construction. An external aligner
(mafft, clustalo) can be used instead via `backend=`; the germline goes in as
the first sequence and is pulled back out of the result the same way.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

GAP = "-"

# BLOSUM62, as a flat table. Antibody V regions against their own germline are
# high-identity, so the substitution matrix matters much less than the gap
# model -- but scoring a conservative change as cheaper than a radical one
# still keeps the aligner from opening a gap to dodge a mismatch.
_ORDER = "ARNDCQEGHILKMFPSTWYVBZX*"
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def _load_b62() -> dict[tuple[str, str], int]:
    rows = [r.split() for r in _B62_ROWS.strip().splitlines()]
    return {(a, b): int(rows[i][j])
            for i, a in enumerate(_ORDER) for j, b in enumerate(_ORDER)}


BLOSUM62 = _load_b62()


def score_pair(a: str, b: str) -> int:
    return BLOSUM62.get((a, b), BLOSUM62.get(("X", "X"), -1))


def align_pair(ref: str, seq: str, gap_open: float = -11, gap_extend: float = -1,
               end_gap_free: bool = True,
               gap_profile: list[tuple[float, float]] | None = None
               ) -> tuple[str, str]:
    """Global pairwise alignment, affine gaps (Gotoh).

    `end_gap_free` makes leading and trailing gaps in the *query* free. Real
    OAS reads are 5' truncated -- a sequence that starts at germline residue
    16 is not paying for 15 deletions it never had, and without this the
    aligner would rather shift the whole read than eat the opening penalty.

    `gap_profile` gives a (open, extend) pair per reference position and
    governs INSERTIONS only -- gaps in the reference, where the query carries
    residues the germline does not. Deletions (gaps in the query, a germline
    residue the read does not have) keep the flat penalty everywhere, because
    "this sequence is missing framework" is a strong claim at any position and
    must stay expensive. Relaxing both is how the junction ends up eating its
    flanking anchors.

    The cost of an insertion depends on WHERE it is. One penalty for a whole antibody is
    wrong in a specific way: framework regions are near-invariant and a gap in
    one is almost certainly an alignment error, while CDR3 is the product of
    V(D)J recombination plus N-additions and has no germline to be indel-free
    against. A uniform penalty makes the aligner pay framework prices to open
    the junction, so it buys mismatches instead and smears CDR3 across
    framework columns. See `region_gap_profile`.
    """
    n, m = len(ref), len(seq)
    if n == 0 or m == 0:
        return ref or GAP * m, seq or GAP * n

    def gaps(i: int) -> tuple[float, float]:
        """(open, extend) for a gap at reference position i."""
        if not gap_profile:
            return gap_open, gap_extend
        return gap_profile[min(max(i, 0), len(gap_profile) - 1)]

    NEG = float("-inf")
    # M: ends aligned. X: gap in seq (ref consumed). Y: gap in ref (seq consumed).
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    pM = [[0] * (m + 1) for _ in range(n + 1)]
    pX = [[0] * (m + 1) for _ in range(n + 1)]
    pY = [[0] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = 0 if end_gap_free else gap_open + gap_extend * (i - 1)
        pX[i][0] = 1
    o0, e0 = gaps(0)
    for j in range(1, m + 1):
        Y[0][j] = o0 + e0 * (j - 1)
        pY[0][j] = 2

    for i in range(1, n + 1):
        ri = ref[i - 1]
        Mi, Mp = M[i], M[i - 1]
        Xi, Xp = X[i], X[i - 1]
        Yi = Y[i]
        ins_open, ins_extend = gaps(i)
        for j in range(1, m + 1):
            s = score_pair(ri, seq[j - 1])
            best, src = Mp[j - 1], 0
            if X[i - 1][j - 1] > best:
                best, src = X[i - 1][j - 1], 1
            if Y[i - 1][j - 1] > best:
                best, src = Y[i - 1][j - 1], 2
            Mi[j] = best + s
            pM[i][j] = src

            # gap in seq: consume a ref residue
            free = end_gap_free and j == m
            open_x = Mp[j] + (0 if free else gap_open)
            ext_x = Xp[j] + (0 if free else gap_extend)
            if open_x >= ext_x:
                Xi[j], pX[i][j] = open_x, 0
            else:
                Xi[j], pX[i][j] = ext_x, 1

            # gap in ref: consume a seq residue (an insertion)
            open_y = Mi[j - 1] + ins_open
            ext_y = Yi[j - 1] + ins_extend
            if open_y >= ext_y:
                Yi[j], pY[i][j] = open_y, 0
            else:
                Yi[j], pY[i][j] = ext_y, 2

    i, j = n, m
    state = max(((M[i][j], 0), (X[i][j], 1), (Y[i][j], 2)))[1]
    out_ref: list[str] = []
    out_seq: list[str] = []
    while i > 0 or j > 0:
        if state == 0:
            out_ref.append(ref[i - 1])
            out_seq.append(seq[j - 1])
            state = pM[i][j]
            i, j = i - 1, j - 1
        elif state == 1:
            out_ref.append(ref[i - 1])
            out_seq.append(GAP)
            state = pX[i][j]
            i -= 1
        else:
            out_ref.append(GAP)
            out_seq.append(seq[j - 1])
            state = pY[i][j]
            j -= 1
        if i == 0 and j > 0:
            state = 2
        elif j == 0 and i > 0:
            state = 1
    return "".join(reversed(out_ref)), "".join(reversed(out_seq))


@dataclass
class Msa:
    """A finished alignment: the reference row, the rows below it, numbering."""
    ref_label: str
    ref: str                       # gapped reference, one char per column
    labels: list[str]
    rows: list[str]                # gapped, all exactly len(ref)
    numbers: list[int | None]      # per column: germline residue number, or None
    ref_source: str = "unknown"
    backend: str = "reference-anchored"
    ref_numbering: list = field(default_factory=list)
    regions: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def column_regions(self) -> list[str | None]:
        """The IMGT region of every column, insertion columns included.

        A column holding a reference residue takes that residue's region. An
        insertion column takes the region of the position that would come
        NEXT after the preceding reference residue -- so the block between the
        conserved cysteine (104) and tryptophan (118) is CDR3, which is the
        whole point, since the reference has no residues of its own in there.
        """
        numbers = self.ref_numbering or [None] * len(self.ref.replace("-", ""))
        out: list[str | None] = []
        at, previous = 0, None
        for aa in self.ref:
            if aa not in ".-":
                num = numbers[at] if at < len(numbers) else None
                out.append(region_of(num))
                previous = num
                at += 1
            else:
                nxt = numbers[at] if at < len(numbers) else None
                if previous is not None:
                    out.append(region_of(previous + 1))
                else:
                    out.append(region_of(nxt))
        return out

    def region_columns(self) -> list[tuple[str, int, int]]:
        """IMGT regions as (name, first column, last column + 1), tiling."""
        per_column = self.column_regions()
        if not any(per_column):
            return []
        spans, current, start = [], per_column[0], 0
        for c, name in enumerate(per_column[1:], start=1):
            if name != current:
                if current is not None:
                    spans.append((current, start, c))
                current, start = name, c
        if current is not None:
            spans.append((current, start, len(per_column)))
        return spans

    @property
    def width(self) -> int:
        return len(self.ref)

    @property
    def n_sequences(self) -> int:
        return len(self.rows)

    def identity(self, row: str) -> float:
        """Fraction of germline positions this row matches. Gaps count against."""
        pairs = [(r, q) for r, q in zip(self.ref, row) if r != GAP]
        return sum(r == q for r, q in pairs) / len(pairs) if pairs else 0.0

    def differences(self, row: str) -> list[tuple[int, str, str]]:
        """(germline number, germline residue, observed) for every mismatch."""
        out = []
        for num, r, q in zip(self.numbers, self.ref, row):
            if num is not None and r != q:
                out.append((num, r, q))
        return out

    def as_fasta(self) -> str:
        parts = [f">{self.ref_label}\n{self.ref}"]
        parts += [f">{lab}\n{row}" for lab, row in zip(self.labels, self.rows)]
        return "\n".join(parts) + "\n"

    def as_text(self, block: int = 60, max_label: int = 20) -> str:
        """Clustal-ish blocks, with germline numbering over the reference only."""
        labels = [self.ref_label] + list(self.labels)
        rows = [self.ref] + list(self.rows)
        pad = min(max(len(x) for x in labels), max_label)
        out: list[str] = []
        for start in range(0, self.width, block):
            stop = min(start + block, self.width)
            ticks = [" "] * (stop - start)
            for col in range(start, stop):
                num = self.numbers[col]
                if num is not None and (num % 10 == 0 or num == 1):
                    text = str(num)
                    at = col - start - len(text) + 1
                    if at >= 0:
                        ticks[at:at + len(text)] = list(text)
            out.append(" " * (pad + 1) + "".join(ticks[:stop - start]))
            for lab, row in zip(labels, rows):
                out.append(f"{lab[:pad]:<{pad}} {row[start:stop]}")
            out.append("")
        return "\n".join(out)


def build(ref: str, seqs: dict[str, str], ref_label: str = "germline",
          ref_source: str = "unknown", backend: str = "reference",
          region_aware: bool = True, **kwargs) -> Msa:
    """Align every sequence in `seqs` to `ref` and merge into one column space.

    `region_aware` charges framework prices for framework gaps and junction
    prices for junction gaps, using the reference's IMGT numbering. It is on
    by default and silently inert when the reference carries no numbering.
    """
    ref = "".join(c for c in ref.upper() if c not in ".-")
    clean = {lab: "".join(c for c in s.upper() if c not in ".-")
             for lab, s in seqs.items() if isinstance(s, str) and s.strip()}
    if not clean:
        raise ValueError("no sequences to align")

    fn = BACKENDS.get(backend)
    if fn is None:
        raise KeyError(f"unknown backend {backend!r}; known: {sorted(BACKENDS)}")

    numbers = numbering_for(ref)
    if region_aware and backend == "reference" and any(x is not None for x in numbers):
        kwargs.setdefault("gap_profile", region_gap_profile(numbers))
        kwargs.setdefault("insertion_layout", insertion_layout(numbers))
        kwargs.setdefault("anchors", _anchor_indices(numbers))

    msa = fn(ref, clean, **kwargs)
    msa.ref_label = ref_label
    msa.ref_source = ref_source
    msa.ref_numbering = numbers
    msa.regions = regions_of(ref)
    return msa


# --- backends --------------------------------------------------------------

BACKENDS: dict[str, object] = {}


def register_backend(name):
    def deco(fn):
        BACKENDS[name] = fn
        return fn
    return deco


def _decompose(gapped_ref: str, gapped_seq: str) -> tuple[list[str], list[str]]:
    """Split a pairwise result into per-germline-position residues + insertions.

    Returns (at, before) where at[i] is what sits on germline residue i and
    before[i] is the run of residues inserted immediately before it. before has
    one extra slot at the end for anything past the germline's last residue.
    """
    n_ref = sum(1 for c in gapped_ref if c != GAP)
    at = [GAP] * n_ref
    before = [""] * (n_ref + 1)
    i = 0
    for r, q in zip(gapped_ref, gapped_seq):
        if r != GAP:
            at[i] = q
            i += 1
        elif q != GAP:
            before[i] += q
    return at, before


def place(segment: str, width: int, mode: str = "left") -> str:
    """Put an insertion segment into a column block of fixed width.

    `left` pads on the right. Fine for the odd one- or two-residue insertion
    in a framework.

    `center` fills OUTWARD FROM BOTH ENDS and leaves the gap in the middle.
    This is the only sane way to lay out CDR3, and it is what IMGT itself
    does -- its CDR3 insertion order runs 111, 111A, 111B ... 112B, 112A, 112,
    accumulating at both ends and meeting in the middle.

    The reason is that a CDR3 has two anchors, not one: the conserved cysteine
    it starts after and the conserved tryptophan it ends before. Its first
    residues are homologous to other CDR3s' first residues and its last to
    their last; what varies is the length of the middle. Left-justifying makes
    a 10-mer and a 25-mer share their beginning and nothing else, so the whole
    column block reads as ragged expansion rather than alignment. An odd
    residue goes to the left, as IMGT does.
    """
    pad = width - len(segment)
    if pad <= 0:
        return segment[:width] if width else ""
    if mode != "center":
        return segment + GAP * pad
    left = (len(segment) + 1) // 2
    return segment[:left] + GAP * pad + segment[left:]


def _anchor_indices(numbers: list[int | None]) -> tuple[int | None, int | None]:
    """Reference indices of the residues that bracket CDR3, if numbered."""
    before = after = None
    for i, num in enumerate(numbers):
        if num == ANCHOR_BEFORE_CDR3:
            before = i
        elif num == ANCHOR_AFTER_CDR3:
            after = i
    return before, after


def consolidate_junction(at: list[str], before: list[str],
                         ref: str, c_index: int, w_index: int) -> None:
    """Force one sequence's junction into the single block between the anchors.

    Modifies `at` and `before` in place.

    The failure this exists for: when a read's conserved cysteine is missing
    or mutated, the pairwise aligner has a free column and puts the junction's
    first residue in it. That sequence's junction then starts one column to
    the LEFT of every other sequence's, and is one residue short, so it also
    gains a gap. The mirror image happens at the tryptophan, shifting a
    junction right. Two outliers in an otherwise stacked block, which is
    exactly what it looks like.

    Scoring cannot fix this. Forbidding substitution at the anchors by residue
    identity is worse -- a junction containing its own W will then capture the
    tryptophan column from across the block. So it is fixed structurally
    instead: an anchor column holds its own residue or nothing, and everything
    between the anchors is junction, by definition rather than by inference.

    There is a second way in, which the first version of this missed. When the
    cysteine is *substituted* rather than deleted, the residue count does not
    change, and because IMGT 104 takes the relaxed junction gap penalty (a
    seam takes the more permissive side), deleting the germline C is cheaper
    than mismatching it. The aligner therefore leaves the anchor column empty
    and puts the whole junction in the insertion slot *before* it -- so
    nothing lands at the anchor for the check above to catch, and that row's
    entire junction renders one block to the left with a single residue left
    behind in the CDR3 columns. Measured on 150 synthetic reads across three
    V families, three rows did exactly this.

    An insertion before an intact anchor is a genuine FR3 insertion and is
    left alone. An insertion before an *empty* anchor is a junction that has
    slid: its first residue is the substituted 104 and belongs in the anchor
    column, and the rest is junction. Checked against `cdr3_aa`, which the
    aligner never sees, this reconstructs the junction exactly.
    """
    seam = w_index
    gathered = []
    if 0 <= c_index < len(at) and c_index < len(before):
        held = at[c_index] if at[c_index] not in ("", GAP) else ""
        spilled = before[c_index]
        if held != ref[c_index] and (spilled or held):
            # Two different failures, and the insertion slot before the anchor
            # is what tells them apart:
            #
            #   spilled == ""  the cysteine is GONE and the aligner borrowed
            #                  the junction's first residue to fill the column.
            #                  Give it back; the anchor holds nothing.
            #   spilled != ""  the cysteine was SUBSTITUTED, so the residue
            #                  count never changed. Because IMGT 104 takes the
            #                  relaxed junction gap penalty, deleting the
            #                  germline C beat mismatching it, and the whole
            #                  junction slid into the slot before the anchor.
            #                  In sequence order that run starts at 104, so
            #                  its first residue belongs in the anchor column
            #                  -- substituted, but still 104 -- and the rest
            #                  is junction.
            if spilled:
                run = spilled + held
                at[c_index] = run[0]
                gathered.append(run[1:])
                before[c_index] = ""
            else:
                gathered.append(held)
                at[c_index] = GAP
    gathered.append(before[seam])
    tail = ""
    if 0 <= w_index < len(at) and at[w_index] not in ("", GAP) \
            and at[w_index] != ref[w_index]:
        tail = at[w_index]
        at[w_index] = GAP
    before[seam] = "".join(gathered) + tail


@register_backend("reference")
def _reference_backend(ref: str, seqs: dict[str, str],
                       insertion_layout: dict | None = None,
                       anchors: tuple | None = None, **kwargs) -> Msa:
    pairs = {lab: align_pair(ref, s, **kwargs) for lab, s in seqs.items()}
    decomposed = {lab: _decompose(*gp) for lab, gp in pairs.items()}

    if anchors and anchors[0] is not None and anchors[1] is not None:
        for at, before in decomposed.values():
            consolidate_junction(at, before, ref, anchors[0], anchors[1])

    n = len(ref)
    widths = [0] * (n + 1)
    for _, before in decomposed.values():
        for i, ins in enumerate(before):
            widths[i] = max(widths[i], len(ins))

    layout = insertion_layout or {}
    modes = [layout.get(i, "left") for i in range(n + 1)]

    numbers: list[int | None] = []
    for i in range(n):
        numbers += [None] * widths[i]
        numbers.append(i + 1)
    numbers += [None] * widths[n]

    def lay(at: list[str], before: list[str], ref_seq: str | None) -> str:
        out: list[str] = []
        for i in range(n):
            out.append(place(before[i], widths[i], modes[i]))
            out.append(ref_seq[i] if ref_seq else at[i])
        out.append(place(before[n], widths[n], modes[n]))
        return "".join(out)

    ref_row = lay([], [""] * (n + 1), ref)
    labels = list(decomposed)
    rows = [lay(*decomposed[lab], None) for lab in labels]
    return Msa(ref_label="germline", ref=ref_row, labels=labels, rows=rows,
               numbers=numbers, backend="reference-anchored")


def _external(binary: str, argv, ref: str, seqs: dict[str, str]) -> Msa:
    """Run an aligner that reads FASTA on disk and writes aligned FASTA."""
    exe = shutil.which(binary)
    if not exe:
        raise FileNotFoundError(f"{binary} is not on PATH")
    key = {f"s{i}": lab for i, lab in enumerate(seqs)}      # short, safe names
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.fasta"
        src.write_text(">ref\n" + ref + "\n"
                       + "".join(f">{k}\n{seqs[v]}\n" for k, v in key.items()))
        out = subprocess.run(argv(exe, str(src)), capture_output=True, text=True,
                             check=True)
    got: dict[str, list[str]] = {}
    name = None
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            got[name] = []
        elif name:
            got[name].append(line.strip())
    aligned = {k: "".join(v).upper() for k, v in got.items()}
    ref_row = aligned.pop("ref")
    numbers, seen = [], 0
    for c in ref_row:
        if c == GAP:
            numbers.append(None)
        else:
            seen += 1
            numbers.append(seen)
    labels = [key[k] for k in aligned]
    return Msa(ref_label="germline", ref=ref_row, labels=labels,
               rows=list(aligned.values()), numbers=numbers, backend=binary)


@register_backend("mafft")
def _mafft(ref: str, seqs: dict[str, str], **kwargs) -> Msa:
    return _external("mafft", lambda exe, f: [exe, "--quiet", "--auto", f], ref, seqs)


@register_backend("clustalo")
def _clustalo(ref: str, seqs: dict[str, str], **kwargs) -> Msa:
    return _external("clustalo",
                     lambda exe, f: [exe, "-i", f, "--outfmt", "fa", "--force"],
                     ref, seqs)


def available_backends() -> list[str]:
    """Backends that would actually run on this machine."""
    out = ["reference"]
    for name in ("mafft", "clustalo"):
        if shutil.which(name):
            out.append(name)
    return out


# --- getting a germline to align against ------------------------------------

def consensus_reference(germlines_seen: list[str]) -> tuple[str, int]:
    """One germline sequence from the many OAS reports for a gene.

    OAS gives a germline per *read*, so they differ in span (5' truncation)
    and occasionally in allele. The longest is the anchor; the rest are
    aligned onto it and each column takes the majority residue. Columns no
    read covered cannot appear, so the result is honest about span: it is the
    germline as far as the data saw it, not as far as IMGT defines it.
    """
    from collections import Counter

    seqs = [s for s in {("".join(c for c in g.upper() if c not in ".-"))
                        for g in germlines_seen if isinstance(g, str) and g.strip()} if s]
    if not seqs:
        return "", 0
    anchor = max(seqs, key=len)
    if len(seqs) == 1:
        return anchor, 1
    columns: list[Counter] = [Counter() for _ in anchor]
    for s in seqs:
        ga, gs = align_pair(anchor, s)
        i = 0
        for r, q in zip(ga, gs):
            if r == GAP:
                continue
            if q != GAP:
                columns[i][q] += 1
            i += 1
    out = []
    for i, col in enumerate(columns):
        out.append(col.most_common(1)[0][0] if col else anchor[i])
    return "".join(out), len(seqs)


def reference_for(rows, gene: str, germline_column: str = "germline_aa",
                  fasta_path=None, use_imgt: bool = True) -> tuple[str, str]:
    """(sequence, provenance label) for the germline this gene aligns against.

    The real IMGT reference wins whenever it can be had: the tables bundled
    with anarci first, then a GENE-DB FASTA on disk. A V gene stops at the
    conserved cysteine-plus-two (...YYCAK), so the modal J gene these reads
    call is appended -- otherwise CDR3 and FR4 would have no germline above
    them at all and every sequence's C-terminus would read as one long
    insertion.

    The D/N region between them is genuinely not germline-encoded, so it stays
    an insertion. That is the correct picture of a CDR3.

    Failing all that, the germline is reconstructed from what OAS reported for
    these very sequences, and the label says so, so no figure can imply a
    reference it did not have.
    """
    if use_imgt:
        found = _imgt_reference(rows, gene, fasta_path)
        if found:
            return found
    sub = rows[rows["v_gene"] == gene] if "v_gene" in getattr(rows, "columns", []) else rows
    if germline_column not in getattr(sub, "columns", []):
        return "", "no germline column"
    seq, n = consensus_reference(list(sub[germline_column].dropna().astype(str)))
    if not seq:
        return "", "no germline column"
    return seq, f"OAS germline consensus (n={n} distinct, read span)"


# IMGT region boundaries, by IMGT position number. These are definitions, not
# measurements -- IMGT fixes them so that the same number means the same
# structural position in every antibody.
IMGT_REGIONS: list[tuple[str, int, int]] = [
    ("FR1", 1, 26), ("CDR1", 27, 38), ("FR2", 39, 55), ("CDR2", 56, 65),
    ("FR3", 66, 104), ("CDR3", 105, 117), ("FR4", 118, 128),
]

# (gap open, gap extend) by region. Frameworks keep the standard BLOSUM62
# penalty: they are near-invariant, and a gap in one is almost always the
# aligner making a mistake rather than biology. CDR1 and CDR2 are germline-
# encoded but length-variable across genes, so they are loosened. CDR3 is the
# product of V(D)J recombination, exonuclease chew-back and N-addition -- it
# has no germline to be indel-free against, and charging framework prices to
# open it makes the aligner buy mismatches instead, smearing the junction
# across FR3 and FR4.
REGION_GAP_PENALTY: dict[str, tuple[float, float]] = {
    "FR1": (-11.0, -1.0),
    "CDR1": (-6.0, -0.8),
    "FR2": (-11.0, -1.0),
    "CDR2": (-6.0, -0.8),
    "FR3": (-11.0, -1.0),
    "CDR3": (-1.5, -0.2),
    "FR4": (-11.0, -1.0),
    None: (-11.0, -1.0),
}


# Which insertion blocks are laid out from both ends rather than left-
# justified. CDR3 is the one that matters; CDR1 and CDR2 are length-variable
# too and read better centred for the same reason.
CENTRED_REGIONS = {"CDR1", "CDR2", "CDR3"}

# The two conserved residues that bracket CDR3: IMGT 104 is the second
# cysteine of the domain's disulphide and 118 the tryptophan of the W-G-x-G
# motif. Every CDR3 definition is stated relative to them, which is why the
# junction block sits between them and why nothing else may occupy their
# columns.
ANCHOR_BEFORE_CDR3, ANCHOR_AFTER_CDR3 = 104, 118


def insertion_layout(imgt_numbers: list[int | None]) -> dict[int, str]:
    """Placement mode per seam: 'center' inside a CDR, 'left' elsewhere.

    A seam's region is the region of the position that would come NEXT after
    the residue before it -- the same rule `Msa.column_regions` uses. That
    matters at the one seam that counts: the reference has no CDR3 residues,
    so the gap between the conserved cysteine (104) and tryptophan (118) is
    flanked by FR3 and FR4. Asking which region its neighbours are in gives
    the wrong answer; asking what comes after 104 gives CDR3.
    """
    n = len(imgt_numbers)
    modes = {}
    for i in range(n + 1):
        previous = imgt_numbers[i - 1] if i > 0 else None
        nxt = imgt_numbers[i] if i < n else None
        region = (region_of(previous + 1) if previous is not None
                  else region_of(nxt))
        modes[i] = "center" if region in CENTRED_REGIONS else "left"
    return modes


def region_of(imgt: int | None) -> str | None:
    if imgt is None:
        return None
    for name, lo, hi in IMGT_REGIONS:
        if lo <= imgt <= hi:
            return name
    return None


def region_gap_profile(imgt_numbers: list[int | None],
                       penalties: dict | None = None) -> list[tuple[float, float]]:
    """(open, extend) per reference position, from IMGT region membership.

    Position i in the profile governs a gap *at* reference index i -- for an
    insertion that is the seam between residue i-1 and residue i. A seam takes
    the more permissive of the two regions it sits between, so the V/J junction
    is relaxed from its first column rather than one residue late.

    The junction is where this matters most. An IMGT V gene stops at position
    106 and a J gene starts at 115, so positions 107-114 -- the D and N-addition
    zone -- have no germline residue at all. Everything a sequence carries there
    is an insertion by construction, and it should be cheap.
    """
    penalties = penalties or REGION_GAP_PENALTY
    n = len(imgt_numbers)
    default = penalties.get(None, (-11.0, -1.0))
    at = [penalties.get(region_of(x), default) for x in imgt_numbers]

    profile = []
    for i in range(n + 1):
        left = at[i - 1] if i > 0 else None
        right = at[i] if i < n else None
        options = [x for x in (left, right) if x is not None] or [default]
        profile.append(max(options))          # least negative == most permissive
    return profile


def _modal_j(rows, column: str = "j_gene") -> str | None:
    """The J gene most of these reads call, if the frame carries the column."""
    if column not in getattr(rows, "columns", []):
        return None
    calls = rows[column].dropna().astype(str)
    calls = calls[calls.str.strip().ne("") & calls.str.lower().ne("nan")]
    return calls.value_counts().index[0] if len(calls) else None


def _imgt_reference(rows, gene: str, fasta_path=None) -> tuple[str, str] | None:
    """V (+ modal J) from the IMGT tables, with both alleles named."""
    try:
        from . import germline_db
    except Exception:                                       # noqa: BLE001
        return None

    v = None
    if fasta_path:
        try:
            v = germline_db.from_fasta(gene, fasta_path)
        except Exception:                                   # noqa: BLE001
            v = None
    if v is None:
        v = germline_db.from_anarci(gene)
    if v is None or not v.sequence:
        return None

    label = f"IMGT {v.allele or gene}"
    numbers = _numbers_for(v)
    seq = v.sequence
    j_gene = _modal_j(rows)
    if j_gene:
        j = germline_db.from_anarci(j_gene)
        if j and j.sequence:
            seq += j.sequence
            numbers += _numbers_for(j)
            label += f" + {j.allele or j_gene}"

    seq, numbers = _drop_cdr3(seq, numbers)
    _NUMBERING[seq] = numbers
    return seq, label


def _drop_cdr3(seq: str, numbers: list[int | None]) -> tuple[str, list]:
    """Remove the germline's own CDR3-span residues from the reference.

    A V gene contributes two or three residues to CDR3 (the A-R of C-A-R) and
    a J gene contributes its last few (the F-D-Y). Left in the reference they
    are anchors the aligner will match arbitrary junction residues against,
    which chops one CDR3 into three or four separate insertion blocks -- each
    then laid out independently, which is how an alignment turns into
    confetti.

    Dropping them leaves exactly one seam between the conserved cysteine
    (IMGT 104) and the conserved tryptophan (118). Every sequence's whole
    junction lands in that one block and can be laid out as a unit. Nothing
    is lost that the germline really defined: those residues ARE part of the
    CDR3, and each sequence shows its own.
    """
    keep = [(aa, num) for aa, num in zip(seq, numbers)
            if num is None or not (105 <= num <= 117)]
    return "".join(aa for aa, _ in keep), [num for _, num in keep]


# IMGT numbers for the reference, keyed by the reference sequence itself.
# reference_for's two-value return is used in several places; threading a
# third value through all of them to carry numbering that only the aligner
# wants would be worse than looking it up.
_NUMBERING: dict[str, list[int | None]] = {}


def _numbers_for(germline) -> list[int | None]:
    """IMGT number per residue of an ungapped germline, in order."""
    by_pos = getattr(germline, "by_position", None) or {}
    if not by_pos:
        return [None] * len(germline.sequence)
    return sorted(int(k) for k in by_pos)


def numbering_for(ref: str) -> list[int | None]:
    """IMGT numbers for a reference built by `reference_for`, if it had any."""
    return _NUMBERING.get(ref) or [None] * len(ref)


def regions_of(ref: str) -> list[tuple[str, int, int]]:
    """(region, first index, last index + 1) over an IMGT-numbered reference."""
    numbers = numbering_for(ref)
    spans, current, start = [], None, 0
    for i, num in enumerate(numbers):
        name = region_of(num)
        if name != current:
            if current is not None:
                spans.append((current, start, i))
            current, start = name, i
    if current is not None:
        spans.append((current, start, len(numbers)))
    return spans
