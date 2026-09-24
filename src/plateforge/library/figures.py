"""Figures for the library module.

These are regenerated from the pool, not screenshots. They double as a visual
acceptance test: if the sampler stops working, figure 3 shows it immediately.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from . import diversity, germlines

# Validated categorical palette (CVD-checked; see docs/formats/README.md).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#d8d7d2"


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def gene_composition(pool: pd.DataFrame, picked: pd.DataFrame,
                     panel: germlines.Panel | None = None, out: str | Path = "fig1_genes.png"):
    """Pool composition beside what the sampler drew, as proportions."""
    panel = panel or germlines.DEFAULT
    genes = panel.genes
    pool_frac = [(pool["v_gene"] == g).mean() for g in genes]
    pick_frac = [(picked["v_gene"] == g).mean() for g in genes]
    target = [panel.weight_of(g) / sum(panel.weight_of(x) for x in genes) for g in genes]

    fig, ax = plt.subplots(figsize=(7, 3.6), dpi=150)
    x = range(len(genes))
    w = 0.27
    for i, (vals, label, color) in enumerate([
        (pool_frac, "pool", SERIES[0]),
        (pick_frac, "sampled", SERIES[1]),
        (target, "panel target", SERIES[2]),
    ]):
        pos = [p + (i - 1) * w for p in x]
        bars = ax.bar(pos, vals, width=w - 0.02, color=color, label=label, zorder=3)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.0%}",
                    ha="center", va="bottom", fontsize=8, color=MUTED)
    ax.set_xticks(list(x))
    ax.set_xticklabels(genes, color=INK)
    ax.set_ylabel("fraction of sequences", color=MUTED, fontsize=9)
    ax.set_ylim(0, max(max(pool_frac), max(pick_frac), max(target)) * 1.15)
    ax.set_title("Germline composition: pool, sample, and panel target",
                 color=INK, fontsize=11, loc="left", pad=28)
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED, ncol=3,
              loc="lower right", bbox_to_anchor=(1, 1.02))
    _style(ax)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return Path(out)


def cdr3_lengths(pool: pd.DataFrame, picked: pd.DataFrame,
                 out: str | Path = "fig2_cdr3.png"):
    """CDRH3 length per germline, with the sampled members marked."""
    genes = sorted(pool["v_gene"].dropna().unique())
    fig, ax = plt.subplots(figsize=(7, 3.6), dpi=150)
    for i, g in enumerate(genes):
        sub = pool[pool["v_gene"] == g]["cdr3_len"]
        parts = ax.violinplot([sub], positions=[i], widths=0.7, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(SERIES[i % len(SERIES)])
            body.set_alpha(0.30)
        sel = picked[picked["v_gene"] == g]["cdr3_len"]
        if len(sel):
            ax.scatter([i] * len(sel), sel, s=34, color=SERIES[i % len(SERIES)],
                       edgecolor="white", linewidth=1.4, zorder=4)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, color=INK)
    ax.set_ylabel("CDRH3 length (aa)", color=MUTED, fontsize=9)
    ax.set_title("CDRH3 length: pool distribution, sampled members marked",
                 color=INK, fontsize=11, loc="left", pad=12)
    _style(ax)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return Path(out)


def identity_check(pool: pd.DataFrame, picked: pd.DataFrame, n_random: int = 30,
                   seed: int = 0, out: str | Path = "fig3_identity.png"):
    """Worst pairwise CDRH3 identity: the sample against random draws.

    This is the figure that proves the sampler is doing something.
    """
    import random

    def worst(xs: list[str]) -> float:
        return max(diversity.identity(a, b)
                   for i, a in enumerate(xs) for b in xs[i + 1:])

    k = len(picked)
    rng = random.Random(seed)
    seqs = list(pool["cdr3_aa"])
    randoms = [worst(rng.sample(seqs, k)) for _ in range(n_random)]
    chosen = worst(list(picked["cdr3_aa"]))

    fig, ax = plt.subplots(figsize=(7, 3.2), dpi=150)
    ax.hist(randoms, bins=12, color=SERIES[0], alpha=0.85, zorder=3,
            label=f"random draws of {k} (n={n_random})")
    ax.axvline(chosen, color=SERIES[1], linewidth=2.5, zorder=4,
               label=f"diversity sampler ({chosen:.2f})")
    ax.set_xlabel("worst pairwise CDRH3 identity in the set", color=MUTED, fontsize=9)
    ax.set_ylabel("draws", color=MUTED, fontsize=9)
    ax.set_title("Lower is more diverse", color=INK, fontsize=11, loc="left", pad=28)
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED, ncol=2,
              loc="lower left", bbox_to_anchor=(0, 1.02))
    _style(ax)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return Path(out)


def all_figures(pool: pd.DataFrame, picked: pd.DataFrame, out_dir: str | Path,
                panel: germlines.Panel | None = None) -> list[Path]:
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    return [
        gene_composition(pool, picked, panel, d / "fig1_germline_composition.png"),
        cdr3_lengths(pool, picked, d / "fig2_cdr3_lengths.png"),
        identity_check(pool, picked, out=d / "fig3_identity_check.png"),
        family_distribution(pool, d / "fig4_family_distribution.png"),
    ]


# --- composition across every family in the pool ---------------------------

# One hue per residue, banded by physico-chemical property so families read as
# colour neighbourhoods: hydrophobics are blues, aromatics teals, positives
# reds, negatives purples, polars greens, and C/G/P get their own.
#
# Twenty categories cannot be told apart by colour alone -- that is a limit of
# perception, not of palette design. These hues were chosen by maximising the
# minimum pairwise OKLab distance under normal, deuteranopic and protanopic
# vision within each property band (worst pair ~5.0 dE100, every hue above the
# gray chroma floor). The residue letter is drawn in every cell, so colour
# assists recognition and never carries it alone.
AA_COLOR = {
    "A": "#a3cfec", "V": "#95bbe9", "L": "#86a3e6", "I": "#7888e3", "M": "#6a6ae1",
    "F": "#a5ddcb", "W": "#74cabe", "Y": "#45b5b7", "H": "#337787",
    "K": "#ed3157", "R": "#a0520e",
    "D": "#6b03f3", "E": "#f803ae",
    "N": "#c5e1ac", "Q": "#94d07b", "S": "#57be4a", "T": "#34913a",
    "C": "#c84386", "G": "#f89028", "P": "#ebde61",
}
AA_GROUPS = {
    "hydrophobic": "AVLIM",
    "aromatic": "FWYH",
    "positive": "KR",
    "negative": "DE",
    "polar": "NQST",
    "cysteine": "C",
    "glycine": "G",
    "proline": "P",
}
GAP_COLOR = "#f0efe9"
UNKNOWN_COLOR = "#c9c8c1"
# Three greys, in decreasing weight: the germline reference, a residue that
# matches it, and a gap. The germline is the darkest because it is the thing
# being read against; nothing on that row is ever filled.
REF_INK = "#3a3935"
MATCH_INK = "#9c9a92"
GAP_INK = "#b0aea6"


def _text_on(hex_color: str) -> str:
    """Dark or light ink, whichever is readable on this fill."""
    r, g, b = (int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#1a1a18" if lum > 0.48 else "#ffffff"


def family_distribution(pool: pd.DataFrame, out: str | Path = "fig4_families.png",
                        top: int = 15):
    """Every V gene present in the pool, grouped by family.

    fig1 answers "did the sampler hit its quota". This answers "what is
    actually in the pool", which is the question you ask of a real download.
    """
    counts = pool["v_gene"].value_counts().head(top)
    fams = [germlines.family(g) for g in counts.index]
    order = sorted(set(fams))
    colors = [SERIES[order.index(f) % len(SERIES)] for f in fams]

    fig, ax = plt.subplots(figsize=(7, 0.34 * len(counts) + 1.6), dpi=150)
    ypos = range(len(counts))
    ax.barh(list(ypos), counts.values, color=colors, height=0.72, zorder=3)
    for y, v in zip(ypos, counts.values):
        ax.text(v + counts.max() * 0.012, y, f"{v:,}", va="center",
                fontsize=8, color=MUTED)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(counts.index, color=INK, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("sequences in pool", color=MUTED, fontsize=9)
    ax.set_xlim(0, counts.max() * 1.12)
    ax.set_title("V gene usage in the pool", color=INK, fontsize=11, loc="left", pad=28)
    handles = [plt.Rectangle((0, 0), 1, 1, color=SERIES[i % len(SERIES)])
               for i in range(len(order))]
    ax.legend(handles, order, frameon=False, fontsize=9, labelcolor=MUTED,
              ncol=len(order), loc="lower left", bbox_to_anchor=(0, 1.02))
    _style(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.yaxis.grid(False)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return Path(out)


# --- alignment / divergence view -------------------------------------------

def _consensus(seqs: list[str]) -> str:
    from collections import Counter
    width = max(len(s) for s in seqs)
    padded = [s.ljust(width, ".") for s in seqs]
    return "".join(Counter(col).most_common(1)[0][0] for col in zip(*padded))


def reference_row(rows: pd.DataFrame, seqs: list[str],
                  germline_column: str = "germline_aa") -> tuple[str, str]:
    """The row everything is compared against, and a label for it.

    Prefers the germline reference OAS ships per sequence (IgBlast's own, and
    already column-aligned). Falls back to the observed consensus only when the
    data unit carries no germline column, and says so in the label so a figure
    can never silently imply a germline it did not have.
    """
    if germline_column in rows.columns:
        germ = rows[germline_column].dropna().astype(str)
        germ = germ[germ.str.strip().ne("") & germ.str.lower().ne("nan")]
        if len(germ):
            uniq = germ.unique()
            if len(uniq) == 1:
                return uniq[0], "germline"
            # Several alleles in one gene: use the most common, and say so.
            top = germ.value_counts().index[0]
            return top, f"germline (modal of {len(uniq)})"
    return _consensus(seqs), "consensus (no germline column)"


REGION_ORDER = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa"]
REGION_LABEL = {"fwr1_aa": "FR1", "cdr1_aa": "CDR1", "fwr2_aa": "FR2",
                "cdr2_aa": "CDR2", "fwr3_aa": "FR3", "cdr3_aa": "CDR3"}
REGION_BAND = {True: "#ded9f0", False: "#eceae3"}   # CDR shaded, FR plain


def region_spans(row: pd.Series, gapped: str) -> list[tuple[str, int, int]]:
    """Locate the IMGT regions as column ranges in a gapped sequence.

    OAS gives each region as its own amino acid string (fwr1_aa, cdr1_aa, ...).
    Rather than trusting cumulative lengths, each region is *searched for* in
    the ungapped sequence starting after the previous one, then mapped back to
    gapped columns. That survives the places where region strings overlap or
    do not tile the sequence exactly.
    """
    col_of, ungapped = [], []
    for col, ch in enumerate(gapped):
        if ch not in ".-":
            col_of.append(col)
            ungapped.append(ch)
    seq = "".join(ungapped)

    spans, cursor = [], 0
    for key in REGION_ORDER:
        frag = row.get(key)
        if not isinstance(frag, str) or not frag.strip():
            continue
        i = seq.find(frag, cursor)
        if i < 0:
            i = seq.find(frag)
        if i < 0:
            continue
        spans.append((REGION_LABEL[key], col_of[i], col_of[i + len(frag) - 1] + 1))
        cursor = i + len(frag)
    return spans


def build_msa(rows: pd.DataFrame, gene: str, *, max_rows: int = 50,
              seq_column: str = "aa_gapped", germline_column: str = "germline_aa",
              backend: str = "reference", fasta_path=None,
              order_by: str | None = "cdr3_len"):
    """Align one V gene's sequences against that gene's germline.

    Returns (Msa, source_rows) where source_rows is the frame the alignment was
    built from, in the same order as the Msa's rows -- callers need it for
    region bands and for labelling.

    Nothing is excluded for being short. A 5' truncated read aligns with
    leading gaps, which is what the gaps are for.
    """
    from . import msa as _msa

    sub = rows[rows["v_gene"] == gene] if "v_gene" in rows.columns else rows
    if sub.empty:
        raise ValueError(f"no sequences for {gene}")
    if seq_column not in sub.columns:
        raise KeyError(
            f"{seq_column!r} not in the frame; pass pool.fetch(...) output, "
            "which reads the parquet side where the gapped sequence lives")

    ref, ref_source = _msa.reference_for(sub, gene, germline_column, fasta_path)
    if not ref:
        return None, sub

    if order_by and order_by in sub.columns:
        sub = sub.sort_values([order_by, "seq_id"], kind="stable")
    sub = sub.head(max_rows)

    seqs, seen = {}, set()
    for _, row in sub.iterrows():
        label = str(row["seq_id"])
        seq = str(row.get(seq_column) or "")
        if not seq.strip() or label in seen:
            continue
        seen.add(label)
        seqs[label] = seq
    if not seqs:
        return None, sub

    try:
        m = _msa.build(ref, seqs, ref_label="germline", ref_source=ref_source,
                       backend=backend)
    except (FileNotFoundError, subprocess.CalledProcessError, KeyError):
        m = _msa.build(ref, seqs, ref_label="germline", ref_source=ref_source,
                       backend="reference")
    ordered = sub.set_index("seq_id").loc[m.labels].reset_index()
    return m, ordered


def draw_msa(ax, m, source_rows=None, *, highlight_only_differences=True,
             show_regions=True, show_letters=None, label_size=7,
             letter_size=5.5, row_labels=True, number_step=None):
    """Draw an Msa: germline on top and numbered, variants below, unnumbered.

    The numbering belongs to the reference and only to the reference. An
    inserted column is not a germline position, so it gets no number -- and a
    variant's own residues are never numbered, because numbering them would
    assert a position the germline does not define.

    Colour means one thing and one thing only: *this differs from the
    germline*. So the germline row itself is never filled -- it is grey letters
    on white, the baseline the eye reads everything else against. A filled cell
    below it is a substitution, or an insertion the germline has no residue
    for, which is why an inserted residue is coloured even though there is
    nothing above it to differ from.
    """
    width, n = m.width, m.n_sequences
    if show_letters is None:
        show_letters = width <= 300

    rows_to_draw = [(m.ref_label, m.ref)] + list(zip(m.labels, m.rows))
    for r, (label, seq) in enumerate(rows_to_draw):
        y = n - r
        reference_row = r == 0
        for c, aa in enumerate(seq):
            if aa in ".-":
                # A gap on the reference row is an insertion column: no
                # germline residue exists there, so tint it like the rest of
                # the column rather than marking the germline as divergent.
                ax.add_patch(plt.Rectangle((c, y), 1, 1, facecolor=GAP_COLOR,
                                           edgecolor="white", linewidth=0.25))
                if show_letters:
                    ax.text(c + 0.5, y + 0.5, "-", ha="center", va="center",
                            fontsize=letter_size, color=GAP_INK)
                continue
            differs = (not reference_row) and aa != m.ref[c]
            if differs or (not highlight_only_differences and not reference_row):
                color = AA_COLOR.get(aa, UNKNOWN_COLOR)
                ax.add_patch(plt.Rectangle((c, y), 1, 1, facecolor=color,
                                           edgecolor="white", linewidth=0.25))
                ink = _text_on(color)
            else:
                ink = REF_INK if reference_row else MATCH_INK
            if show_letters:
                ax.text(c + 0.5, y + 0.5, aa, ha="center", va="center",
                        fontsize=letter_size, color=ink,
                        fontweight="bold" if reference_row else "normal")
        if row_labels:
            ax.text(-1.5, y + 0.5, str(label)[:22], ha="right", va="center",
                    fontsize=label_size, color=INK if r == 0 else MUTED,
                    fontweight="bold" if r == 0 else "normal")

    top = n + 1
    spans = []
    if show_regions and source_rows is not None and len(source_rows):
        spans = region_spans(source_rows.iloc[0], m.rows[0])
    for name, start, end in spans:
        is_cdr = name.startswith("CDR")
        ax.add_patch(plt.Rectangle((start, top + 0.15), end - start, 0.62,
                                   facecolor=REGION_BAND[is_cdr],
                                   edgecolor="white", linewidth=0.6))
        ax.text((start + end) / 2, top + 0.46, name, ha="center", va="center",
                fontsize=label_size, color="#4a4a46",
                fontweight="bold" if is_cdr else "normal")

    ax.set_xlim(-0.5, width + 0.5)
    ax.set_ylim(-0.4, top + (1.1 if spans else 0.2))

    step = number_step or (10 if width <= 200 else 20)
    ticks, tick_labels = [], []
    for c, num in enumerate(m.numbers):
        if num is not None and (num % step == 0 or num == 1):
            ticks.append(c + 0.5)
            tick_labels.append(str(num))
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=label_size, color=MUTED)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2)
    return m


def alignment(rows: pd.DataFrame, gene: str, out: str | Path = "fig5_alignment.png",
              max_rows: int = 50, highlight_only_differences: bool = True,
              seq_column: str = "aa_gapped", germline_column: str = "germline_aa",
              show_regions: bool = True, backend: str = "reference",
              fasta_path=None, also_write_fasta: bool = True):
    """Alignment figure for one V gene, germline at the top.

    Every sequence calling this gene is aligned to the gene's germline with
    gaps (see `library.msa`). The germline is the first row and carries the
    numbering; nothing below it is numbered. With
    `highlight_only_differences` only positions that depart from the germline
    get a coloured fill, so divergence is what catches the eye.

    Writes the alignment beside the figure as FASTA, so it can be opened in
    Geneious, Jalview, or anything else, and so the picture is checkable.
    """
    m, source = build_msa(rows, gene, max_rows=max_rows, seq_column=seq_column,
                          germline_column=germline_column, backend=backend,
                          fasta_path=fasta_path)
    out = Path(out)
    if m is None:
        fig, ax = plt.subplots(figsize=(9, 2.2), dpi=150)
        ax.axis("off")
        ax.text(0.02, 0.6, f"No alignment for {gene}.", fontsize=12, color=INK)
        ax.text(0.02, 0.32, f"Needs a germline reference ({germline_column}) and "
                            f"sequences in {seq_column}.", fontsize=9, color=MUTED)
        fig.savefig(out, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        return out

    width, n = m.width, m.n_sequences
    fig_w = min(26, max(8, width * 0.115))
    fig, ax = plt.subplots(figsize=(fig_w, 0.26 * (n + 2) + 1.8), dpi=150)
    draw_msa(ax, m, source, highlight_only_differences=highlight_only_differences,
             show_regions=show_regions)
    ax.set_xlabel("germline residue number — insertions are not germline positions",
                  fontsize=8, color=MUTED)
    mode = ("colour = differs from germline" if highlight_only_differences
            else "every residue coloured")
    ax.set_title(f"{gene} — {n} sequences vs germline [{m.ref_source}]; "
                 f"{mode}; {m.backend}",
                 color=INK, fontsize=11, loc="left", pad=10)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)

    if also_write_fasta:
        out.with_suffix(".fasta").write_text(m.as_fasta())
        out.with_suffix(".aln.txt").write_text(m.as_text())
    return out


def alignments_per_gene(rows: pd.DataFrame, out_dir: str | Path,
                        genes: list[str] | None = None, *, max_rows: int = 50,
                        prefix: str = "alignment", **kwargs) -> dict[str, Path]:
    """One alignment per V gene present. This is the per-family deliverable."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if genes is None:
        genes = list(rows["v_gene"].dropna().value_counts().index)
    made = {}
    for gene in genes:
        if not len(rows[rows["v_gene"] == gene]):
            continue
        made[gene] = alignment(rows, gene, out_dir / f"{prefix}_{gene}.png",
                               max_rows=max_rows, **kwargs)
    return made


def aa_legend(out: str | Path = "fig5b_aa_legend.png"):
    """Key for the residue colouring used by `alignment`, grouped by property."""
    fig, ax = plt.subplots(figsize=(9, 1.5), dpi=150)
    x = 0.0
    for group, residues in AA_GROUPS.items():
        ax.text(x, 1.30, group, fontsize=8, color=MUTED, ha="left")
        for i, aa in enumerate(residues):
            color = AA_COLOR[aa]
            ax.add_patch(plt.Rectangle((x + i * 0.62, 0.55), 0.56, 0.58,
                                       facecolor=color, edgecolor="white",
                                       linewidth=0.6))
            ax.text(x + i * 0.62 + 0.28, 0.84, aa, ha="center", va="center",
                    fontsize=9, color=_text_on(color))
        x += max(len(residues), 1) * 0.62 + 0.75
    ax.set_xlim(-0.2, x)
    ax.set_ylim(0.2, 1.7)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return Path(out)


# --- one figure summarising the input data ---------------------------------

def summary_panel(pool_df: pd.DataFrame, picked: pd.DataFrame,
                  alignment_rows: pd.DataFrame | None = None,
                  gene: str | None = None,
                  out: str | Path = "input_summary.png",
                  panel: germlines.Panel | None = None,
                  max_alignment_rows: int = 50):
    """Germline usage, CDRH3 lengths and an example alignment in one figure.

    The alignment is drawn as a difference map against the germline, the same
    way the per-family figures are, so one reading applies to every alignment
    this repo produces: grey is germline, colour is departure from it.

    50 sequences by default: enough to stand in for the list that will be
    marched through the downstream steps, few enough that the rows are still
    distinguishable at figure size.
    """
    from matplotlib import gridspec

    panel = panel or germlines.DEFAULT
    genes = [g for g in panel.genes if (pool_df["v_gene"] == g).any()]
    gene = gene or (picked["v_gene"].value_counts().index[0] if len(picked) else genes[0])

    have_aln = alignment_rows is not None and len(alignment_rows)
    n_rows = max_alignment_rows
    # Height scales with the alignment: 50 rows need room to stay legible, and
    # a fixed figure height would either squash them or leave a field of white
    # under a 10-row one.
    aln_h = 0.20 * (n_rows + 3) if have_aln else 0
    fig = plt.figure(figsize=(17, 3.4 + aln_h), dpi=150)
    gs = gridspec.GridSpec(2 if have_aln else 1, 2,
                           height_ratios=[3.0, aln_h] if have_aln else [1],
                           hspace=0.18, wspace=0.16,
                           left=0.05, right=0.985, top=0.93, bottom=0.05)

    ax1 = fig.add_subplot(gs[0, 0])
    pool_frac = [(pool_df["v_gene"] == g).mean() for g in genes]
    pick_frac = [(picked["v_gene"] == g).mean() for g in genes]
    x = range(len(genes))
    for i, (vals, label, color) in enumerate([(pool_frac, "pool", SERIES[0]),
                                              (pick_frac, "picked", SERIES[1])]):
        pos = [p + (i - 0.5) * 0.36 for p in x]
        bars = ax1.bar(pos, vals, width=0.34, color=color, label=label, zorder=3)
        for b, v in zip(bars, vals):
            ax1.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.0%}",
                     ha="center", va="bottom", fontsize=8, color=MUTED)
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(genes, color=INK, fontsize=9)
    ax1.set_ylabel("fraction", color=MUTED, fontsize=9)
    ax1.set_ylim(0, max(pool_frac + pick_frac) * 1.22)
    ax1.set_title("Germline usage", color=INK, fontsize=11, loc="left", pad=22)
    ax1.legend(frameon=False, fontsize=9, labelcolor=MUTED, ncol=2,
               loc="lower left", bbox_to_anchor=(0, 1.0))
    _style(ax1)

    ax2 = fig.add_subplot(gs[0, 1])
    for i, g in enumerate(genes):
        sub = pool_df[pool_df["v_gene"] == g]["cdr3_len"]
        parts = ax2.violinplot([sub], positions=[i], widths=0.7, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(SERIES[i % len(SERIES)])
            body.set_alpha(0.30)
        sel = picked[picked["v_gene"] == g]["cdr3_len"]
        if len(sel):
            ax2.scatter([i] * len(sel), sel, s=16, color=SERIES[i % len(SERIES)],
                        edgecolor="white", linewidth=0.8, zorder=4)
    ax2.set_xticks(range(len(genes)))
    ax2.set_xticklabels(genes, color=INK, fontsize=9)
    ax2.set_ylabel("CDRH3 length (aa)", color=MUTED, fontsize=9)
    ax2.set_title("CDRH3 length — pool, picked marked", color=INK, fontsize=11,
                  loc="left", pad=22)
    _style(ax2)

    if have_aln:
        ax3 = fig.add_subplot(gs[1, :])
        m, source = build_msa(alignment_rows, gene, max_rows=n_rows)
        if m is None:
            # Never hand back a panel of bare axes: say why it is empty.
            ax3.text(0.5, 0.5,
                     f"no germline reference for {gene}, so no alignment\n"
                     "(needs germline_aa on these rows — see pool.schema_drift)",
                     ha="center", va="center", fontsize=10, color=MUTED,
                     transform=ax3.transAxes)
            ax3.set_xticks([]); ax3.set_yticks([])
            for side in ("top", "right", "left", "bottom"):
                ax3.spines[side].set_visible(False)
        else:
            draw_msa(ax3, m, source, highlight_only_differences=True,
                     show_letters=True, letter_size=4.2, label_size=6)
            ax3.set_title(f"{gene} — {m.n_sequences} sequences vs germline "
                          f"[{m.ref_source}]; colour = differs from germline",
                          color=INK, fontsize=11, loc="left", pad=20)
            ax3.set_xlabel("germline residue number — only the germline is numbered",
                           fontsize=8, color=MUTED)

    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return Path(out)


# --- what we record about each clone ---------------------------------------

TABLE_COLUMNS = [
    ("well", "Well", 5),
    ("clone_id", "Clone ID", 14),
    ("v_gene", "V gene", 9),
    ("j_gene", "J gene", 7),
    ("cdr3_aa", "CDRH3", 16),
    ("insert_aa_length", "VH aa", 6),
    ("order_length", "Order bp", 9),
    ("order_gc", "GC", 6),
    ("fusion_site_5p", "OH 5'", 6),
    ("fusion_site_3p", "OH 3'", 6),
    ("aa_length", "IgG aa", 7),
    ("germline_source", "Germline src", 13),
    ("synthesis_warnings", "Warnings", 10),
]


# Clone ids share a long prefix, so truncating from the front prints the same
# string in every row. Their unique tail is what distinguishes them.
KEEP_TAIL = {"clone_id", "seq_id"}


def _shorten(value, width: int, keep_tail: bool = False) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    text = f"{value:.3f}" if isinstance(value, float) else str(value)
    if len(text) <= width:
        return text
    return ("…" + text[-(width - 1):]) if keep_tail else text[:width - 1] + "…"


def metadata_table(clones: pd.DataFrame, n: int = 10,
                   out: str | Path = "fig6_metadata.png",
                   columns: list | None = None, title: str | None = None):
    """The clone metadata table as a figure: what is recorded, per clone.

    A CSV with 38 columns is not something anyone reads. This is the same data
    with the long fields truncated, so the *shape* of the record is visible at
    a glance -- which is the question people actually have when they ask what
    the pipeline tracks.
    """
    columns = columns or [c for c in TABLE_COLUMNS if c[0] in clones.columns]
    rows = clones.head(n)
    n_cols, n_rows = len(columns), len(rows)

    widths = [max(len(label), width) for _, label, width in columns]
    total = sum(widths)
    fig_w = min(22, max(9, total * 0.115))
    fig, ax = plt.subplots(figsize=(fig_w, 0.32 * (n_rows + 3)), dpi=150)
    ax.axis("off")
    ax.set_xlim(0, total)
    ax.set_ylim(-1.2, n_rows + 1.4)

    edges, at = [], 0
    for w in widths:
        edges.append(at)
        at += w

    for c, ((_, label, _), x, w) in enumerate(zip(columns, edges, widths)):
        ax.add_patch(plt.Rectangle((x, n_rows), w, 1.0, facecolor="#eceae3",
                                   edgecolor="white", linewidth=1.2))
        ax.text(x + 0.4, n_rows + 0.5, label, ha="left", va="center",
                fontsize=8, color=INK, fontweight="bold")

    for r in range(n_rows):
        y = n_rows - 1 - r
        if r % 2 == 0:
            ax.add_patch(plt.Rectangle((0, y), total, 1.0, facecolor="#faf9f6",
                                       edgecolor="none"))
        row = rows.iloc[r]
        for (key, _, width), x in zip(columns, edges):
            ax.text(x + 0.4, y + 0.5,
                    _shorten(row.get(key), width, key in KEEP_TAIL),
                    ha="left", va="center", fontsize=7.2, color=MUTED,
                    family="DejaVu Sans Mono")

    ax.plot([0, total], [n_rows, n_rows], color=GRID, linewidth=0.8)
    ax.plot([0, total], [0, 0], color=GRID, linewidth=0.8)
    ax.set_title(title or f"Clone metadata — first {n_rows} of {len(clones)} wells "
                          f"({len(clones.columns)} columns recorded per clone)",
                 color=INK, fontsize=11, loc="left", pad=12)
    ax.text(0, -0.9, "Long fields truncated for display; clones.csv holds the "
                     "full values.", fontsize=7.5, color=MUTED, ha="left")
    fig.tight_layout()
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return Path(out)


def stamp(out_dir: str | Path, run_id: str = "", note: str = "") -> Path:
    """Record which run wrote the figures in this directory.

    "Did my figures update?" is otherwise unanswerable by looking at them, and
    a figure whose output path moved between versions sits there forever
    looking current. This makes the answer one `cat` away.
    """
    from datetime import datetime, timezone
    from ..core import version

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "FIGURES_FROM.txt"
    lines = [
        f"written   {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"code      {version.code_version()}",
        f"run       {run_id or '(not recorded)'}",
    ]
    if note:
        lines.append(f"note      {note}")
    lines.append("")
    lines.append("Any PNG in this directory older than the timestamp above was")
    lines.append("left behind by an earlier version and is NOT current.")
    path.write_text("\n".join(lines) + "\n")
    return path


def stale(out_dir: str | Path, made: list) -> list[Path]:
    """PNGs in this directory that the current run did not write."""
    out_dir = Path(out_dir)
    if not out_dir.exists():
        return []
    fresh = {Path(p).resolve() for p in made}
    return sorted(p for p in out_dir.glob("*.png") if p.resolve() not in fresh)
