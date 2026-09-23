"""Figures for the library module.

These are regenerated from the pool, not screenshots. They double as a visual
acceptance test: if the sampler stops working, figure 3 shows it immediately.
"""
from __future__ import annotations

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


def alignment(rows: pd.DataFrame, gene: str, out: str | Path = "fig5_alignment.png",
              max_rows: int = 24, highlight_only_differences: bool = True,
              seq_column: str = "aa_gapped", germline_column: str = "germline_aa",
              show_regions: bool = True, allow_ragged: bool = False):
    """Geneious-style alignment for one germline, against the germline reference.

    Every residue letter is printed. With `highlight_only_differences` (the
    default) only positions that differ from the germline get a coloured fill;
    matching positions show the letter in muted ink on no fill, so the sequence
    stays readable while divergence is what catches the eye.

    Sequences from a single V gene arrive IMGT-gapped from OAS, so the columns
    already line up and no aligner is needed. The reference row is the germline
    OAS carries per sequence; where that is absent the observed consensus is
    used and the row is labelled as such.
    """
    sub = rows[rows["v_gene"] == gene]
    if sub.empty:
        raise ValueError(f"no sequences for {gene}")
    if seq_column not in sub.columns:
        raise KeyError(
            f"{seq_column!r} not in the frame; pass pool.fetch(...) output, "
            "which reads the parquet side where the gapped sequence lives")
    # Columns only mean anything if the sequences share them. OAS units are
    # IMGT-gapped so one germline's sequences are equal length, but reads that
    # do not cover the whole domain come back shorter. Padding those on the
    # right would slide every residue out of register and draw a confident,
    # wrong picture -- so restrict to the modal length and say so.
    lengths = sub[seq_column].astype(str).str.len()
    modal = int(lengths.mode().iloc[0])
    n_total = len(sub)
    aligned = sub[lengths == modal]
    dropped = n_total - len(aligned)
    if dropped and not allow_ragged:
        sub = aligned
    elif dropped and allow_ragged:
        sub = sub.assign(**{seq_column: sub[seq_column].astype(str).str.ljust(modal, ".")})

    sub = sub.head(max_rows)
    seqs = [str(s) for s in sub[seq_column]]
    labels = [str(s) for s in sub["seq_id"]]
    width = max(len(s) for s in seqs)
    seqs = [s.ljust(width, ".") for s in seqs]

    ref, ref_label = reference_row(sub, seqs, germline_column)
    ref = ref.ljust(width, ".")[:width]

    spans = region_spans(sub.iloc[0], seqs[0]) if show_regions else []

    n = len(seqs)
    fig_w = min(24, max(8, width * 0.115))
    band_h = 0.9 if spans else 0.0
    fig, ax = plt.subplots(figsize=(fig_w, 0.26 * (n + 2) + 1.5 + band_h * 0.25), dpi=150)
    show_letters = width <= 260

    for r, (seq, label) in enumerate([(ref, ref_label)] + list(zip(seqs, labels))):
        y = n - r
        for c, aa in enumerate(seq):
            differs = r == 0 or aa != ref[c]
            if aa in ".-":
                if r == 0:
                    ax.add_patch(plt.Rectangle((c, y), 1, 1, facecolor=GAP_COLOR,
                                               edgecolor="white", linewidth=0.25))
                continue
            if differs or not highlight_only_differences:
                color = AA_COLOR.get(aa, UNKNOWN_COLOR)
                ax.add_patch(plt.Rectangle((c, y), 1, 1, facecolor=color,
                                           edgecolor="white", linewidth=0.25))
                ink = _text_on(color)
            else:
                ink = "#9c9a92"          # matches germline: present, recessive
            if show_letters:
                ax.text(c + 0.5, y + 0.5, aa, ha="center", va="center",
                        fontsize=5.5, color=ink)
        ax.text(-1.5, y + 0.5, label[:22], ha="right", va="center",
                fontsize=7, color=INK if r == 0 else MUTED,
                fontweight="bold" if r == 0 else "normal")

    top = n + 1
    for name, start, end in spans:
        is_cdr = name.startswith("CDR")
        ax.add_patch(plt.Rectangle((start, top + 0.15), end - start, 0.62,
                                   facecolor=REGION_BAND[is_cdr],
                                   edgecolor="white", linewidth=0.6))
        ax.text((start + end) / 2, top + 0.46, name, ha="center", va="center",
                fontsize=7, color="#4a4a46",
                fontweight="bold" if is_cdr else "normal")

    ax.set_xlim(-0.5, width + 0.5)
    ax.set_ylim(-0.4, top + 1.1)
    ax.set_xticks(range(0, width, 10))
    ax.set_xticklabels(range(0, width, 10), fontsize=7, color=MUTED)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, length=2)
    mode = "differences highlighted vs" if highlight_only_differences else "full colour vs"
    note = ""
    if dropped:
        note = (f"; {dropped} of {n_total} padded to fit" if allow_ragged
                else f"; {dropped} of {n_total} excluded as not column-aligned")
    ax.set_title(f"{gene} — {mode} {ref_label}; IMGT regions above{note}",
                 color=INK, fontsize=11, loc="left", pad=10)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return Path(out)


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
