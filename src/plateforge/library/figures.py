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
    ]
