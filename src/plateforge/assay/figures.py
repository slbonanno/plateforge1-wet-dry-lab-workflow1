"""Figures for plates: what was read, and what happened to the samples.

Three pictures, each answering a question that a table answers badly.

    elisa_pair      is this clone binding the target, or everything?
    clone_journey   where has this sample been, and what was it each time?
    pipeline_map    what did the whole run actually do?

Anything drawn from `assay.elisa` is simulated and says so on its face, not
only in its manifest.

Colour follows `core.style`, so these and the `library` figures read as one
system. The plate heatmaps use the sequential ramp -- one hue, light to dark
-- because a well's value is a magnitude. A categorical scheme on a plate,
and the rainbow ramps plate software tends to ship, both invent boundaries
that are not in the data.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..core import style, wells as wellmod

CHANNEL = "od450"
TARGET, CONTROL = "od450_target", "od450_control"

# One colour per content kind, fixed. The sample is the entity being tracked,
# and what it is at each hop is the thing worth encoding -- never its position
# in the chain, which the reader can already see.
CONTENT_COLOUR = {
    "dna": style.SEQUENTIAL[1],
    "cells": style.SEQUENTIAL[2],
    "supernatant": style.SEQUENTIAL[3],
    "beads": style.SEQUENTIAL[5],
    "eluate": style.SEQUENTIAL[4],
    "igg_prep": style.SEQUENTIAL[6],
    "dilution": style.SEQUENTIAL[4],
    "assay": style.SEQUENTIAL[6],
}


def _grid(frame: pd.DataFrame, plate_format: int, column: str,
          samples_only: bool = False) -> np.ndarray:
    """The plate as an array, with wells that hold nothing left as NaN.

    `samples_only` matters more than it looks. A 96-well plate rarely holds
    96 samples, and an empty well still returns a number -- the blank. Drawn
    on the same ramp as the samples it reads as a very weak sample, so a
    12-sample plate looks like a 96-sample plate that mostly failed. NaN
    instead, and the colormap paints those wells as "nothing here".
    """
    n_rows, n_cols = wellmod.dims(plate_format)
    out = np.full((n_rows, n_cols), np.nan)
    for _, row in frame.iterrows():
        well = wellmod.normalize(str(row["well"]))
        if samples_only and "role" in frame.columns and row["role"] != "sample":
            continue
        value = row.get(column)
        if value is not None and value == value:
            out[ord(well[0]) - ord("A"), int(well[1:]) - 1] = float(value)
    return out


def _plate_axes(ax, plate_format: int) -> None:
    n_rows, n_cols = wellmod.dims(plate_format)
    ax.set_xticks(range(n_cols))
    # Two digits, always (rule 3): the axis of a plate figure reads the same
    # way as every well ID in the repo, so A01 is findable by eye.
    ax.set_xticklabels([f"{c:02d}" for c in range(1, n_cols + 1)],
                       fontsize=6, color=style.MUTED)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([wellmod.row_label(r + 1) for r in range(n_rows)],
                       fontsize=6, color=style.MUTED)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)


def elisa_pair(pair, out: str | Path = "elisa_pair.png",
               plate_format: int = 96, mark_truth: bool = True):
    """One experiment: the same clones on target and on control.

    Both plates share one colour scale, because the entire question is which
    wells are higher *here than there*, and two scales would answer it for
    you. True binders are ringed when the truth is known -- a simulated pair
    knows, a real one does not, and the ring is what makes it obvious at a
    glance whether a calling rule would have found them.

    The difference panel is the one people actually read: target minus
    control, on a diverging scale centred at zero, so a well that is higher
    on the control plate is visibly the wrong sign rather than merely small.
    """
    out = Path(out)
    target = _grid(pair.target, plate_format, CHANNEL, samples_only=True)
    control = _grid(pair.control, plate_format, CHANNEL, samples_only=True)
    delta = target - control
    ceiling = float(np.nanmax([np.nanmax(target), np.nanmax(control)]))
    reach = float(np.nanmax(np.abs(delta))) or 1.0

    # `with_extremes` rather than `set_bad`: the setter is on its way out, and
    # it mutates a colormap in place, which on a registered one would change
    # it for every other figure drawn in the process.
    ramp = style.sequential_cmap().with_extremes(bad=style.BAND)
    # RdBu, not RdBu_r: positive delta -- the thing you want -- comes out
    # BLUE, the same hue the two magnitude panels use for signal. Red for
    # "higher on the control plate" then reads as the warning it is. The
    # reverse map is the more common default and says the opposite.
    diverging = plt.get_cmap("RdBu").with_extremes(bad=style.BAND)

    truth = pair.truth() if mark_truth and hasattr(pair, "truth") else None
    binders = set(truth[truth["is_binder"]]["well"]) if truth is not None else set()
    cross = (set(truth[truth["is_control_binder"]]["well"])
             if truth is not None else set())

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 3.15), dpi=150)
    panels = [
        (axes[0], target, "Target antigen", ramp, 0, ceiling),
        (axes[1], control, "Control antigen (HLA-His)", ramp, 0, ceiling),
        (axes[2], delta, "Target − control", diverging, -reach, reach),
    ]
    for ax, data, title, cmap, low, high in panels:
        image = ax.imshow(data, vmin=low, vmax=high, cmap=cmap)
        _plate_axes(ax, plate_format)
        ax.set_title(title, fontsize=10, color=style.INK, loc="left", pad=8)
        fig.colorbar(image, ax=ax, shrink=0.78,
                     label="OD 450" if high > 0 and low == 0 else "Δ OD 450")

    # Ring the wells that really are binders, on all three panels.
    for ax, _data, *_rest in panels:
        for well in binders:
            position = wellmod.normalize(well)
            ax.add_patch(plt.Circle(
                (int(position[1:]) - 1, ord(position[0]) - ord("A")), 0.42,
                fill=False, edgecolor=style.STATUS["good"], linewidth=1.4))
        for well in cross - binders:
            position = wellmod.normalize(well)
            ax.add_patch(plt.Circle(
                (int(position[1:]) - 1, ord(position[0]) - ord("A")), 0.42,
                fill=False, edgecolor=style.STATUS["critical"], linewidth=1.4,
                linestyle=":"))

    handles = []
    if np.isnan(target).any():
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=style.BAND,
                                     edgecolor="white",
                                     label=f"no sample ({int(np.isnan(target).sum())})"))
    if binders:
        handles.append(plt.Line2D([], [], marker="o", linestyle="none",
                                  markerfacecolor="none",
                                  markeredgecolor=style.STATUS["good"],
                                  markersize=9,
                                  label=f"binds the target ({len(binders)})"))
    if cross - binders:
        handles.append(plt.Line2D([], [], marker="o", linestyle="none",
                                  markerfacecolor="none",
                                  markeredgecolor=style.STATUS["critical"],
                                  markersize=9,
                                  label=f"binds the control ({len(cross - binders)})"))
    if handles:
        axes[0].legend(handles=handles, loc="upper left",
                       bbox_to_anchor=(0, -0.13), frameon=False, fontsize=8,
                       labelcolor=style.MUTED, ncol=2, handletextpad=0.4,
                       columnspacing=1.2)

    filled = int((pair.target["role"] == "sample").sum()) \
        if "role" in pair.target.columns else len(pair.target)
    fig.suptitle(f"{pair.pair_id} — {pair.scenario.replace('_', ' ')} — "
                 f"{filled} of {plate_format} wells used",
                 fontsize=11, color=style.INK, x=0.008, ha="left", y=1.0)
    fig.text(0.992, 1.0, "SIMULATED — not real measurements", ha="right",
             va="baseline", fontsize=8.5, color=style.MUTED, style="italic")
    fig.subplots_adjust(wspace=0.32)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def clone_journey(chain: list[dict], out: str | Path = "clone_journey.png",
                  title: str | None = None):
    """One sample, through every container it passed through.

    `chain` is what `plates.trace` returns, newest first; this reads it
    oldest first, because that is the direction the liquid went.

    The point of drawing it rather than printing it is that the volumes and
    concentrations change along the way, and a table of eight rows does not
    make "50 µL of DNA became 100 µL at 10 µg/mL" feel like one object.
    """
    out = Path(out)
    steps = list(reversed(chain))
    n = len(steps)
    fig, ax = plt.subplots(figsize=(max(8.0, 1.85 * n), 2.25), dpi=150)
    ax.axis("off")
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-1.0, 1.0)

    for i, hop in enumerate(steps):
        # Colour follows the CONTENT, not the position: a beads plate looks
        # like a beads plate wherever it appears, so the repeat in the middle
        # of this chain is visible as a repeat.
        colour = CONTENT_COLOUR.get(str(hop.get("content_kind")),
                                    style.SEQUENTIAL[3])
        ax.add_patch(plt.Rectangle((i - 0.38, -0.3), 0.76, 0.62,
                                   facecolor=colour, edgecolor="white",
                                   linewidth=1.5))
        ax.text(i, 0.01, str(hop.get("content_kind") or "?"),
                ha="center", va="center", fontsize=8.5,
                color=style.text_on(colour))
        ax.text(i, 0.46, str(hop.get("barcode") or "no barcode"),
                ha="center", va="bottom", fontsize=7.5, color=style.MUTED)
        ax.text(i, 0.72, str(hop.get("well") or ""), ha="center", va="bottom",
                fontsize=8, color=style.INK, fontweight="bold")

        details = []
        if hop.get("volume_ul") is not None:
            details.append(f"{float(hop['volume_ul']):g} µL")
        if hop.get("concentration") is not None:
            details.append(f"{float(hop['concentration']):g} "
                           f"{hop.get('conc_units') or ''}".strip())
        ax.text(i, -0.38, "\n".join(details) or "—", ha="center", va="top",
                fontsize=7.5, color=style.MUTED)

        if i < n - 1:
            ax.annotate("", xy=(i + 0.62, 0.01), xytext=(i + 0.4, 0.01),
                        arrowprops={"arrowstyle": "-|>", "color": style.GRID,
                                    "linewidth": 1.4})

    clone = next((h.get("clone_id") for h in steps if h.get("clone_id")), "?")
    ax.set_title(title or f"{clone} — {n} containers, oldest first",
                 fontsize=11, color=style.INK, loc="left", pad=8)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def pipeline_map(stages: list[dict], out: str | Path = "pipeline.png",
                 title: str = "What this run did"):
    """The whole path, with the number that actually shrank at each step.

    Each stage is {name, count, unit, note}.

    Drawn on a **log axis with visible ticks**, not with log-scaled bar
    widths. Those are the same arithmetic and not the same picture: a bar's
    length is read as proportional to its value, so quietly compressing it
    makes 96 look like two thirds of 616,809. Declaring the scale on the axis
    lets the reader see the six orders of magnitude instead of being walked
    past them, which is the whole point of the figure -- the story here is
    where the numbers fall off.
    """
    out = Path(out)
    counts = [max(float(s.get("count") or 0), 1.0) for s in stages]
    n = len(stages)

    fig, ax = plt.subplots(figsize=(10.5, 0.62 * n + 1.5), dpi=150)
    positions = list(range(n - 1, -1, -1))
    for (stage, count, y) in zip(stages, counts, positions):
        colour = style.SEQUENTIAL[min(2 + (positions.index(y) * 4) // max(n - 1, 1),
                                      len(style.SEQUENTIAL) - 1)]
        ax.barh(y, count, height=0.62, color=colour, zorder=3,
                left=1.0, edgecolor="white", linewidth=1.2)
        label = f"{int(count):,} {stage.get('unit', '')}".strip()
        ax.text(count * 1.3, y, label, ha="left", va="center", fontsize=9,
                color=style.INK)
        if stage.get("note"):
            # To the RIGHT of the count, never inside the bar: faint ink on a
            # saturated fill is unreadable, and that is where it landed first.
            ax.text(count * 1.3, y - 0.34, stage["note"], ha="left", va="top",
                    fontsize=7.5, color=style.FAINT)

    ax.set_xscale("log")
    ax.set_xlim(1, max(counts) * 60)
    ax.set_yticks(positions)
    ax.set_yticklabels([s["name"] for s in stages], fontsize=9.5,
                       color=style.INK)
    ax.set_xlabel("count (log scale)", fontsize=8.5, color=style.MUTED)
    style.style_axes(ax, grid="x")
    ax.tick_params(axis="y", length=0)

    ax.set_title(title, fontsize=12, color=style.INK, loc="left", pad=12)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def hit_calls(pair, callset, out: str | Path = "hit_calls.png",
              plate_format: int = 96):
    """What a caller got right and wrong, on the plate.

    A confusion matrix says 4 false positives. This says *which wells*, which
    is the version you can act on -- a cluster of false positives along one
    edge is a different problem from four scattered ones.

    Status colours, each with a label, because colour alone must never carry
    identity.
    """
    out = Path(out)
    truth = pair.truth()
    merged = callset.calls.merge(truth[["well", "is_binder"]], on="well",
                                 how="left")
    merged["is_binder"] = merged["is_binder"].fillna(False)

    outcome = np.where(merged["is_hit"] & merged["is_binder"], "found",
              np.where(merged["is_hit"] & ~merged["is_binder"], "false call",
              np.where(~merged["is_hit"] & merged["is_binder"], "missed",
                       "correctly ignored")))
    colours = {"found": style.STATUS["good"],
               "false call": style.STATUS["critical"],
               "missed": style.STATUS["warning"],
               "correctly ignored": style.BAND}

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(11.5, 3.3), dpi=150,
        gridspec_kw={"width_ratios": [1.45, 1]})

    n_rows, n_cols = wellmod.dims(plate_format)
    # Every well of the plate, faintly, before anything is drawn on it: a
    # partly-filled plate should look partly filled rather than look like a
    # small plate.
    for r in range(n_rows):
        for c in range(n_cols):
            ax.add_patch(plt.Rectangle((c - 0.45, r - 0.45), 0.9, 0.9,
                                       facecolor="white", edgecolor=style.GRID,
                                       linewidth=0.6, linestyle=":"))
    for well, kind in zip(merged["well"], outcome):
        position = wellmod.normalize(str(well))
        ax.add_patch(plt.Rectangle(
            (int(position[1:]) - 1.45, ord(position[0]) - ord("A") - 0.45),
            0.9, 0.9, facecolor=colours[kind], edgecolor="white",
            linewidth=1.1))
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(n_rows - 0.5, -0.5)
    _plate_axes(ax, plate_format)
    counts = pd.Series(outcome).value_counts()
    ax.set_title(f"{callset.caller} — {pair.scenario.replace('_', ' ')} "
                 f"[{callset.verdict}]", fontsize=10, color=style.INK,
                 loc="left", pad=8)
    ax.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, facecolor=colours[k], edgecolor="white",
                      label=f"{k} ({counts.get(k, 0)})")
        for k in ("found", "missed", "false call", "correctly ignored")
        if counts.get(k, 0)],
        loc="upper left", bbox_to_anchor=(0, -0.12), frameon=False,
        fontsize=8, labelcolor=style.MUTED, ncol=2, handlelength=1.1)

    # The same plate as a scatter: control against target, log-log, with the
    # decision boundary drawn. This is where a threshold becomes visible as a
    # line rather than a number.
    binder = merged["is_binder"].values
    bx.scatter(merged[CONTROL][~binder], merged[TARGET][~binder], s=16,
               color=style.FAINT, edgecolor="white", linewidth=0.5,
               label="not a binder", zorder=3)
    if binder.any():
        bx.scatter(merged[CONTROL][binder], merged[TARGET][binder], s=26,
                   color=style.STATUS["good"], edgecolor="white", linewidth=0.6,
                   label="binds the target", zorder=4)
    low = max(1e-2, float(min(merged[CONTROL].min(), merged[TARGET].min())))
    high = float(max(merged[CONTROL].max(), merged[TARGET].max())) * 1.3
    line = np.array([low, high])
    ratio = callset.params.get("min_ratio")
    if ratio:
        bx.plot(line, line * ratio, color=style.INK, linewidth=1.0,
                linestyle="--", zorder=2, label=f"{ratio:g}× ratio")
    bx.plot(line, line, color=style.GRID, linewidth=1.0, zorder=1)
    bx.set_xscale("log"); bx.set_yscale("log")
    bx.set_xlim(low, high); bx.set_ylim(low, high)
    bx.set_xlabel("control OD 450", fontsize=8.5, color=style.MUTED)
    bx.set_ylabel("target OD 450", fontsize=8.5, color=style.MUTED)
    style.style_axes(bx, grid="both")
    bx.legend(frameon=False, fontsize=7.5, labelcolor=style.MUTED, loc="lower right")

    fig.text(0.992, 1.02, "SIMULATED — not real measurements", ha="right",
             va="baseline", fontsize=8.5, color=style.MUTED, style="italic")
    fig.subplots_adjust(wspace=0.28)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out
