"""One visual language for every figure this project makes.

`library` and `assay` both draw, and a palette copied into two modules drifts
within a month. This is the single source, in `core` because both may import
it and neither may import the other (rule 1).

Every categorical set here has been through a CVD validator rather than
eyeballed. The numbers below are the recorded result, so a future edit can be
re-checked against them instead of re-litigated:

    SERIES   (4 slots)  worst adjacent pair ΔE 9.1 protan, 22.9 normal vision
    STATUS   (3 slots)  worst adjacent pair ΔE 9.1 protan, 22.9 normal vision

Both carry a contrast warning against the surface, which is discharged the
way the guidance requires: every figure using them also carries direct
labels, so identity is never colour alone.

`SEQUENTIAL` is a single-hue ramp, light to dark. A ramp is checked for
lightness monotonicity rather than pairwise separation -- adjacent steps are
*supposed* to be close, which is what makes it read as a magnitude.
"""
from __future__ import annotations

# Categorical identity. Assigned in fixed order, never cycled: a fifth series
# folds into "other" or becomes a small multiple rather than inventing a hue.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

# Reserved for state, never reused as "series 5". Always shipped with a label.
STATUS = {
    "good": "#1baf7a",        # a call that was right
    "warning": "#eda100",     # a hit that was missed
    "critical": "#d1495b",    # a call that was wrong
    "neutral": "#9c9a92",     # nothing to report
}

# Single hue, light to dark. Lightness is monotone by construction.
SEQUENTIAL = ["#f2f7fd", "#cfe2f7", "#9cc4ee", "#5f9ee2", "#2a78d6",
              "#1b55a0", "#123a6e"]

# Ink and furniture. Text wears these, never a series colour.
INK = "#0b0b0b"
MUTED = "#52514e"
FAINT = "#9c9a92"
GRID = "#d8d7d2"
SURFACE = "#ffffff"
BAND = "#eceae3"
BAND_HIGHLIGHT = "#ded9f0"


def sequential_cmap(name: str = "plateforge_seq"):
    """The ramp as a matplotlib colormap."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name, SEQUENTIAL)


def style_axes(ax, grid: str = "y") -> None:
    """Recessive furniture: the data should be the loudest thing present."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9, length=2)
    if grid in ("y", "both"):
        ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    if grid in ("x", "both"):
        ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def text_on(hex_color: str) -> str:
    """Ink that stays legible on a given fill."""
    value = hex_color.lstrip("#")
    r, g, b = (int(value[i:i + 2], 16) / 255 for i in (0, 2, 4))
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#111111" if luminance > 0.55 else "#ffffff"


def simulated_banner(fig, text: str = "SIMULATED — not real measurements",
                     y: float = 1.0) -> None:
    """Stamp a figure drawn from synthetic data.

    CLAUDE.md calls presenting synthetic data as real the most damaging
    mistake this repo can make. Every figure built on simulation says so on
    its face, not only in its manifest.
    """
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=9,
             color=MUTED, style="italic")
