"""SI figures for fixed-reference placement scores and rank occupancy.

The entropy-treatment figure is the four-panel prototype: one panel per
treatment, splitA/splitB pairs joined per family, showing that the selected
six-family order is preserved under every low-frequency treatment.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .figure1 import B3LYP_C, CM, WB97_C, save_figure

QRRHO_C = "#CC79A7"
CANON = "#215C38"; NONCAN = "#727272"
CANONICAL = {"adenine", "cytosine", "guanine"}
TREATMENTS = [
    ("harmonic", "Harmonic", "H", "#4D4D4D"),
    ("floor_50", "50 cm$^{-1}$ floor", "50", B3LYP_C),
    ("floor_100", "100 cm$^{-1}$ floor", "100", WB97_C),
    ("grimme_qrrho_100", "Grimme qRRHO (100 cm$^{-1}$)", "q", QRRHO_C),
]
PANEL_TITLES = {
    "harmonic": "Harmonic",
    "floor_50": "50 cm$^{-1}$ floor",
    "floor_100": "100 cm$^{-1}$ floor",
    "grimme_qrrho_100": "Grimme qRRHO (100 cm$^{-1}$)",
}
PLACEMENT_MARKERS = {"splitA": "o", "splitB": "s"}
NAMES = {"tap": "TAP (C5-linked)"}


def _display_order(config: dict) -> list[str]:
    order = list(config["si_microsolvation"]["plot_family_order"])
    if len(order) != 6 or len(set(order)) != 6:
        raise ValueError("SI display order must contain six unique families")
    return order


def render_si_microsolvation(tables: Path, output: Path, dpi: int, config: dict) -> Path:
    """Render the four-panel entropy-treatment placement figure (SI).

    One panel per low-frequency treatment; splitA/splitB scores are joined
    per family. Single-column PNAS width (8.7 cm); all text >= 6 pt.
    """
    figure = config["si_microsolvation"]
    width, height = float(figure["width_cm"]), float(figure["height_cm"])
    x_limits = tuple(float(x) for x in figure["x_limits_kcal_mol"])
    if dpi < 1 or width <= 0 or height <= 0 or x_limits[0] >= x_limits[1]:
        raise ValueError("SI figure dimensions, resolution and limits must be valid")
    order = _display_order(config)
    placements = pd.read_csv(tables / "si" / "05_fixed_reference_placement_scores.csv")
    expected = {(t, f, p) for t, *_ in TREATMENTS for f in order for p in PLACEMENT_MARKERS}
    actual = set(zip(placements.treatment, placements.family, placements.placement))
    if actual != expected or len(placements) != 48:
        raise ValueError("SI figure requires all 48 treatment/family/placement scores")
    # The six families must keep the selected harmonic order under every treatment.
    for treatment, _title in PANEL_TITLES.items():
        minima = (placements[placements.treatment.eq(treatment)]
                  .groupby("family").fixed_reference_score_kcal_mol.min())
        if minima.sort_values().index.tolist() != list(figure["display_family_order"]):
            raise ValueError(f"Selected order changed under {treatment}")

    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.8})
    output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(width * CM, height * CM))
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=NONCAN, mec=NONCAN, mew=0.7, ms=4.0, label="Placement A"),
        Line2D([], [], marker="s", ls="", mfc=NONCAN, mec=NONCAN, mew=0.7, ms=4.0, label="Placement B"),
    ], loc="upper center", bbox_to_anchor=(0.53, 0.979), ncol=2,
        frameon=False, handletextpad=0.45, columnspacing=1.0, borderaxespad=0, fontsize=6.5)

    rows = []
    bottoms = [0.755, 0.552, 0.349, 0.146]
    axes = []
    for (treatment, title), bottom in zip(PANEL_TITLES.items(), bottoms):
        ax = fig.add_axes((0.30, bottom, 0.57, 0.163))
        axes.append(ax)
        block = placements[placements.treatment.eq(treatment)].set_index(["family", "placement"])
        for yi, family in enumerate(order):
            a = float(block.loc[(family, "splitA"), "fixed_reference_score_kcal_mol"])
            b = float(block.loc[(family, "splitB"), "fixed_reference_score_kcal_mol"])
            color = CANON if family in CANONICAL else NONCAN
            ax.plot((a, b), (yi, yi), color=color, alpha=0.64, linewidth=1.05, solid_capstyle="round", zorder=2)
            if treatment == "grimme_qrrho_100" and family == "tap":
                # Both x coordinates are retained. Their separation is below
                # the printed marker width, so use a nested composite glyph.
                ax.scatter(a, yi, marker="o", s=33, facecolors=color,
                           edgecolors=color, linewidths=0.6, zorder=4)
                ax.scatter(b, yi, marker="s", s=12, facecolors="white",
                           edgecolors=color, linewidths=0.6, zorder=5)
            else:
                ax.scatter(a, yi, marker="o", s=16, facecolors=color,
                           edgecolors=color, linewidths=0.6, zorder=4)
                ax.scatter(b, yi, marker="s", s=16, facecolors=color,
                           edgecolors=color, linewidths=0.6, zorder=5)
            for placement, value in (("splitA", a), ("splitB", b)):
                rows.append({"panel": treatment, "family": family, "treatment": treatment,
                             "placement": placement, "series": "fixed_reference_placement_score",
                             "selected_for_family": bool(
                                 block.loc[(family, placement), "selected_for_family"]),
                             "value": value})
        ax.set_yticks(range(6), [NAMES.get(f, f.capitalize()) for f in order])
        for tick, family in zip(ax.get_yticklabels(), order):
            tick.set_color(CANON if family in CANONICAL else NONCAN)
            tick.set_fontweight("bold" if family in CANONICAL else "normal")
        ax.set_ylim(len(order) - 0.55, -0.45)
        ax.set_xlim(*x_limits)
        ax.set_xticks((0, 3, 6, 9, 12))
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#E0E0E0", linewidth=0.6)
        ax.axvline(0, color="#9E9E9E", linestyle=(0, (3, 2)), linewidth=0.8, zorder=1)
        ax.tick_params(axis="y", length=0, pad=4)
        ax.tick_params(axis="x", length=3, width=0.8, labelsize=6.5,
                       labelbottom=treatment == "grimme_qrrho_100")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=7.5, pad=1.5)

    axes[-1].set_xticklabels(("0\nPCM adenine", "3", "6", "9", "12"))
    axes[-1].get_xticklabels()[0].set_fontweight("bold")
    axes[-1].set_xlabel(
        "Relative Nucleotide " + r"$\Delta G$" + "\n"
        + r"with two waters (kcal mol$^{-1}$)",
        fontsize=7.0, labelpad=8,
    )
    if not np.allclose(fig.get_size_inches() * 2.54, [width, height]):
        raise ValueError("Exported SI figure size differs from configured dimensions")
    save_figure(fig, output / "si_microsolvation_sensitivity", dpi)
    plt.close(fig)
    plotted = output / "si_microsolvation_plotted_data.csv"
    pd.DataFrame(rows).to_csv(plotted, index=False)
    return plotted


