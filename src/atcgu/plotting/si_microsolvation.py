"""SI figures for fixed-reference placement scores and rank occupancy.

The entropy-treatment figure is the four-panel prototype: one panel per
treatment, splitA/splitB pairs joined per family, showing that the selected
six-family order is preserved under every low-frequency treatment.
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from .figure1 import B3LYP_C, CM, WB97_C, refline, save_figure, style_labels

QRRHO_C = "#CC79A7"
CANON = "#245A35"; NONCAN = "#777777"
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
    order = list(config["si_microsolvation"]["display_family_order"])
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
        if minima.sort_values().index.tolist() != order:
            raise ValueError(f"Selected order changed under {treatment}")

    output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(width * CM, height * CM))
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=NONCAN, mec=NONCAN, ms=4.2, label="splitA"),
        Line2D([], [], marker="s", ls="", mfc=NONCAN, mec=NONCAN, ms=4.0, label="splitB"),
    ], loc="upper center", bbox_to_anchor=(0.5, 0.985), ncol=2,
        frameon=False, handletextpad=0.3, columnspacing=0.9, fontsize=6.5)

    rows = []
    bottoms = [0.755, 0.552, 0.349, 0.146]
    axes = []
    for (treatment, title), bottom in zip(PANEL_TITLES.items(), bottoms):
        ax = fig.add_axes([0.30, bottom, 0.57, 0.163], sharex=axes[0] if axes else None)
        axes.append(ax)
        block = placements[placements.treatment.eq(treatment)].set_index(["family", "placement"])
        for yi, family in enumerate(order):
            a = float(block.loc[(family, "splitA"), "fixed_reference_score_kcal_mol"])
            b = float(block.loc[(family, "splitB"), "fixed_reference_score_kcal_mol"])
            color = CANON if family in CANONICAL else NONCAN
            ax.plot([a, b], [yi, yi], color=color, lw=0.85, alpha=0.65, zorder=2)
            if treatment == "grimme_qrrho_100" and family == "tap":
                # Both x coordinates are retained. Their separation is below
                # the printed marker width, so use a nested composite glyph.
                ax.scatter(a, yi, marker="o", s=52, color=color,
                           edgecolor="white", linewidth=0.35, zorder=4)
                ax.scatter(b, yi, marker="s", s=19, facecolor="white",
                           edgecolor=color, linewidth=0.8, zorder=5)
            else:
                ax.scatter(a, yi, marker="o", s=25, color=color,
                           edgecolor="white", linewidth=0.35, zorder=4)
                ax.scatter(b, yi, marker="s", s=25, color=color,
                           edgecolor="white", linewidth=0.35, zorder=4)
            for placement, value in (("splitA", a), ("splitB", b)):
                rows.append({"panel": treatment, "family": family, "treatment": treatment,
                             "placement": placement, "series": "fixed_reference_placement_score",
                             "selected_for_family": bool(
                                 block.loc[(family, placement), "selected_for_family"]),
                             "value": value})
        ax.set_yticks(range(6), [NAMES.get(f, f.capitalize()) for f in order])
        style_labels(ax, order)
        ax.tick_params(axis="y", labelsize=6.5)
        ax.set_ylim(5.45, -0.45)
        ax.set_xlim(*x_limits)
        ax.set_xticks(np.arange(0, 13, 3))
        refline(ax)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=7.5, pad=3, color="#1a1a1a")
        if treatment != "grimme_qrrho_100":
            ax.tick_params(labelbottom=False)
        else:
            ax.tick_params(axis="x", labelsize=6.5)

    fig.text(0.585, 0.083,
             "Two-water score relative to harmonic\n"
             r"continuum adenine (kcal mol$^{-1}$)",
             ha="center", va="center", fontsize=6.6)
    if not np.allclose(fig.get_size_inches() * 2.54, [width, height]):
        raise ValueError("Exported SI figure size differs from configured dimensions")
    save_figure(fig, output / "si_microsolvation_sensitivity", dpi)
    plt.close(fig)
    plotted = output / "si_microsolvation_plotted_data.csv"
    pd.DataFrame(rows).to_csv(plotted, index=False)
    return plotted


def render_si_ranking_occupancy(tables: Path, output: Path, dpi: int, config: dict) -> Path:
    """Render counts over the 64 enumerated placement choices for each treatment."""
    figure = config["si_microsolvation_ranking"]
    width, height = float(figure["width_cm"]), float(figure["height_cm"])
    order = _display_order(config)
    occupancy = pd.read_csv(tables / "si/10_rank_occupancy.csv")
    if len(occupancy) != 144 or not occupancy.groupby(["treatment", "family"]).combination_count.sum().eq(64).all():
        raise ValueError("Ranking figure requires 64 rank assignments for every family/treatment")
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(width * CM, height * CM), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.14, right=0.89, bottom=0.12, top=0.91, wspace=0.18, hspace=0.28)
    vmax = int(occupancy.combination_count.max())
    image = None
    for ax, (treatment, label, _, _) in zip(axes.flat, TREATMENTS):
        block = occupancy.loc[occupancy.treatment.eq(treatment)]
        matrix = block.pivot(index="family", columns="rank", values="combination_count").loc[order, range(1, 7)]
        image = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
        for row_i, family in enumerate(order):
            selected_rank = int(block.loc[block.family.eq(family), "selected_lower_g_rank"].iloc[0])
            ax.add_patch(Rectangle((selected_rank - 1.5, row_i - 0.5), 1, 1, fill=False,
                                   edgecolor="#111111", linewidth=1.4))
            for col_i, count in enumerate(matrix.loc[family]):
                colour = "white" if count > vmax * 0.55 else "#1a1a1a"
                ax.text(col_i, row_i, str(int(count)), ha="center", va="center", fontsize=6.3, color=colour)
        ax.set_title(label, loc="left", fontweight="bold")
        ax.set_xticks(range(6), range(1, 7))
        ax.set_yticks(range(6), [NAMES.get(f, f.capitalize()) for f in order])
        style_labels(ax, order)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Rank")
    fig.text(0.015, 0.53, "Family", rotation=90, va="center", fontsize=7)
    cax = fig.add_axes([0.915, 0.18, 0.018, 0.62])
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_label("Count among 64 placement choices")
    fig.legend(handles=[Patch(facecolor="none", edgecolor="#111111", linewidth=1.4,
                              label="Rank from lower-G placement for every family")],
               loc="upper center", bbox_to_anchor=(0.52, 0.985), frameon=False)
    if not np.allclose(fig.get_size_inches() * 2.54, [width, height]):
        raise ValueError("Exported ranking figure size differs from configured dimensions")
    save_figure(fig, output / "si_microsolvation_ranking_occupancy", dpi)
    plt.close(fig)
    plotted = output / "si_microsolvation_ranking_plotted_data.csv"
    occupancy.to_csv(plotted, index=False)
    return plotted
