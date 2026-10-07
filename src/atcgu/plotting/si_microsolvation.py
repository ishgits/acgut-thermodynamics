"""SI figures for fixed-reference placement scores and rank occupancy."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from .figure1 import B3LYP_C, CM, WB97_C, refline, save_figure, style_labels

QRRHO_C = "#CC79A7"
TREATMENTS = [
    ("harmonic", "Harmonic", "H", "#4D4D4D"),
    ("floor_50", "50 cm$^{-1}$ floor", "50", B3LYP_C),
    ("floor_100", "100 cm$^{-1}$ floor", "100", WB97_C),
    ("grimme_qrrho_100", "Grimme qRRHO (100 cm$^{-1}$)", "q", QRRHO_C),
]
PLACEMENT_MARKERS = {"splitA": "o", "splitB": "s"}
NAMES = {"tap": "TAP (C5-linked)"}


def _display_order(config: dict) -> list[str]:
    order = list(config["si_microsolvation"]["display_family_order"])
    if len(order) != 6 or len(set(order)) != 6:
        raise ValueError("SI display order must contain six unique families")
    return order


def render_si_microsolvation(tables: Path, output: Path, dpi: int, config: dict) -> Path:
    """Render the fixed-reference placement pairs and selected-score shifts."""
    figure = config["si_microsolvation"]
    width, height = float(figure["width_cm"]), float(figure["height_cm"])
    x_limits = tuple(float(x) for x in figure["x_limits_kcal_mol"])
    if dpi < 1 or width <= 0 or height <= 0 or x_limits[0] >= x_limits[1]:
        raise ValueError("SI figure dimensions, resolution and limits must be valid")
    order = _display_order(config)
    si = tables / "si"
    placements = pd.read_csv(si / "05_fixed_reference_placement_scores.csv")
    shifts = pd.read_csv(si / "08_selected_fixed_reference_shifts.csv")
    expected = {(t, f, p) for t, *_ in TREATMENTS for f in order for p in PLACEMENT_MARKERS}
    actual = set(zip(placements.treatment, placements.family, placements.placement))
    if actual != expected or len(placements) != 48:
        raise ValueError("SI panel A requires all 48 treatment/family/placement scores")
    if not placements.groupby(["treatment", "family"]).selected_for_family.sum().eq(1).all():
        raise ValueError("Each SI placement pair must contain exactly one selected marker")
    selected = placements.loc[placements.selected_for_family]
    other = placements.loc[~placements.selected_for_family]
    pair_check = selected.merge(other, on=["treatment", "family"], suffixes=("_selected", "_other"), validate="one_to_one")
    if (pair_check.fixed_reference_score_kcal_mol_selected > pair_check.fixed_reference_score_kcal_mol_other + 1e-10).any():
        raise ValueError("A filled marker would not identify the lower endpoint")

    output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(width * CM, height * CM))
    ax_a = fig.add_axes([0.155, 0.105, 0.525, 0.78])
    ax_b = fig.add_axes([0.775, 0.105, 0.21, 0.78])
    row_gap, group_gap = 0.62, 0.68
    y_lookup, centers = {}, []
    y = 0.0
    for family in order:
        ys = []
        for treatment, _, _, _ in TREATMENTS:
            y_lookup[family, treatment] = y
            ys.append(y)
            y += row_gap
        centers.append(float(np.mean(ys)))
        y += group_gap
    rows = []
    for family in order:
        family_rows = placements.loc[placements.family.eq(family)].set_index(["treatment", "placement"])
        for treatment, _, short, colour in TREATMENTS:
            y = y_lookup[family, treatment]
            a = family_rows.loc[treatment, "splitA"]
            b = family_rows.loc[treatment, "splitB"]
            ax_a.plot([a.fixed_reference_score_kcal_mol, b.fixed_reference_score_kcal_mol], [y, y],
                      color=colour, lw=1.25, alpha=0.75, solid_capstyle="round", zorder=1)
            for placement, point in [("splitA", a), ("splitB", b)]:
                filled = bool(point.selected_for_family)
                ax_a.scatter(point.fixed_reference_score_kcal_mol, y, marker=PLACEMENT_MARKERS[placement], s=27,
                             facecolor=colour if filled else "white", edgecolor=colour, linewidth=0.9, zorder=3)
                rows.append({"panel": "A", "family": family, "treatment": treatment,
                             "placement": placement, "series": "fixed_reference_placement_score",
                             "selected_for_family": filled, "value": point.fixed_reference_score_kcal_mol})
            ax_a.text(x_limits[0] + 0.12, y, short, fontsize=5.7, color=colour, va="center", ha="left")
    ax_a.set_yticks(centers, [NAMES.get(f, f.capitalize()) for f in order])
    style_labels(ax_a, order)
    ax_a.set_ylim(max(y_lookup.values()) + row_gap, -row_gap)
    ax_a.set_xlim(*x_limits)
    refline(ax_a)
    ax_a.set_xlabel(r"Two-water score relative to harmonic continuum adenine (kcal mol$^{-1}$)")
    ax_a.set_title("Fixed-reference placement pairs", loc="left", fontweight="bold", pad=42)
    selection_legend = ax_a.legend(
        handles=[Line2D([], [], marker="o", ls="", mfc="#555555", mec="#555555", ms=4.8,
                        label="Lower-G sampled placement"),
                 Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", ms=4.8,
                        label="Other sampled placement")],
        loc="lower left", bbox_to_anchor=(0.0, 1.055), ncol=2, frameon=False,
        handletextpad=0.35, columnspacing=0.9, borderaxespad=0.0)
    ax_a.add_artist(selection_legend)
    ax_a.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", ms=4.4, label="splitA"),
                         Line2D([], [], marker="s", ls="", mfc="white", mec="#555555", ms=4.2, label="splitB")],
                loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=2, frameon=False,
                handletextpad=0.3, columnspacing=0.8, borderaxespad=0.0)

    nonharmonic = [(t, label, colour) for t, label, _, colour in TREATMENTS if t != "harmonic"]
    shift_plot = shifts.loc[~shifts.treatment.eq("harmonic")].copy()
    if len(shift_plot) != 18:
        raise ValueError("SI panel B requires three non-harmonic shifts for each family")
    offsets = {"floor_50": -0.16, "floor_100": 0.0, "grimme_qrrho_100": 0.16}
    markers = {"floor_50": "o", "floor_100": "^", "grimme_qrrho_100": "s"}
    for yi, family in enumerate(order):
        family_rows = shift_plot.loc[shift_plot.family.eq(family)].set_index("treatment")
        for treatment, _, colour in nonharmonic:
            point = family_rows.loc[treatment]
            py = yi + offsets[treatment]
            value = point.selected_fixed_reference_shift_from_harmonic_kcal_mol
            ax_b.scatter(value, py, marker=markers[treatment], s=21, facecolor=colour, edgecolor=colour,
                         linewidth=0.7, zorder=3)
            if bool(point.placement_switched_from_harmonic):
                ax_b.scatter(value, py, marker="o", s=43, facecolor="none", edgecolor="#222222", linewidth=0.75, zorder=4)
            rows.append({"panel": "B", "family": family, "treatment": treatment,
                         "placement": point.selected_placement, "series": "selected_fixed_reference_shift_from_harmonic",
                         "selected_for_family": True, "value": value})
    ax_b.set_yticks(range(len(order)), [NAMES.get(f, f.capitalize()) for f in order])
    style_labels(ax_b, order)
    ax_b.set_ylim(len(order) - 0.5, -0.5)
    refline(ax_b)
    ax_b.set_xlabel("Selected fixed-reference score shift\n" + r"from harmonic (kcal mol$^{-1}$)", fontsize=6.5)
    ax_b.set_title("Entropy shifts", loc="left", fontweight="bold", pad=42)
    treatment_handles = [Line2D([], [], marker=markers[t], ls="", color=c, ms=4.4, label=label)
                         for t, label, c in nonharmonic]
    treatment_handles.append(Line2D([], [], marker="o", ls="", mfc="none", mec="#222222", ms=6,
                                    label="Placement switch"))
    ax_b.legend(handles=treatment_handles, loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=1,
                frameon=False, handletextpad=0.3, borderaxespad=0.0)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, label in [(ax_a, "A"), (ax_b, "B")]:
        top = inv.transform((0, ax._left_title.get_window_extent(renderer).y1))[1]
        fig.text(ax.get_position().x0 - 0.055, top, label, fontsize=10, fontweight="bold", ha="left", va="top")
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
