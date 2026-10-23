"""Manuscript Figure S1, adapted from the recovered final standalone source."""
from __future__ import annotations

import pandas as pd
from .figure1 import save_figure
import csv
import math
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


REACTIONS = (
    "P1a_nucleoside_formation",
    "P1b_nucleotide_from_nucleoside",
    "P2b_nucleotide_from_ribophosphate",
    "P1_total_net_nucleotide_formation",
)
CANONICAL = {"adenine", "cytosine", "guanine", "thymine", "uracil"}
LABEL_OVERRIDES = {
    "5_amino_4_imidazolecarboxamide": "AICA",
    "tap_c5link": "TAP (C5-linked)",
    "barbituric_acid_clink": "BA (C5-linked)",
    "barbituric_acid_nlink": "BA (N-linked)",
}
GREEN = "#245A35"
GRAY = "#777777"
SOFT = "#666666"
GRID = "#E3E3E3"
CM = 1 / 2.54
WIDTH_CM = 17.8
HEIGHT_CM = 15.6


def load_rankings(source: Path) -> dict[str, dict[str, dict[str, object]]]:
    """Read the baseline all-family rankings and verify the path identities."""
    panels: dict[str, dict[str, dict[str, object]]] = {key: {} for key in REACTIONS}
    with source.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if not (
                row["analysis_id"] == "full21_b3lyp_pcm"
                and row["method_id"] == "1_b3lyp_pcm"
                and row["family_panel"] == "all_21"
                and row["reaction_id"] in panels
            ):
                continue
            reaction, family = row["reaction_id"], row["reaction_family_id"]
            if family in panels[reaction]:
                raise ValueError(f"Duplicate baseline ranking: {reaction}, {family}")
            value = float(row["relative_g_kcal_mol"])
            if not math.isfinite(value):
                raise ValueError(f"Nonfinite ranking: {reaction}, {family}")
            panels[reaction][family] = {
                "value": value,
                "name": row["family_display_name"],
                "canonical": row["canonical"].lower() == "true",
                "rank": int(row["rank_within_panel"]),
            }

    families = set(panels[REACTIONS[0]])
    if len(families) != 21 or any(set(panel) != families for panel in panels.values()):
        raise ValueError("Expected the same 21 families in all four baseline reactions")
    for reaction, panel in panels.items():
        if abs(panel["adenine"]["value"]) > 1e-9:
            raise ValueError(f"Adenine is not zero in {reaction}")
        if sorted(record["rank"] for record in panel.values()) != list(range(1, 22)):
            raise ValueError(f"Unexpected ranks in {reaction}")
    for family in families:
        a, b, c, net = (panels[key][family]["value"] for key in REACTIONS)
        if abs(a + b - c) > 1e-7 or abs(c - net) > 1e-7:
            raise ValueError(f"Path closure or net comparison failed for {family}")
        if any(panel[family]["canonical"] != (family in CANONICAL) for panel in panels.values()):
            raise ValueError(f"Canonical status differs from Figure 1 for {family}")
    return panels


def render_decomposition(tables: Path, output: Path, dpi: int | None, config: dict) -> Path:
    """Render manuscript Figure S1; return its compatible PNG path."""
    source = tables / "09_relative_rankings.csv"
    dpi = int(config["figure1"]["png_dpi"]) if dpi is None else dpi
    if dpi < 1:
        raise ValueError("PNG resolution must be positive")
    settings = config["figure_s1"]
    panels = load_rankings(source)
    order = sorted(
        panels["P1_total_net_nucleotide_formation"],
        key=lambda family: panels["P1_total_net_nucleotide_formation"][family]["rank"],
    )
    names = [LABEL_OVERRIDES.get(f, panels[REACTIONS[0]][f]["name"]) for f in order]

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(settings["width_cm"] * CM, settings["height_cm"] * CM), facecolor="white")
    # Fit the long family names close to the left page edge, giving each
    # quantitative panel the space recovered from the former outer margin.
    positions = ((0.220, 0.125, 0.225, 0.690),
                 (0.477, 0.125, 0.225, 0.690),
                 (0.734, 0.125, 0.225, 0.690))
    axes = [fig.add_axes(pos) for pos in positions]
    titles = (
        "Nucleoside\nformation",
        "Nucleoside\nphosphorylation",
        "Base addition to\nribose 5-phosphate",
    )
    limits = settings["x_limits_kcal_mol"]
    ticks = ((0, 4, 8, 12), (-2, 0, 2, 4, 6, 8), (0, 4, 8, 12, 16))

    for ax, reaction, title, xlim, xticks, letter in zip(
        axes, REACTIONS[:3], titles, limits, ticks, "ABC"
    ):
        for y, family in enumerate(order):
            value = panels[reaction][family]["value"]
            if family in CANONICAL:
                ax.scatter(value, y, marker="D", s=19, color=GREEN,
                           edgecolor="white", linewidth=0.4, zorder=3)
            else:
                ax.scatter(value, y, marker="o", s=17, facecolor="white",
                           edgecolor=GRAY, linewidth=0.9, zorder=3)

        ax.set_xlim(*xlim)
        ax.set_xticks(xticks)
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.set_yticks(range(len(order)))
        ax.axvline(0, color="#999999", linewidth=0.6, linestyle="--", zorder=1)
        ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=3)
        ax.tick_params(axis="x", labelsize=6.5)
        fig.text(ax.get_position().x0 - 0.022, 0.946, letter, fontsize=10,
                 fontweight="bold", ha="left", va="top")
        fig.text(ax.get_position().x0 + 0.007, 0.946, title, fontsize=7.4,
                 fontweight="bold", ha="left", va="top", linespacing=1.1)

    axes[0].set_yticklabels(names)
    for label, family in zip(axes[0].get_yticklabels(), order):
        label.set_color(GREEN if family in CANONICAL else SOFT)
        label.set_fontweight("bold" if family in CANONICAL else "normal")
    for ax in axes[1:]:
        ax.set_yticklabels([])

    legend = [
        Line2D([], [], marker="D", linestyle="", color=GREEN,
               markeredgecolor="white", markeredgewidth=0.4,
               markersize=4.5, label="Canonical"),
        Line2D([], [], marker="o", linestyle="", markerfacecolor="white",
               markeredgecolor=GRAY, markeredgewidth=0.9,
               markersize=4.2, label="Noncanonical"),
    ]
    fig.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.214, 0.870),
               ncol=2, frameon=False, handletextpad=0.3, columnspacing=1.0,
               borderaxespad=0, fontsize=6.8)
    fig.text(0.590, 0.070,
             r"Relative reaction $\Delta G$ vs adenine (kcal mol$^{-1}$)",
             fontsize=7, ha="center", va="center")

    output.mkdir(parents=True, exist_ok=True)
    stem = output / "decomposition"
    save_figure(fig, stem, dpi)
    plt.close(fig)
    points = [{"panel": letter, "analysis_id": "full21_b3lyp_pcm", "reaction_id": reaction,
               "reaction_family_id": family, "display_order": i + 1,
               "canonical": panels[reaction][family]["canonical"],
               "value_kcal_mol": panels[reaction][family]["value"]}
              for letter, reaction in zip("ABC", REACTIONS[:3]) for i, family in enumerate(order)]
    pd.DataFrame(points).to_csv(output / "figure_s1_plotted_data.csv", index=False)
    return stem.with_suffix(".png")
