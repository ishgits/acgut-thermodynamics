"""Review prototype: narrow, four-panel SI entropy sensitivity figure.

Run with ``.venv/bin/python workflows/prototype_si_microsolvation.py``.
The published SI figure and renderer are left untouched during review.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from atcgu.plotting.figure1 import CANON, NONCAN, CM, refline, save_figure, style_labels


ORDER = ["cytosine", "adenine", "imidazole", "melamine", "tap", "guanine"]
LABELS = {"tap": "TAP (C5-linked)"}
CANONICAL = {"cytosine", "adenine", "guanine"}
TREATMENTS = [
    ("harmonic", "Harmonic"),
    ("floor_50", r"50 cm$^{-1}$ floor"),
    ("floor_100", r"100 cm$^{-1}$ floor"),
    ("grimme_qrrho_100", r"Grimme qRRHO (100 cm$^{-1}$)"),
]


def main() -> None:
    source = ROOT / "data/published/si_microsolvation/05_fixed_reference_placement_scores.csv"
    out = ROOT / "results/supplement"
    data = pd.read_csv(source)
    expected = {(t, f, p) for t, _ in TREATMENTS for f in ORDER for p in ("splitA", "splitB")}
    actual = set(zip(data.treatment, data.family, data.placement))
    if actual != expected or len(data) != 48:
        raise ValueError("Prototype requires all 48 published placement scores")

    # These six rows follow the selected harmonic Figure 1C ranking.
    for treatment, _ in TREATMENTS:
        block = data[data.treatment.eq(treatment)]
        minima = block.groupby("family").fixed_reference_score_kcal_mol.min()
        if minima.sort_values().index.tolist() != ORDER:
            raise ValueError(f"Selected order changed under {treatment}")

    # One-column physical width from the NAS PNAS research-article template.
    # SI is delivered in an appendix PDF, so this is a design target rather
    # than a claim that SI figures must use manuscript-column dimensions.
    fig = plt.figure(figsize=(8.7 * CM, 17.8 * CM))
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=NONCAN, mec=NONCAN, ms=4.2, label="splitA"),
        Line2D([], [], marker="s", ls="", mfc=NONCAN, mec=NONCAN, ms=4.0, label="splitB"),
    ], loc="upper center", bbox_to_anchor=(0.5, 0.985), ncol=2,
        frameon=False, handletextpad=0.3, columnspacing=0.9, fontsize=6.5)

    bottoms = [0.755, 0.552, 0.349, 0.146]
    axes = []
    for (treatment, title), bottom in zip(TREATMENTS, bottoms):
        ax = fig.add_axes([0.30, bottom, 0.57, 0.163], sharex=axes[0] if axes else None)
        axes.append(ax)
        block = data[data.treatment.eq(treatment)].set_index(["family", "placement"])
        for yi, family in enumerate(ORDER):
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
        ax.set_yticks(range(6), [LABELS.get(f, f.capitalize()) for f in ORDER])
        style_labels(ax, ORDER)
        ax.tick_params(axis="y", labelsize=6.5)
        ax.set_ylim(5.45, -0.45)
        ax.set_xlim(-0.6, 13.4)
        ax.set_xticks(np.arange(0, 13, 3))
        refline(ax)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=7.5, pad=3,
                     color="#1a1a1a")
        if treatment != "grimme_qrrho_100":
            ax.tick_params(labelbottom=False)
        else:
            ax.tick_params(axis="x", labelsize=6.5)

    fig.text(0.585, 0.083,
             "Two-water score relative to harmonic\n"
             r"continuum adenine (kcal mol$^{-1}$)",
             ha="center", va="center", fontsize=6.6)
    out.mkdir(parents=True, exist_ok=True)
    save_figure(fig, out / "si_microsolvation_prototype", 1200)
    plt.close(fig)


if __name__ == "__main__":
    main()
