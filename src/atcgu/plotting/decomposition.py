"""Decomposition figure: the net reaction split into its two steps.

Two panels rank the 21 families in nucleoside formation
(base + ribose -> nucleoside + H2O) and in phosphorylation
(nucleoside + phosphoric acid -> nucleotide + H2O), each relative to
adenine. A family can look favored in one step and lose it in the other
(2,6-diaminopurine is 2nd at the first step and 21st at the second).
"""
from pathlib import Path
import matplotlib as mpl
# A GUI backend may quantize the requested size to its display pixels.
# Export headlessly to preserve the exact centimeter dimensions in vector files.
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from .figure1 import CANON, NONCAN, SOFT, save_figure, style_labels, refline, CM, FS

PANELS = [
    ("P1a_nucleoside_formation", "Nucleoside formation\nbase + ribose → nucleoside + H₂O"),
    ("P1b_nucleotide_from_nucleoside", "Phosphorylation\nnucleoside + H₃PO₄ → nucleotide + H₂O"),
]


def render_decomposition(tables: Path, output: Path, dpi: int | None, config: dict) -> Path:
    """Render the two-step decomposition figure from a reproduction's tables."""
    dpi = int(config["figure1"]["png_dpi"]) if dpi is None else dpi
    if dpi < 1:
        raise ValueError("PNG resolution must be positive")
    output.mkdir(parents=True, exist_ok=True)
    r09 = pd.read_csv(tables / "09_relative_rankings.csv")
    p = r09[r09["panel"] == "all_21"]
    frames = []
    for rxn_id, _ in PANELS:
        t = (p[p["reaction_id"] == rxn_id].drop_duplicates("family_display_name")
             .sort_values("relative_g_kcal_mol").reset_index(drop=True))
        if len(t) != 21:
            raise ValueError(f"Decomposition panel {rxn_id} requires 21 families, found {len(t)}")
        frames.append(t)

    fig, axes = plt.subplots(1, 2, figsize=(19 * CM, 10 * CM), sharex=False)
    xmax = max(t["relative_g_kcal_mol"].max() for t in frames)
    for ax, (rxn_id, title), t in zip(axes, PANELS, frames):
        fams = t["family_display_name"].tolist()
        vals = t["relative_g_kcal_mol"].tolist()
        colors = [CANON if str(c) == "True" else NONCAN for c in t["canonical"]]
        y = range(len(fams))
        ax.barh(list(y), vals, height=0.62, color=colors, zorder=2)
        ax.set_yticks(list(y)); ax.set_yticklabels(fams, fontsize=FS - 0.5)
        style_labels(ax, t["reaction_family_id"].tolist())
        ax.set_title(title, fontsize=FS, pad=6)
        ax.set_xlim(right=xmax * 1.02)
        refline(ax)
        ax.set_xlabel(r"Relative $\Delta G$ vs adenine (kcal mol$^{-1}$)", fontsize=FS)
        ax.invert_yaxis()
    fig.tight_layout(pad=0.8, w_pad=3.0)
    stem = output / "decomposition"
    save_figure(fig, stem, dpi)
    plt.close(fig)
    return stem.with_suffix(".png")
