"""Figure 1: continuum rankings, method sensitivity and product-side microsolvation.

A and B plot the continuum relative rankings. C plots
``fixed_reference_score_kcal_mol`` from the primary microsolvation table,
x_i(2W) = r_i(PCM) + G(N_i.2W) - G(N_i) - 2G(W) + C, computed in
``analysis.microsolvation``. Both optimized water arrangements are displayed;
no reactant hydration is added.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
# A GUI backend may quantize the requested size to its display pixels.
# Export headlessly to preserve the exact centimeter dimensions in vector files.
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

CM = 1 / 2.54

CANON = "#245A35"; NONCAN = "#777777"
BAND = "#EDEDED"; GRID = "#E3E3E3"; INK = "#1a1a1a"; SOFT = "#666666"
B3LYP_C = "#0072B2"; WB97_C = "#D55E00"
METHODS = [("1_b3lyp_pcm", "B3LYP / PCM", "o", B3LYP_C, True),
           ("2_b3lyp_smd", "B3LYP / SMD", "o", B3LYP_C, False),
           ("3_wb97xd_pcm", "ωB97X-D / PCM", "^", WB97_C, True),
           ("4_wb97xd_smd", "ωB97X-D / SMD", "^", WB97_C, False)]
CANONICAL = {"adenine", "cytosine", "guanine", "thymine", "uracil"}
FS = 7

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": FS, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": FS + 0.5,
                     "legend.fontsize": FS - 0.5, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "svg.hashsalt": "atcgu-study"})
# Undated vector exports so identical inputs give identical PDF bytes.
SAVE_METADATA = {"png": None, "pdf": {"CreationDate": None}}

def style_labels(ax, fams):
    for t, f in zip(ax.get_yticklabels(), fams):
        c = f in CANONICAL
        t.set_fontweight("bold" if c else "normal"); t.set_color(CANON if c else SOFT)

def refline(ax):
    ax.axvline(0, color="#999999", lw=0.6, ls="--", zorder=1)
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0, pad=3)

XL = r"Relative reaction $\Delta G$ vs adenine (kcal mol$^{-1}$)"

def save_figure(fig, stem: Path, dpi: int) -> None:
    """Write PNG/PDF/SVG with undated, reproducible vector metadata."""
    for ext, metadata in SAVE_METADATA.items():
        fig.savefig(stem.with_name(f"{stem.name}.{ext}"), dpi=dpi if ext == "png" else None, metadata=metadata)


def render_figure1(tables: Path, output: Path, dpi: int | None, config: dict) -> Path:
    """Render Figure 1 from a reproduction's tables."""
    figure = config["figure1"]
    WIDTH_CM, HEIGHT_CM = float(figure["width_cm"]), float(figure["height_cm"])
    dpi = int(figure["png_dpi"]) if dpi is None else dpi
    RXN, REFERENCE = figure["reaction_id"], figure["reference_family"]
    if dpi < 1 or WIDTH_CM <= 0 or HEIGHT_CM <= 0:
        raise ValueError("Figure dimensions and resolution must be positive")
    OUT = output
    OUT.mkdir(parents=True, exist_ok=True)
    SRC_A = tables / 'figure1_panel_a.csv'
    SRC_B = tables / '09_relative_rankings.csv'
    SRC_C = tables / 'micro_16_panel_c_plotted_data.csv'
    SRC_C_CAND = tables / 'micro_05_primary_candidate_energies.csv'
    # ---------- data: A ----------
    a = pd.read_csv(SRC_A).sort_values("rank_within_panel").reset_index(drop=True)
    if not (a.analysis_id.eq(figure["panel_a_analysis_id"]).all() and a.reaction_id.eq(RXN).all()):
        raise ValueError("Panel A table does not match the configured analysis and reaction")
    if len(a) != 21 or a.relative_g_kcal_mol.iloc[0] != 0.0 or a.reaction_family_id.iloc[0] != REFERENCE:
        raise ValueError("Panel A requires 21 families with the reference family at zero")
    names = {**dict(zip(a.reaction_family_id, a.family_display_name)), **figure.get("short_names", {})}
    chk = dict(zip(a.reaction_family_id, a.relative_g_kcal_mol))

    # ---------- data: B ----------
    b = pd.read_csv(SRC_B)
    b = b[(b.family_panel == figure["panel_b_family_panel"]) & (b.reaction_id == RXN)]
    bp = b.pivot_table(index="reaction_family_id", columns="method_id", values="relative_g_kcal_mol")
    if bp.shape != (10, 4):
        raise ValueError(f"Panel B requires 10 families x 4 methods; found {bp.shape}")
    if not all(abs(chk[f] - bp.loc[f, "1_b3lyp_pcm"]) < 1e-9 for f in bp.index):
        raise ValueError("Panel B B3LYP/PCM values disagree with panel A")
    b_order = bp["1_b3lyp_pcm"].sort_values().index.tolist()
    targeted = set(bp.index)

    # ---------- data: C ----------
    c = pd.read_csv(SRC_C).sort_values("continuum_rank").reset_index(drop=True)
    FAMMAP = {"tap": "tap_c5link"}
    c_fam = [FAMMAP.get(f, f) for f in c.family]
    cont = dict(zip(c.family, c.continuum_relative_g_kcal_mol))
    if not all(abs(chk[FAMMAP.get(f, f)] - v) < 1e-6 for f, v in cont.items()):
        raise ValueError("Panel C continuum values disagree with panel A")
    e = pd.read_csv(SRC_C_CAND)
    e = e[(e.treatment == "harmonic") & (e.panel_role == "primary_matched_split_site")].copy()
    if len(e) != 12 or not e.is_usable.all():
        raise ValueError("Panel C requires 12 usable primary candidates")
    # The common reference retains the reference family's own association shift.
    e["x"] = e.fixed_reference_score_kcal_mol
    best = e.loc[e.groupby("family").x.idxmin()].set_index("family")
    other = e[~e.candidate_id.isin(best.candidate_id)].set_index("family")
    if len(other) != 6:
        raise ValueError("Panel C requires one unselected placement per family")
    # The selected points must reproduce the reference-relative values of micro_16.
    for f, v in zip(c.family, c.microsolvated_relative_g_kcal_mol):
        if abs((best.x[f] - best.x[REFERENCE]) - v) >= 1e-6:
            raise ValueError(f"Panel C selected value disagrees with the analysis table: {f}")

    # ---------- figure ----------
    fig = plt.figure(figsize=(WIDTH_CM * CM, HEIGHT_CM * CM))
    # Keep the 21-family column readable while enlarging both right-hand plots.
    # Explicit rectangles prevent title/legend changes from altering export size.
    axA = fig.add_axes([0.205, 0.090, 0.275, 0.785])
    axB = fig.add_axes([0.610, 0.595, 0.375, 0.280])
    axC = fig.add_axes([0.610, 0.090, 0.375, 0.280])

    # Panel A
    y = list(range(len(a)))
    for yi, f in zip(y, a.reaction_family_id):
        if f in targeted: axA.axhspan(yi - 0.42, yi + 0.42, color=BAND, zorder=0, lw=0)
    for yi, r in zip(y, a.itertuples()):
        if r.canonical:
            axA.scatter(r.relative_g_kcal_mol, yi, marker="D", s=16, color=CANON, edgecolor="white", lw=0.4, zorder=3)
        else:
            axA.scatter(r.relative_g_kcal_mol, yi, marker="o", s=14, facecolor="white", edgecolor=NONCAN, lw=0.9, zorder=3)
        axA.text(r.relative_g_kcal_mol + 0.45, yi, f"{r.relative_g_kcal_mol:.2f}", va="center", fontsize=6, color=SOFT)
    axA.set_yticks(y, [names[f] for f in a.reaction_family_id]); style_labels(axA, a.reaction_family_id)
    axA.set_ylim(len(a) - 0.5, -0.6); axA.set_xlim(-0.8, 16.8); axA.set_xticks(range(0, 17, 4))
    refline(axA); axA.set_xlabel(XL)
    axA.set_title("21 reaction families, B3LYP/PCM", loc="left", fontweight="bold", pad=24)
    axA.legend(handles=[Line2D([], [], marker="D", ls="", color=CANON, mec="white", ms=4.5, label="Canonical"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec=NONCAN, mew=0.9, ms=4.2, label="Noncanonical"),
                        Patch(color=BAND, label="Families in B")],
               loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, handlelength=1.1,
               handletextpad=0.4, columnspacing=0.9, borderaxespad=0.2)

    # Panel B
    yb = list(range(len(b_order)))
    for yi, f in zip(yb, b_order):
        vals = bp.loc[f]
        axB.plot([vals.min(), vals.max()], [yi, yi], color=CANON if f in CANONICAL else "#C4C4C4",
                 lw=1.2, alpha=0.65 if f in CANONICAL else 1.0, zorder=1, solid_capstyle="round")
        for mid, lab, mk, col, filled in METHODS:
            axB.scatter(vals[mid], yi, marker=mk, s=16, facecolor=col if filled else "white",
                        edgecolor=col, lw=0.9, zorder=3 if filled else 4)
    axB.set_yticks(yb, [names[f] for f in b_order]); style_labels(axB, b_order)
    axB.set_ylim(len(b_order) - 0.5, -0.6); axB.set_xlim(-1.8, 10.8); axB.set_xticks(range(0, 11, 2))
    refline(axB); axB.set_xlabel(XL)
    axB.set_title("Method sensitivity", loc="left", fontweight="bold", pad=24)
    axB.legend(handles=[Line2D([], [], marker=mk, ls="", mfc=col if filled else "white", mec=col, mew=0.9, ms=4.3, label=lab)
                        for _, lab, mk, col, filled in METHODS],
               loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, handletextpad=0.3,
               columnspacing=0.9, borderaxespad=0.2)

    # Panel C: one row per family. The thin connector joins two separately
    # optimized two-water arrangements, not a statistical interval. The filled
    # circle is the lower sampled G, the open circle is the other arrangement, and
    # the open diamond retains the continuum score.
    yc = list(range(len(c)))
    for yi, fc, f in zip(yc, c.family, c_fam):
        x0, x1, x2 = cont[fc], best.x[fc], other.x[fc]
        col = CANON if f in CANONICAL else NONCAN
        axC.plot([x1, x2], [yi, yi], color=CANON if f in CANONICAL else "#C4C4C4",
                 lw=1.2, alpha=0.65 if f in CANONICAL else 1.0,
                 solid_capstyle="round", zorder=1)
        axC.scatter(x0, yi, marker="D", s=28, facecolor="white", edgecolor=col, lw=0.9, zorder=5)
        axC.scatter(x1, yi, s=40, color=col, edgecolor="white", lw=0.45, zorder=6)
        # Keep the other placement visible where it nearly coincides with the
        # PCM-only diamond (notably for imidazole).
        axC.scatter(x2, yi, marker="o", s=23, facecolor="white", edgecolor=col,
                    linewidths=0.9, zorder=7)
    axC.set_yticks(yc, [names[f] for f in c_fam]); style_labels(axC, c_fam)
    axC.set_ylim(len(c) - 0.5, -0.8); axC.set_xlim(-0.6, 8.2); axC.set_xticks(range(0, 9, 2))
    refline(axC); axC.set_xlabel(r"Relative Nucleotide $\Delta G$ vs PCM adenine (kcal mol$^{-1}$)", fontsize=FS-0.5)
    axC.set_title("Nucleotide (product) local water sensitivity", loc="left", fontweight="bold", pad=47)
    axC.legend(handles=[Line2D([], [], marker="D", ls="", mfc="white", mec="#555555", mew=0.9, ms=4.3, label="PCM only (no explicit water)"),
                        Line2D([], [], marker="o", ls="", color="#555555", mec="white", ms=6.3, label="+ 2 H$_2$O: lowest sampled G"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", mew=0.9, ms=4.8,
                               label="+ 2 H$_2$O: other placement")],
               loc="lower left", bbox_to_anchor=(0.0, 1.075), ncol=1, frameon=False, handletextpad=0.3,
               columnspacing=0.9, borderaxespad=0.2)

    # Panel letters align with the left-aligned titles, just left of each plot.
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    for ax, label in [(axA, "A"), (axB, "B"), (axC, "C")]:
        x = ax.get_position().x0 - 0.035
        y = inv.transform((0, ax._left_title.get_window_extent(rend).y1))[1]
        fig.text(x, y, label, fontsize=10, fontweight="bold", ha="left", va="top")
    # Divide the right-hand panels at the exact midpoint of the visible text gap.
    b_label_bottom = axB.xaxis.label.get_window_extent(rend).y0
    c_title_top = axC._left_title.get_window_extent(rend).y1
    divider_y = inv.transform((0, (b_label_bottom + c_title_top) / 2))[1]
    fig.add_artist(Line2D([axB.get_position().x0 - 0.035, axB.get_position().x1],
                          [divider_y, divider_y], transform=fig.transFigure,
                          color="#B8B8B8", lw=0.55))
    if not np.allclose(fig.get_size_inches() * 2.54, [WIDTH_CM, HEIGHT_CM]):
        raise ValueError("Exported figure size differs from the configured dimensions")
    save_figure(fig, OUT / "figure1", dpi)

    rows = [dict(panel="A", family=r.reaction_family_id, series="B3LYP/PCM", value=r.relative_g_kcal_mol) for r in a.itertuples()]
    rows += [dict(panel="B", family=f, series=m, value=bp.loc[f, m]) for f in b_order for m, *_ in METHODS]
    for fc, f in zip(c.family, c_fam):
        rows += [dict(panel="C", family=f, series="pcm_only", value=cont[fc]),
                 dict(panel="C", family=f, series=f"pcm_plus_2w_lower_{best.placement[fc]}", value=best.x[fc]),
                 dict(panel="C", family=f, series=f"pcm_plus_2w_other_{other.placement[fc]}", value=other.x[fc])]
    pd.DataFrame(rows).to_csv(OUT / "figure1_plotted_data.csv", index=False)
    plt.close(fig)
    return OUT / "figure1_plotted_data.csv"
