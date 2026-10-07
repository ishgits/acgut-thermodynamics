"""Shared numerical conventions and structural cut-offs; one definition each."""
from __future__ import annotations

# CODATA hartree-to-kcal/mol conversion used by every table.
HARTREE_TO_KCAL_MOL = 627.5094740631
# Display convention of the Figure 1C water association shift
# (config/figures.yaml panel_c.hartree_to_kcal_mol); not used for tables.
FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL = 627.509474

# Harmonic thermochemistry temperature and the accepted deviation of a record from it.
TEMPERATURE_K = 298.0
TEMPERATURE_TOLERANCE_K = 0.25

SYMBOLS = {1: "H", 6: "C", 7: "N", 8: "O", 15: "P"}
# Covalent radii (angstrom) by atomic number; bonded when d <= BOND_SCALE * sum.
COVALENT_RADII = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 15: 1.07}
BOND_SCALE = 1.25
MIN_BOND_DISTANCE = 0.1
# Van der Waals radii (angstrom) used for water-placement clearance.
VDW_RADII = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 15: 1.80}
# Hydrogen bond: H...A, D...A (angstrom) and D-H...A angle (degrees).
HB_HA, HB_DA, HB_ANGLE = 2.20, 3.20, 130.0
# Near contact used for the near-fold flag.
NEAR_HA, NEAR_DA, NEAR_ANGLE = 2.50, 3.50, 120.0


def check_study_conventions(study: dict, figures: dict) -> None:
    """Fail when config/study.json or config/figures.yaml drift from these constants."""
    if study["hartree_to_kcal_mol"] != HARTREE_TO_KCAL_MOL:
        raise ValueError("config/study.json hartree_to_kcal_mol disagrees with atcgu.constants")
    for key, value in study["panel_c_display_convention"].items():
        if figures["panel_c"][key] != value:
            raise ValueError("Study/figure numerical conventions disagree")
    if figures["panel_c"]["hartree_to_kcal_mol"] != FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL:
        raise ValueError("config/figures.yaml panel_c conversion disagrees with atcgu.constants")
    if figures["si_microsolvation"]["display_family_order"] != study["si_microsolvation"]["display_family_order"]:
        raise ValueError("Study/SI figure display orders disagree")
