"""Balanced reaction construction and fixed-reference relative rankings."""

from __future__ import annotations

import json
import math

import pandas as pd

from ..constants import HARTREE_TO_KCAL_MOL
from ..geometry import formula_counts

REACTION_ORDER = (
    "P1a_nucleoside_formation",
    "P1b_nucleotide_from_nucleoside",
    "P1_total_net_nucleotide_formation",
    "P2a_ribose_phosphate_formation",
    "P2b_nucleotide_from_ribophosphate",
    "P2_total_net_nucleotide_formation",
)

REACTION_METADATA = {
    "P1a_nucleoside_formation": ("Path 1", 1),
    "P1b_nucleotide_from_nucleoside": ("Path 1", 2),
    "P1_total_net_nucleotide_formation": ("Path 1", 3),
    "P2a_ribose_phosphate_formation": ("Path 2", 1),
    "P2b_nucleotide_from_ribophosphate": ("Path 2", 2),
    "P2_total_net_nucleotide_formation": ("Path 2", 3),
}


def reaction_definitions(
    base: str,
    nucleoside: str,
    nucleotide: str,
) -> dict[str, tuple[dict[str, int], dict[str, int]]]:
    return {
        "P1a_nucleoside_formation": ({base: 1, "ribose": 1}, {nucleoside: 1, "water": 1}),
        "P1b_nucleotide_from_nucleoside": ({nucleoside: 1, "phosphate": 1}, {nucleotide: 1, "water": 1}),
        "P1_total_net_nucleotide_formation": (
            {base: 1, "ribose": 1, "phosphate": 1},
            {nucleotide: 1, "water": 2},
        ),
        "P2a_ribose_phosphate_formation": (
            {"ribose": 1, "phosphate": 1},
            {"ribo_phosphate": 1, "water": 1},
        ),
        "P2b_nucleotide_from_ribophosphate": (
            {"ribo_phosphate": 1, base: 1},
            {nucleotide: 1, "water": 1},
        ),
        "P2_total_net_nucleotide_formation": (
            {base: 1, "ribose": 1, "phosphate": 1},
            {nucleotide: 1, "water": 2},
        ),
    }


def required_species_for_reactions(
    families: pd.DataFrame,
    reaction_ids: list[str] | tuple[str, ...],
) -> set[str]:
    required: set[str] = set()
    for family in families.to_dict("records"):
        definitions = reaction_definitions(
            family["base_species_id"],
            family["nucleoside_species_id"],
            family["nucleotide_species_id"],
        )
        for reaction_id in reaction_ids:
            reactants, products = definitions[reaction_id]
            required.update(reactants)
            required.update(products)
    return required


def _residual(
    reactants: dict[str, int],
    products: dict[str, int],
    formulas: dict[str, str],
) -> dict[str, int]:
    result: dict[str, int] = {}
    for sign, side in ((-1, reactants), (1, products)):
        for species, coefficient in side.items():
            counts = formula_counts(formulas[species])
            if not counts:
                raise ValueError(f"Unreadable formula for {species}: {formulas[species]!r}")
            for element, count in counts.items():
                result[element] = result.get(element, 0) + sign * coefficient * count
    return {key: value for key, value in sorted(result.items()) if value}


def build_reactions(
    selected: pd.DataFrame,
    registry: pd.DataFrame,
    families: pd.DataFrame,
    reaction_ids: list[str] | tuple[str, ...] = REACTION_ORDER,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate requested balanced reactions and any applicable path-closure audit."""
    unknown = sorted(set(reaction_ids) - set(REACTION_ORDER))
    if unknown:
        raise ValueError(f"Unknown reaction ids: {unknown}")
    if selected["species_id"].duplicated().any():
        duplicates = sorted(selected.loc[selected["species_id"].duplicated(False), "species_id"].unique())
        raise ValueError(f"Selected-energy table contains duplicate species: {duplicates}")

    energies = selected.set_index("species_id")["computed_gibbs_hartree"].to_dict()
    if not all(math.isfinite(float(v)) for v in energies.values()):
        raise ValueError("Selected-energy table contains nonfinite values")
    charges = registry.set_index("species_id")["charge"].astype(int).to_dict()
    formulas = registry.set_index("species_id")["formula"].to_dict()
    names = registry.set_index("species_id")["display_name"].to_dict()
    records: list[dict[str, object]] = []
    for family in families.sort_values("display_order").to_dict("records"):
        definitions = reaction_definitions(
            family["base_species_id"],
            family["nucleoside_species_id"],
            family["nucleotide_species_id"],
        )
        for reaction_id in REACTION_ORDER:
            if reaction_id not in reaction_ids:
                continue
            path_name, step = REACTION_METADATA[reaction_id]
            reactants, products = definitions[reaction_id]
            required = set(reactants) | set(products)
            missing = sorted(required - set(energies))
            balance = _residual(reactants, products, formulas)
            charge_residual = sum(coef * charges[s] for s, coef in products.items()) - sum(coef * charges[s] for s, coef in reactants.items())
            delta = None if missing else (
                sum(value * energies[key] for key, value in products.items())
                - sum(value * energies[key] for key, value in reactants.items())
            ) * HARTREE_TO_KCAL_MOL

            def format_side(side: dict[str, int]) -> str:
                return " + ".join(
                    ("" if coefficient == 1 else f"{coefficient} ") + names[species]
                    for species, coefficient in side.items()
                )

            records.append(
                {
                    "display_order": int(family["display_order"]),
                    "reaction_family_id": family["reaction_family_id"],
                    "family_display_name": family["display_name"],
                    "canonical": str(family["canonical"]).strip().casefold() == "true",
                    "reaction_id": reaction_id,
                    "path": path_name,
                    "step_order": step,
                    "reaction": f"{format_side(reactants)} -> {format_side(products)}",
                    "reactants": json.dumps(reactants, sort_keys=True),
                    "products": json.dumps(products, sort_keys=True),
                    "missing_species": "; ".join(missing),
                    "element_balance_residual": json.dumps(balance, sort_keys=True),
                    "charge_balance_residual": charge_residual,
                    "balanced": not balance and charge_residual == 0,
                    "delta_g_kcal_mol": delta,
                }
            )
    reaction_table = pd.DataFrame(records)

    closure_columns = [
        "display_order",
        "reaction_family_id",
        "display_name",
        "path1_step_sum_error_kcal_mol",
        "path2_step_sum_error_kcal_mol",
        "net_path_difference_kcal_mol",
        "passes",
    ]
    if not set(REACTION_ORDER).issubset(reaction_ids):
        return reaction_table, pd.DataFrame(columns=closure_columns)

    wide = reaction_table.pivot(index="reaction_family_id", columns="reaction_id", values="delta_g_kcal_mol")
    closure = families[["display_order", "reaction_family_id", "display_name"]].merge(
        wide.reset_index(), on="reaction_family_id"
    )
    closure["path1_step_sum_error_kcal_mol"] = (
        closure["P1a_nucleoside_formation"]
        + closure["P1b_nucleotide_from_nucleoside"]
        - closure["P1_total_net_nucleotide_formation"]
    )
    closure["path2_step_sum_error_kcal_mol"] = (
        closure["P2a_ribose_phosphate_formation"]
        + closure["P2b_nucleotide_from_ribophosphate"]
        - closure["P2_total_net_nucleotide_formation"]
    )
    closure["net_path_difference_kcal_mol"] = (
        closure["P1_total_net_nucleotide_formation"]
        - closure["P2_total_net_nucleotide_formation"]
    )
    error_columns = [
        "path1_step_sum_error_kcal_mol",
        "path2_step_sum_error_kcal_mol",
        "net_path_difference_kcal_mol",
    ]
    closure["passes"] = closure[error_columns].abs().le(1e-7).all(axis=1)
    return reaction_table, closure[closure_columns]


def add_relative_rankings(
    reactions: pd.DataFrame,
    *,
    panel_name: str,
    reference_family: str,
) -> pd.DataFrame:
    """Reference each reaction to ``reference_family`` and rank within the supplied family panel.

    Ranks use ``rank(method="min")``: tied energies share the lower rank.
    """
    parts: list[pd.DataFrame] = []
    for reaction_id, group in reactions.groupby("reaction_id", sort=False):
        reference = group.loc[group["reaction_family_id"].eq(reference_family), "delta_g_kcal_mol"]
        if len(reference) != 1 or pd.isna(reference.iloc[0]):
            raise ValueError(f"Reference family {reference_family!r} is unavailable for {reaction_id}")
        panel = group.copy()
        panel["panel"] = panel_name
        panel["reference_family"] = reference_family
        panel["relative_g_kcal_mol"] = panel["delta_g_kcal_mol"] - float(reference.iloc[0])
        panel["rank_within_panel"] = panel["delta_g_kcal_mol"].rank(method="min").astype("Int64")
        panel["panel_family_count"] = len(panel)
        parts.append(panel)
    return pd.concat(parts, ignore_index=True)
