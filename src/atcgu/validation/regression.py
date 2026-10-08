"""Compare computed tables with the published results by scientific keys."""
from __future__ import annotations
from pathlib import Path
from typing import Callable
import numpy as np
import pandas as pd

HARTREE_TOLERANCE = 1e-9
KCAL_TOLERANCE = 1e-6


def compare_table(actual: pd.DataFrame, published: Path, keys: list[str], columns: list[str], tolerance: float = HARTREE_TOLERANCE,
                  rows: Callable[[pd.DataFrame], pd.Series] | None = None, label: str = "") -> dict:
    """Outer-join on ``keys`` and compare ``columns``; ``rows`` optionally restricts both sides."""
    expected = pd.read_csv(published, float_precision="round_trip", dtype={key: str for key in keys})
    check_id = f"{published.parent.name}/{published.name}" + (f" [{label}]" if label else "")
    if rows is not None:
        actual, expected = actual.loc[rows(actual)], expected.loc[rows(expected)]
    if actual.duplicated(keys).any() or expected.duplicated(keys).any():
        raise ValueError(f"Nonunique comparison keys: {check_id}")
    left, right = actual[keys + columns].copy(), expected[keys + columns].copy()
    for key in keys:
        left[key], right[key] = left[key].astype(str), right[key].astype(str)
    merged = left.merge(right, on=keys, how="outer", suffixes=("_actual", "_expected"), indicator=True, validate="one_to_one")
    if not merged._merge.eq("both").all():
        raise ValueError(f"Missing/extra benchmark rows in {check_id}")
    max_difference = 0.0
    for column in columns:
        a, b = merged[column + "_actual"], merged[column + "_expected"]
        if pd.api.types.is_bool_dtype(a) or pd.api.types.is_bool_dtype(b):
            if not a.astype(str).str.lower().eq(b.astype(str).str.lower()).all():
                raise ValueError(f"Benchmark mismatch {check_id}:{column}")
        elif pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            if not np.isfinite(a).all() or not np.isfinite(b).all():
                raise ValueError(f"Nonfinite benchmark field {check_id}:{column}")
            error = float((a - b).abs().max()) if len(a) else 0.0
            max_difference = max(max_difference, error)
            if error > tolerance:
                raise ValueError(f"Benchmark mismatch {check_id}:{column}: {error}")
        elif not a.fillna("").astype(str).eq(b.fillna("").astype(str)).all():
            raise ValueError(f"Benchmark mismatch {check_id}:{column}")
    return {"check_id": check_id, "passed": True, "rows": len(merged), "maximum_numeric_difference": max_difference, "tolerance": tolerance}


def _not_failed(frame: pd.DataFrame) -> pd.Series:
    return frame.status.ne("failed")


def compare_run(root: Path, outputs: dict[str, pd.DataFrame], plotted: Path | None = None,
                si_plotted: Path | None = None) -> pd.DataFrame:
    """Every table comparison runs with or without figures; figures add plotted-point checks."""
    continuum = root / "data/published/continuum"
    micro = root / "data/published/figure1"
    si = root / "data/published/si_microsolvation"
    identity = outputs["01_input_inventory"]
    usable = set(identity.loc[identity.is_usable, "calculation_id"])
    counts = ["base_candidate_count", "nucleoside_candidate_count", "nucleotide_candidate_count", "ribose_candidate_count",
              "phosphate_candidate_count", "water_candidate_count", "ribo_phosphate_candidate_count", "required_species_count",
              "complete_required_coverage"]
    structure = ["family", "panel_role", "placement", "structural_status", "structural_failure", "fragment_sizes",
                 "intact_water_count", "direct_fold", "near_fold", "fold_state", "phosphate_to_water_hb_count",
                 "water_to_base_hb_count", "bridging_water_count", "structure_matches_manifest",
                 "association_shift_kcal_mol", "common_reference_value_kcal_mol",
                 "fixed_reference_score_kcal_mol", "selected_for_family"]
    checks = [
        compare_table(outputs["07_reaction_energies"], continuum / "07_reaction_energies.csv", ["analysis_id", "reaction_family_id", "reaction_id"], ["delta_g_kcal_mol"]),
        compare_table(outputs["09_relative_rankings"], continuum / "09_relative_rankings.csv", ["analysis_id", "reaction_family_id", "reaction_id"], ["relative_g_kcal_mol", "rank_within_panel"]),
        compare_table(outputs["05_selected_species"], continuum / "05_selected_species.csv", ["analysis_id", "species_id"], ["calculation_id", "computed_gibbs_hartree", "candidate_count"]),
        compare_table(outputs["04_analysis_candidates"], continuum / "04_analysis_candidates.csv", ["analysis_id", "calculation_id"], ["candidate_rank_within_species", "selected_for_analysis", "computed_gibbs_hartree"]),
        compare_table(outputs["10_legacy_vs_transferred_audit"], continuum / "10_legacy_vs_transferred_audit.csv", ["method_id", "species_id"], ["transferred_minus_legacy_kcal_mol", "primary_policy_selected_file"]),
        compare_table(identity, continuum / "01_input_inventory.csv", ["calculation_id"], ["method_id", "analysis_role", "primary_include", "manifest_species_id", "local_path", "expected_sha256", "parsed_sha256", "hash_matches_manifest"]),
        compare_table(outputs["02_log_validation"], continuum / "02_log_validation.csv", ["calculation_id"], ["status", "is_usable"]),
        # `failures` is compared on all rows (no failed-source records remain in the release).
        compare_table(outputs["02_log_validation"], continuum / "02_log_validation.csv", ["calculation_id"], ["failures"], rows=_not_failed, label="failures, non-failed rows"),
        # Energies are compared on usable rows (all records in the release are usable).
        compare_table(outputs["03_thermochemistry"], continuum / "03_thermochemistry.csv", ["calculation_id"],
                      ["electronic_energy_hartree", "thermal_gibbs_correction_hartree", "computed_gibbs_hartree", "printed_gibbs_hartree", "gibbs_residual_hartree"],
                      rows=lambda frame: frame.calculation_id.isin(usable), label="usable rows"),
        compare_table(outputs["06_sampling_coverage"], continuum / "06_sampling_coverage.csv", ["analysis_id", "reaction_family_id"], counts),
        compare_table(outputs["08_path_closure"], continuum / "08_path_closure.csv", ["analysis_id", "reaction_family_id"], ["path1_step_sum_error_kcal_mol", "path2_step_sum_error_kcal_mol", "net_path_difference_kcal_mol", "passes"]),
        compare_table(outputs["micro_02_reference_manifest_and_validation"], micro / "02_reference_manifest_and_validation.csv", ["family", "role"], ["calculation_id", "species_id", "expected_sha256", "computed_gibbs_hartree"]),
        compare_table(outputs["micro_03_fixed_reference_convention"], micro / "03_fixed_reference_convention.csv", ["reference_family"],
                      ["harmonic_continuum_nucleotide_g_hartree", "harmonic_base_g_hartree",
                       "harmonic_continuum_nucleotide_minus_base_hartree", "harmonic_water_g_hartree",
                       "explicit_water_count", "temperature_k", "display_hartree_to_kcal_mol",
                       "standard_state_water_reservoir_correction_kcal_mol", "convention"], KCAL_TOLERANCE),
        compare_table(outputs["micro_04_structure_validation"], micro / "04_structure_validation.csv", ["candidate_id"], structure, KCAL_TOLERANCE),
        compare_table(outputs["micro_05_primary_candidate_energies"], micro / "05_primary_candidate_energies.csv", ["candidate_id"], ["g_hartree", "fold_state", "intact_water_count", "phosphate_to_water_hb_count", "water_to_base_hb_count", "fixed_reference_score_kcal_mol"], KCAL_TOLERANCE),
        compare_table(outputs["micro_06_primary_selected_candidates"], micro / "06_primary_selected_candidates.csv", ["family"], ["candidate_id", "relative_g_kcal_mol", "fixed_reference_score_kcal_mol", "rank"], KCAL_TOLERANCE),
        compare_table(outputs["micro_16_panel_c_plotted_data"], micro / "16_panel_c_plotted_data.csv", ["family"], ["continuum_relative_g_kcal_mol", "microsolvated_relative_g_kcal_mol", "selected_candidate"], KCAL_TOLERANCE),
        compare_table(outputs["si_03_treatment_energies"], si / "03_treatment_energies.csv", ["treatment", "calculation_id"], ["treated_g_hartree", "entropy_shift_hartree"]),
        compare_table(outputs["si_05_fixed_reference_placement_scores"], si / "05_fixed_reference_placement_scores.csv", ["treatment", "family", "placement"],
                      ["candidate_id", "calculation_id", "score_hartree", "fixed_reference_score_kcal_mol",
                       "selected_for_family", "selection_rank_within_pair"], KCAL_TOLERANCE),
        compare_table(outputs["si_06_selected_candidates"], si / "06_selected_candidates.csv", ["treatment", "family"], ["candidate_id", "relative_g_kcal_mol", "fixed_reference_score_kcal_mol", "rank"], KCAL_TOLERANCE),
        compare_table(outputs["si_07_continuum_rankings"], si / "07_continuum_rankings.csv", ["treatment", "family"], ["relative_g_kcal_mol", "rank"], KCAL_TOLERANCE),
        compare_table(outputs["si_08_selected_fixed_reference_shifts"], si / "08_selected_fixed_reference_shifts.csv", ["treatment", "family"],
                      ["selected_candidate_id", "selected_placement", "harmonic_selected_candidate_id",
                       "harmonic_selected_placement", "selected_fixed_reference_score_kcal_mol",
                       "harmonic_selected_fixed_reference_score_kcal_mol",
                       "selected_fixed_reference_shift_from_harmonic_kcal_mol", "placement_switched_from_harmonic"], KCAL_TOLERANCE),
        compare_table(outputs["si_09_placement_combinations"], si / "09_placement_combinations.csv", ["treatment", "choice_bits"], ["rank_order", "c_minus_a_kcal_mol", "pairwise_retained", "ac_top_two"], KCAL_TOLERANCE),
        compare_table(outputs["si_10_rank_occupancy"], si / "10_rank_occupancy.csv", ["treatment", "family", "rank"],
                      ["combination_count", "selected_lower_g_rank", "is_selected_lower_g_rank"]),
        compare_table(outputs["si_thermochemical_summary"], si / "thermochemical_summary.csv", ["treatment"], ["selected_order", "c_minus_a", "pairwise_retained", "unique_placement_orders", "ac_top_two_choices"], HARTREE_TOLERANCE),
    ]
    if plotted:
        checks.append(compare_table(pd.read_csv(plotted), micro / "figure1_plotted_data.csv", ["panel", "family", "series"], ["value"], KCAL_TOLERANCE))
    if si_plotted:
        checks.append(compare_table(pd.read_csv(si_plotted), root / "results/supplement/si_microsolvation_plotted_data.csv",
                                    ["panel", "family", "treatment", "placement", "series"],
                                    ["selected_for_family", "value"], KCAL_TOLERANCE))
    return pd.DataFrame(checks)
