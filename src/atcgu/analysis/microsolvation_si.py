"""SI fixed-reference placement scores and the 64 placement combinations.

Every displayed placement score uses Figure 1C's fixed harmonic continuum-
adenine, isolated-water and standard-state/water-reservoir terms. Only the
cluster and its corresponding base receive the requested entropy treatment.
Treatment-relative adenine scores are retained as diagnostics because the
64-vector rankings are invariant to a common zero.
"""
from __future__ import annotations
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..constants import FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL, HARTREE_TO_KCAL_MOL
from .microsolvation import FAMILIES
from .thermochemistry_sensitivity import TREATMENTS, treated_energies

PLACEMENTS = ("splitA", "splitB")


def _settings(root: Path) -> tuple[list[str], list[str], list[str], str]:
    block = json.loads((root / "config/study.json").read_text())["si_microsolvation"]
    treatments, order = list(block["treatments"]), list(block["family_order"])
    display_order = list(block["display_family_order"])
    if tuple(treatments) != TREATMENTS:
        raise ValueError(f"SI treatments must be {TREATMENTS}")
    if set(order) != FAMILIES or len(order) != len(FAMILIES):
        raise ValueError("SI family order must list each hydration family once")
    if set(display_order) != FAMILIES or len(display_order) != len(FAMILIES):
        raise ValueError("SI display family order must list each hydration family once")
    if block["placement_bits"] != list(PLACEMENTS):
        raise ValueError("SI placement bits must be [splitA, splitB]")
    return treatments, order, display_order, block["reference_family"]


def _energy_rows(primary: dict) -> pd.DataFrame:
    records = primary["records"]
    candidates = primary["micro_05_primary_candidate_energies"]
    references = primary["micro_02_reference_manifest_and_validation"]
    items = [("candidate", r["family"], "two_water_cluster", r["candidate_id"], r) for r in candidates.to_dict("records")]
    items += [("reference", r["family"], r["role"], "", r) for r in references.to_dict("records")]
    treated = {r["calculation_id"]: treated_energies(records[r["calculation_id"]]) for *_, r in items}
    rows = []
    for treatment in TREATMENTS:
        for record_type, family, role, candidate_id, row in items:
            record = records[row["calculation_id"]]
            energies = treated[row["calculation_id"]]
            rows.append({"treatment": treatment, "record_type": record_type, "family": family, "role": role,
                         "candidate_id": candidate_id, "calculation_id": row["calculation_id"],
                         "record_path": row["record_path"], "harmonic_g_hartree": energies["harmonic"],
                         "entropy_shift_hartree": energies[treatment] - energies["harmonic"],
                         "treated_g_hartree": energies[treatment], "temperature_k": record["temperature_k"],
                         "frequency_count": record["frequency_count"]})
    frame = pd.DataFrame(rows)
    if frame.record_path.str.startswith("/").any() or frame.duplicated(["treatment", "calculation_id"]).any():
        raise ValueError("SI energies require unique, repository-relative records")
    return frame


def _rank(frame: pd.DataFrame, reference_family: str) -> pd.DataFrame:
    reference = float(frame.loc[frame.family.eq(reference_family), "score_hartree"].iloc[0])
    frame["relative_g_kcal_mol"] = (frame.score_hartree - reference) * HARTREE_TO_KCAL_MOL
    frame["rank"] = frame.score_hartree.rank(method="first").astype(int)
    return frame.sort_values("rank")


def run_si(root: Path, tables: Path, primary: dict) -> dict[str, pd.DataFrame]:
    treatments, order, display_order, reference_family = _settings(root)
    energies = _energy_rows(primary)
    placements = primary["micro_05_primary_candidate_energies"].set_index("candidate_id").placement
    convention = primary["micro_03_fixed_reference_convention"]
    if len(convention) != 1 or convention.reference_family.iloc[0] != reference_family:
        raise ValueError("SI requires one fixed Figure 1C reference matching the configured family")
    fixed = convention.iloc[0]
    if float(fixed.display_hartree_to_kcal_mol) != FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL:
        raise ValueError("SI fixed-reference conversion disagrees with Figure 1C")
    fixed_zero = float(fixed.harmonic_continuum_nucleotide_minus_base_hartree)
    fixed_water = float(fixed.harmonic_water_g_hartree)
    fixed_correction = float(fixed.standard_state_water_reservoir_correction_kcal_mol)
    selected_parts, dry_parts, placement_parts, vector_rows, occupancy_rows, summary_rows = [], [], [], [], [], []
    for treatment in treatments:
        group = energies.loc[energies.treatment.eq(treatment)]
        refs = group.loc[group.record_type.eq("reference")].set_index(["family", "role"])
        base = refs.xs("base", level="role")
        dry = refs.xs("continuum_nucleotide", level="role")
        clusters = group.loc[group.record_type.eq("candidate")].copy()
        clusters["placement"] = clusters.candidate_id.map(placements)
        clusters["score_hartree"] = clusters.treated_g_hartree - clusters.family.map(base.treated_g_hartree)
        clusters["fixed_reference_score_kcal_mol"] = (
            clusters.score_hartree - fixed_zero - 2 * fixed_water
        ) * FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL + fixed_correction
        if len(clusters) != 12 or len(base) != 6 or len(dry) != 6:
            raise ValueError(f"Incomplete SI records for {treatment}")
        chosen = clusters.sort_values(["family", "score_hartree", "candidate_id"]).groupby("family").head(1)
        chosen = _rank(chosen.copy(), reference_family)
        chosen["treatment"] = treatment
        selected_parts.append(chosen[["treatment", "family", "candidate_id", "calculation_id", "placement",
                                      "score_hartree", "relative_g_kcal_mol", "fixed_reference_score_kcal_mol", "rank"]])
        chosen_ids = set(chosen.candidate_id)
        clusters["selected_for_family"] = clusters.candidate_id.isin(chosen_ids)
        clusters = clusters.sort_values(["family", "score_hartree", "candidate_id"])
        clusters["selection_rank_within_pair"] = clusters.groupby("family").cumcount() + 1
        clusters["harmonic_continuum_adenine_minus_base_hartree"] = fixed_zero
        clusters["harmonic_water_g_hartree"] = fixed_water
        clusters["standard_state_water_reservoir_correction_kcal_mol"] = fixed_correction
        placement_parts.append(clusters[["treatment", "family", "candidate_id", "calculation_id", "placement",
                                         "treated_g_hartree", "score_hartree", "fixed_reference_score_kcal_mol",
                                         "selected_for_family", "selection_rank_within_pair",
                                         "harmonic_continuum_adenine_minus_base_hartree", "harmonic_water_g_hartree",
                                         "standard_state_water_reservoir_correction_kcal_mol"]])
        continuum = pd.DataFrame({"family": order, "base_calculation_id": [base.calculation_id[f] for f in order],
                                  "dry_calculation_id": [dry.calculation_id[f] for f in order],
                                  "score_hartree": [dry.treated_g_hartree[f] - base.treated_g_hartree[f] for f in order]})
        continuum = _rank(continuum, reference_family)
        continuum.insert(0, "treatment", treatment)
        dry_parts.append(continuum)
        dry_score = continuum.set_index("family").score_hartree
        options = {f: clusters.loc[clusters.family.eq(f)].set_index("placement") for f in order}

        def retained(scores: dict) -> int:
            return int(sum(np.sign(dry_score[a] - dry_score[b]) == np.sign(scores[a] - scores[b])
                           for a, b in combinations(order, 2)))

        rank_counts = {(family, rank): 0 for family in order for rank in range(1, len(order) + 1)}
        for bits in product([0, 1], repeat=len(order)):
            picks = {f: options[f].loc[PLACEMENTS[bit]] for f, bit in zip(order, bits)}
            scores = {f: float(pick.score_hartree) for f, pick in picks.items()}
            ranked = sorted(order, key=scores.get)
            for rank, family in enumerate(ranked, start=1):
                rank_counts[family, rank] += 1
            vector_rows.append({"treatment": treatment, "choice_bits": "".join(map(str, bits)),
                                **{f"{f}_candidate_id": picks[f].candidate_id for f in order},
                                "rank_order": " < ".join(ranked),
                                "c_minus_a_kcal_mol": (scores["cytosine"] - scores["adenine"]) * HARTREE_TO_KCAL_MOL,
                                "pairwise_retained": retained(scores),
                                "ac_top_two": set(ranked[:2]) == {"adenine", "cytosine"}})
        vectors = pd.DataFrame(vector_rows[-2 ** len(order):])
        if vectors.choice_bits.nunique() != 64:
            raise ValueError("Expected 64 unique placement vectors")
        picked = chosen.set_index("family")
        selected_ranks = {family: rank for rank, family in enumerate(chosen.family, start=1)}
        for family in display_order:
            for rank in range(1, len(order) + 1):
                occupancy_rows.append({"treatment": treatment, "family": family, "rank": rank,
                                       "combination_count": rank_counts[family, rank],
                                       "selected_lower_g_rank": selected_ranks[family],
                                       "is_selected_lower_g_rank": rank == selected_ranks[family]})
        if any(sum(rank_counts[family, rank] for rank in range(1, len(order) + 1)) != 64 for family in order):
            raise ValueError(f"Rank occupancy does not sum to 64 for {treatment}")
        summary_rows.append({"treatment": treatment, "selected_order": " < ".join(chosen.family),
                             "c_minus_a": float(picked.relative_g_kcal_mol["cytosine"]),
                             "pairwise_retained": retained(picked.score_hartree.to_dict()),
                             "unique_placement_orders": int(vectors.rank_order.nunique()),
                             "ac_top_two_choices": int(vectors.ac_top_two.sum())})
    selected_all = pd.concat(selected_parts, ignore_index=True)
    placement_all = pd.concat(placement_parts, ignore_index=True)
    harmonic_selected = selected_all.loc[selected_all.treatment.eq("harmonic")].set_index("family")
    shift_rows = []
    for row in selected_all.itertuples(index=False):
        harmonic = harmonic_selected.loc[row.family]
        shift_rows.append({"treatment": row.treatment, "family": row.family,
                           "selected_candidate_id": row.candidate_id, "selected_placement": row.placement,
                           "harmonic_selected_candidate_id": harmonic.candidate_id,
                           "harmonic_selected_placement": harmonic.placement,
                           "selected_fixed_reference_score_kcal_mol": row.fixed_reference_score_kcal_mol,
                           "harmonic_selected_fixed_reference_score_kcal_mol": harmonic.fixed_reference_score_kcal_mol,
                           "selected_fixed_reference_shift_from_harmonic_kcal_mol":
                               row.fixed_reference_score_kcal_mol - harmonic.fixed_reference_score_kcal_mol,
                           "placement_switched_from_harmonic": row.candidate_id != harmonic.candidate_id})
    shifts = pd.DataFrame(shift_rows)
    harmonic_primary = primary["micro_05_primary_candidate_energies"].set_index("candidate_id").fixed_reference_score_kcal_mol
    harmonic_si = placement_all.loc[placement_all.treatment.eq("harmonic")].set_index("candidate_id").fixed_reference_score_kcal_mol
    if not np.allclose(harmonic_si.loc[harmonic_primary.index], harmonic_primary, atol=1e-8, rtol=0):
        raise ValueError("Harmonic SI fixed-reference endpoints disagree with Figure 1C")
    outputs = {"si_03_treatment_energies": energies,
               "si_05_fixed_reference_placement_scores": placement_all,
               "si_06_selected_candidates": selected_all,
               "si_07_continuum_rankings": pd.concat(dry_parts, ignore_index=True),
               "si_08_selected_fixed_reference_shifts": shifts,
               "si_09_placement_combinations": pd.DataFrame(vector_rows),
               "si_10_rank_occupancy": pd.DataFrame(occupancy_rows),
               "si_thermochemical_summary": pd.DataFrame(summary_rows).sort_values("treatment").reset_index(drop=True)}
    expected = {"si_03_treatment_energies": 96, "si_05_fixed_reference_placement_scores": 48,
                "si_06_selected_candidates": 24, "si_07_continuum_rankings": 24,
                "si_08_selected_fixed_reference_shifts": 24, "si_09_placement_combinations": 256,
                "si_10_rank_occupancy": 144, "si_thermochemical_summary": 4}
    if any(len(outputs[name]) != rows for name, rows in expected.items()):
        raise ValueError("SI tables have unexpected row counts")
    folder = tables / "si"
    folder.mkdir(exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(folder / f"{name.removeprefix('si_')}.csv", index=False)
    return outputs
