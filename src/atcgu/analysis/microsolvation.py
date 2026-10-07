"""Primary Figure 1 two-water comparison from portable scientific records."""
from __future__ import annotations
from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
import yaml

from ..constants import FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL, HARTREE_TO_KCAL_MOL
from ..geometry import bond_graph, components, fingerprint
from ..records import read_bundle
from ..validation.records import check_record
from .structure import analyze_microsolvation_structure

FAMILIES = {"adenine", "cytosine", "guanine", "tap", "melamine", "imidazole"}
STRUCTURE_COLUMNS = ["family", "candidate_id", "panel_role", "placement", "structural_status", "structural_failure",
                     "fragment_sizes", "intact_water_count", "direct_fold", "near_fold", "fold_state",
                     "phosphate_to_water_hb_count", "water_to_base_hb_count", "bridging_water_count",
                     "structure_matches_manifest", "association_shift_kcal_mol", "common_reference_value_kcal_mol",
                     "fixed_reference_score_kcal_mol",
                     "selected_for_family"]


def _load(root: Path, row: dict, inventory: pd.DataFrame, method: dict) -> dict:
    matches = inventory.loc[inventory.calculation_id.eq(row["calculation_id"])]
    if len(matches) != 1 or matches.iloc[0].record_path != row["record_path"]:
        raise ValueError("Primary hydration record identity/path disagrees with the released bundle inventory")
    record = read_bundle(root, row["record_path"], matches.iloc[0].record_sha256)
    if record["calculation_id"] != row["calculation_id"]:
        raise ValueError("Primary hydration calculation ID mismatch")
    failures = check_record(record, method)
    if record["source_log_sha256"] != row["expected_sha256"]:
        failures.append("source_hash_mismatch")
    if failures:
        raise ValueError(f"Invalid primary hydration record {row['calculation_id']}: {failures}")
    return record


def run_primary(root: Path, tables: Path, selected: pd.DataFrame) -> dict:
    """Write the primary tables; also return each validated record under ``records``.

    ``records`` maps calculation_id to its checked bundle so later steps
    (the SI treatments) reuse these identities instead of re-validating.
    """
    spec = pd.read_csv(root / "config/microsolvation_candidates.csv", dtype=str, keep_default_na=False)
    if set(spec.family) != FAMILIES or len(spec) != 12:
        raise ValueError("Primary hydration panel requires exactly six families and twelve candidates")
    if spec.candidate_id.duplicated().any() or not spec.groupby("family").placement.apply(lambda s: set(s) == {"splitA", "splitB"}).all():
        raise ValueError("Each hydration family requires unique splitA and splitB placements")
    references = pd.read_csv(root / "config/microsolvation_references.csv", dtype=str, keep_default_na=False)
    if references.duplicated(["family", "role"]).any() or set(references.family) != FAMILIES:
        raise ValueError("Hydration references are incomplete or nonunique")
    inventory = pd.read_csv(root / "data/calculations/manifest.csv", dtype=str)
    method = yaml.safe_load((root / "config/methods.yaml").read_text())["methods"]["1_b3lyp_pcm"]
    reference_family = json.loads((root / "config/study.json").read_text())["reference_family"]
    figure_config = yaml.safe_load((root / "config/figures.yaml").read_text())
    loaded = {}
    refrows = []
    for row in references.to_dict("records"):
        record = loaded[row["calculation_id"]] = _load(root, row, inventory, method)
        identity = {k: v for k, v in row.items() if k != "folded_start"}  # folded_start is a preparation setting
        refrows.append({**identity, "computed_gibbs_hartree": record["computed_gibbs_hartree"], "is_usable": True})
    refs = pd.DataFrame(refrows)
    for role in ["base", "continuum_nucleotide"]:
        if set(refs.loc[refs.role.eq(role), "family"]) != FAMILIES:
            raise ValueError(f"Incomplete hydration reference role: {role}")
    base = refs.loc[refs.role.eq("base")].set_index("family").computed_gibbs_hartree
    dry = refs.loc[refs.role.eq("continuum_nucleotide")].set_index("family").computed_gibbs_hartree
    full = selected.loc[selected.analysis_id.eq(figure_config["figure1"]["panel_a_analysis_id"])].set_index("species_id")
    for row in references.to_dict("records"):
        if full.loc[row["species_id"], "calculation_id"] != row["calculation_id"]:
            raise ValueError("Primary hydration references must match the selected continuum structures")
    water_g = float(full.loc["water", "computed_gibbs_hartree"])
    rows = []
    for row in spec.to_dict("records"):
        record = loaded[row["calculation_id"]] = _load(root, row, inventory, method)
        geometry = record["geometry"]
        z = np.array([atom[0] for atom in geometry])
        xyz = np.array([atom[1:] for atom in geometry])
        _, adjacency = bond_graph(z, xyz)
        solutes = [c for c in components(adjacency) if any(z[i] == 15 for i in c)]
        dry_ref = references.loc[references.family.eq(row["family"]) & references.role.eq("continuum_nucleotide")].iloc[0]
        dry_record = loaded[dry_ref.calculation_id]
        if len(solutes) != 1 or fingerprint([geometry[i] for i in sorted(solutes[0])]) != dry_record["connectivity_fingerprint"]:
            raise ValueError("Primary cluster solute does not match its intended continuum nucleotide identity")
        metrics = analyze_microsolvation_structure(root / row["record_path"].replace("record.json", record["geometry_file"]))
        if metrics.get("structural_status") != "valid":
            raise ValueError(f"Cannot classify primary candidate {row['candidate_id']}: {metrics}")
        checks = [metrics["intact_water_count"] == int(row["expected_water_count"]),
                  metrics["fold_state"] == row["expected_fold_state"],
                  metrics["direct_fold"] == (row["expected_direct_fold"].lower() == "true"),
                  metrics["near_fold"] == (row["expected_near_fold"].lower() == "true"),
                  metrics["phosphate_to_water_hb_count"] >= int(row["expected_phosphate_water_hb_min"]),
                  metrics["water_to_base_hb_count"] >= int(row["expected_base_water_hb_min"])]
        if not all(checks):
            raise ValueError(f"Primary structure violates declared criteria: {row['candidate_id']}")
        rows.append({**row, **metrics, "structure_matches_manifest": True, "is_usable": True, "treatment": "harmonic", "g_hartree": record["computed_gibbs_hartree"]})
    energies = pd.DataFrame(rows)
    energies["base_g_hartree"] = energies.family.map(base)
    energies["score_hartree"] = energies.g_hartree - energies.base_g_hartree
    chosen = energies.loc[energies.groupby("family").g_hartree.idxmin()].copy()
    if set(chosen.family) != FAMILIES or len(chosen) != 6:
        raise ValueError("Incomplete harmonic treatment selection")
    reference = float(chosen.loc[chosen.family.eq(reference_family), "score_hartree"].iloc[0])
    chosen["relative_g_kcal_mol"] = (chosen.score_hartree - reference) * HARTREE_TO_KCAL_MOL
    chosen["rank"] = chosen.relative_g_kcal_mol.rank(method="min").astype(int)
    dry_score = dry - base
    dry_relative = (dry_score - dry_score[reference_family]) * HARTREE_TO_KCAL_MOL
    config = figure_config["panel_c"]
    rt = float(config["gas_constant_kcal_mol_k"]) * float(config["temperature_k"])
    offset = -2 * rt * math.log(float(config["gas_to_solution_factor"])) - 2 * rt * math.log(float(config["water_reservoir_molarity"]))
    fixed_reference_hartree = float(dry[reference_family] - base[reference_family])
    fixed_reference = pd.DataFrame([{
        "reference_family": reference_family,
        "harmonic_continuum_nucleotide_g_hartree": float(dry[reference_family]),
        "harmonic_base_g_hartree": float(base[reference_family]),
        "harmonic_continuum_nucleotide_minus_base_hartree": fixed_reference_hartree,
        "harmonic_water_g_hartree": water_g,
        "explicit_water_count": 2,
        "temperature_k": float(config["temperature_k"]),
        "display_hartree_to_kcal_mol": FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL,
        "standard_state_water_reservoir_correction_kcal_mol": offset,
        "convention": "fixed harmonic continuum-adenine nucleotide-minus-base score; fixed harmonic isolated-water reservoir",
    }])
    energies["association_shift_kcal_mol"] = (energies.g_hartree - energies.family.map(dry) - 2 * water_g) * FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL + offset
    energies["common_reference_value_kcal_mol"] = energies.family.map(dry_relative) + energies.association_shift_kcal_mol
    energies["fixed_reference_score_kcal_mol"] = (
        energies.score_hartree - fixed_reference_hartree - 2 * water_g
    ) * FIGURE1C_DISPLAY_HARTREE_TO_KCAL_MOL + offset
    if not np.allclose(energies.fixed_reference_score_kcal_mol, energies.common_reference_value_kcal_mol, atol=1e-6, rtol=0):
        raise ValueError("Fixed-reference score disagrees with the Figure 1C common-zero construction")
    energies["selected_for_family"] = energies.candidate_id.isin(chosen.candidate_id)
    chosen = chosen.drop(columns=["fixed_reference_score_kcal_mol"], errors="ignore").merge(
        energies[["candidate_id", "fixed_reference_score_kcal_mol"]], on="candidate_id", validate="one_to_one"
    )
    panel = pd.DataFrame({"family": dry_relative.index, "continuum_relative_g_kcal_mol": dry_relative.values})
    panel["continuum_rank"] = panel.continuum_relative_g_kcal_mol.rank(method="min").astype(int)
    panel = panel.merge(chosen[["family", "candidate_id", "placement", "fold_state", "relative_g_kcal_mol", "rank"]], on="family", validate="one_to_one")
    panel = panel.rename(columns={"candidate_id": "selected_candidate", "placement": "selected_placement", "relative_g_kcal_mol": "microsolvated_relative_g_kcal_mol", "rank": "microsolvated_rank"})
    panel = panel.sort_values("continuum_rank").reset_index(drop=True)
    outputs = {"micro_02_reference_manifest_and_validation": refs, "micro_03_fixed_reference_convention": fixed_reference,
               "micro_04_structure_validation": energies[STRUCTURE_COLUMNS],
               "micro_05_primary_candidate_energies": energies, "micro_06_primary_selected_candidates": chosen,
               "micro_16_panel_c_plotted_data": panel}
    for name, frame in outputs.items():
        frame.to_csv(tables / f"{name}.csv", index=False)
    outputs["records"] = loaded
    return outputs
