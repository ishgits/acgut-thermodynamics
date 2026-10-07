"""Continuum analysis from checked portable calculation records, without raw logs."""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd
import yaml

from ..constants import HARTREE_TO_KCAL_MOL, check_study_conventions
from ..records import read_bundle
from ..validation.records import check_record
from . import selection as policies


def load_inventory(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Return (identity, analysis sets, registry, reaction families, family panels).

    Study panels are validated before any record is read, so a shrunken registry fails fast.
    """
    sets = pd.read_csv(root / "config/analysis_sets.csv", dtype=str, keep_default_na=False)
    families = pd.read_csv(root / "config/reaction_families.csv", dtype=str, keep_default_na=False)
    study = json.loads((root / "config/study.json").read_text())
    check_study_conventions(study, yaml.safe_load((root / "config/figures.yaml").read_text()))
    missing = {panel: sorted(set(members) - set(families.reaction_family_id)) for panel, members in study["family_panels"].items()}
    missing = {panel: members for panel, members in missing.items() if members}
    if missing:
        raise ValueError(f"config/study.json family panels name families absent from config/reaction_families.csv: {missing}")
    panels = {"all_21": families.reaction_family_id.tolist(), **study["family_panels"]}
    for definition in sets.to_dict("records"):
        size = len(panels[definition["family_panel"]])
        if size != int(definition["family_count"]):
            raise ValueError(f"{definition['analysis_id']} expects {definition['family_count']} families in panel "
                             f"{definition['family_panel']}; found {size}")
    manifest = pd.read_csv(root / "config/calculation_manifest.csv", dtype=str, keep_default_na=False)
    if manifest.calculation_id.duplicated().any() or manifest.local_path.duplicated().any():
        raise ValueError("Continuum manifest contains duplicate calculation records")
    methods = yaml.safe_load((root / "config/methods.yaml").read_text())["methods"]
    registry = pd.read_csv(root / "config/molecules.csv", dtype=str, keep_default_na=False)
    if registry.species_id.duplicated().any():
        raise ValueError("Species registry contains duplicate IDs")
    states = registry.set_index("species_id")[["charge", "multiplicity"]].to_dict("index")
    parsed = []
    for row in manifest.to_dict("records"):
        record = read_bundle(root, row["local_path"], row["record_sha256"])
        if record["calculation_id"] != row["calculation_id"]:
            raise ValueError("Continuum calculation ID disagrees with the configured record")
        expected = dict(methods[row["method_id"]])
        if row["species_id"] not in states:
            raise KeyError(f"Species {row['species_id']!r} in config/calculation_manifest.csv is absent from config/molecules.csv")
        expected.update(states[row["species_id"]])
        errors = check_record(record, expected)
        flat = {k: v for k, v in record.items() if not isinstance(v, (list, dict))}
        flat["parsed_sha256"] = flat.pop("sha256")
        flat["parsed_file_name"] = flat.pop("file_name")
        if errors:
            flat["failures"] = "; ".join(filter(None, [flat.get("failures", ""), *errors]))
            flat["status"] = "failed"
        flat["local_path"] = row["local_path"]
        parsed.append(flat)
    manifest["primary_include"] = policies._as_bool(manifest["primary_include"])
    manifest = manifest.rename(columns={"species_id": "manifest_species_id", "file_name": "manifest_file_name",
                                       "sha256": "expected_sha256"})
    frame = manifest.merge(pd.DataFrame(parsed).drop(columns="calculation_id"), on="local_path", validate="one_to_one")
    frame = policies._resolve_and_validate(frame, registry)
    return frame, sets, registry, families, panels


def legacy_audit(identity: pd.DataFrame) -> pd.DataFrame:
    """Compare each transferred structure with its diagnostic-only alternative; the alternative is never selected."""
    alternate = identity.loc[identity.method_id.ne("1_b3lyp_pcm") & identity.is_usable]
    new = alternate.loc[alternate.analysis_role.isin(["required_transferred_nucleotide", "transferred_nucleoside"])]
    old = alternate.loc[alternate.analysis_role.eq("diagnostic_only_legacy_product")]
    audit = new[["method_id", "species_id", "molecule_class", "manifest_file_name", "computed_gibbs_hartree"]].rename(
        columns={"manifest_file_name": "transferred_file", "computed_gibbs_hartree": "transferred_g_hartree"})
    audit = audit.merge(old[["method_id", "species_id", "manifest_file_name", "computed_gibbs_hartree"]].rename(
        columns={"manifest_file_name": "legacy_file", "computed_gibbs_hartree": "legacy_g_hartree"}),
        on=["method_id", "species_id"], validate="one_to_one")
    audit["transferred_minus_legacy_kcal_mol"] = (audit.transferred_g_hartree - audit.legacy_g_hartree) * HARTREE_TO_KCAL_MOL
    audit["numerically_lower_file"] = audit.apply(lambda r: r.transferred_file if r.transferred_minus_legacy_kcal_mol < 0 else r.legacy_file, axis=1)
    audit["primary_policy_selected_file"] = audit.transferred_file
    audit["primary_policy"] = "required_transferred_structure"
    return audit.sort_values(["molecule_class", "method_id", "species_id"]).reset_index(drop=True)


def run_continuum(root: Path, tables: Path) -> dict[str, pd.DataFrame]:
    identity, sets, registry, families, panels = load_inventory(root)
    candidates, selected, reaction_map = policies._analysis_inputs(identity, sets, registry, families, panels)
    for name, frame in [("candidates", candidates), ("selected", selected)]:
        for definition in sets.to_dict("records"):
            mask = frame.analysis_id.eq(definition["analysis_id"])
            frame.loc[mask, "method_label"] = definition["method_label"]
    coverage = policies._build_coverage(candidates, selected, sets, families, panels, reaction_map)
    reference_family = json.loads((root / "config/study.json").read_text())["reference_family"]
    full_panel_id = yaml.safe_load((root / "config/figures.yaml").read_text())["figure1"]["panel_a_analysis_id"]
    reactions, closure, rankings = policies._build_reaction_outputs(selected, registry, families, panels, sets, reaction_map, reference_family)
    if not coverage.complete_required_coverage.all():
        raise ValueError("Incomplete study coverage")
    summary = policies._run_summary(identity, candidates, selected, reactions, sets, families, panels, reaction_map, full_panel_id)
    outputs = {"00_run_summary": summary, "01_input_inventory": identity,
               "02_log_validation": identity[["calculation_id", "method_id", "status", "failures", "warnings", "is_usable"]],
               "03_thermochemistry": identity[["calculation_id", "method_id", "species_id", "electronic_energy_hartree", "thermal_gibbs_correction_hartree", "computed_gibbs_hartree", "printed_gibbs_hartree", "gibbs_residual_hartree"]],
               "04_analysis_candidates": candidates, "05_selected_species": selected,
               "06_sampling_coverage": coverage, "07_reaction_energies": reactions,
               "08_path_closure": closure, "09_relative_rankings": rankings,
               "10_legacy_vs_transferred_audit": legacy_audit(identity)}
    for name, frame in outputs.items():
        frame.to_csv(tables / f"{name}.csv", index=False)
    return outputs
