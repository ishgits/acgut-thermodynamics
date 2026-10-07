"""Continuum selection policies: lowest-G conformer per species, or one fixed input per species."""
from __future__ import annotations
import pandas as pd
from .reactions import REACTION_ORDER, add_relative_rankings, build_reactions, required_species_for_reactions
from ..constants import HARTREE_TO_KCAL_MOL
from ..geometry import formula_key
NET_NUCLEOTIDE_REACTION = ['P1_total_net_nucleotide_formation']
NUCLEOSIDE_REACTION = ['P1a_nucleoside_formation']


def _allowed_roles(definition: dict[str, str]) -> set[str]:
    """Fixed-input roles for an alternate-method analysis set."""
    roles = {"method_matched_reference"}
    if definition["family_panel"] == "targeted_10":
        roles.add("required_transferred_nucleotide")
    elif definition["family_panel"] == "canonical_5":
        roles.add("transferred_nucleoside")
    return roles

def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.casefold().eq("true")

def _append_failure(frame: pd.DataFrame, mask: pd.Series, reason: str) -> None:
    frame.loc[mask, "failures"] = frame.loc[mask, "failures"].map(
        lambda value: "; ".join(filter(None, [str(value).strip(), reason]))
    )

def _resolve_and_validate(parsed: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    registry = registry.copy()
    registry["formula_key"] = registry["formula"].map(formula_key)
    keys = ["formula_key", "connectivity_fingerprint"]
    ambiguous = registry.groupby(keys, dropna=False)["species_id"].nunique()
    if (ambiguous > 1).any():
        raise ValueError("Registry formula/connectivity keys are not unique")
    identity_columns = keys + [
        "species_id",
        "molecule_class",
        "display_name",
        "base_species_id",
        "canonical",
        "display_order",
    ]
    registry_subset = registry[identity_columns].rename(
        columns={
            "species_id": "resolved_species_id",
            "molecule_class": "resolved_molecule_class",
            "display_name": "resolved_display_name",
            "base_species_id": "resolved_base_species_id",
            "canonical": "resolved_canonical",
            "display_order": "resolved_display_order",
        }
    )
    result = parsed.merge(registry_subset, on=keys, how="left", validate="many_to_one")
    result["failures"] = result["failures"].fillna("")
    result["warnings"] = result["warnings"].fillna("")
    result["hash_matches_manifest"] = result["parsed_sha256"].eq(result["expected_sha256"])
    _append_failure(result, ~result["hash_matches_manifest"], "sha256_mismatch")

    expected_functional = result["expected_functional"].str.upper()
    expected_model = result["expected_solvent_model"].str.upper()
    result["method_matches_manifest"] = (
        result["functional"].eq(expected_functional)
        & result["solvent_model"].eq(expected_model)
        & result["solvent"].eq("WATER")
    )
    _append_failure(result, ~result["method_matches_manifest"], "configured_method_mismatch")

    result["identity_status"] = result["resolved_species_id"].notna().map(
        {True: "resolved", False: "unresolved"}
    )
    manifest_claim = result["manifest_species_id"].ne("")
    mismatch = manifest_claim & result["resolved_species_id"].fillna("").ne(result["manifest_species_id"])
    result["identity_matches_manifest"] = ~mismatch
    _append_failure(result, mismatch, "manifest_identity_mismatch")
    unresolved = result["resolved_species_id"].isna()
    unresolved_claimed = manifest_claim & unresolved
    _append_failure(result, unresolved_claimed, "configured_species_unresolved")
    _append_failure(result, unresolved & ~manifest_claim, "unresolved_identity")

    result.loc[result["failures"].ne(""), "status"] = "failed"
    result["species_id"] = result["resolved_species_id"]
    result["molecule_class"] = result["resolved_molecule_class"]
    result["display_name"] = result["resolved_display_name"]
    result["base_species_id"] = result["resolved_base_species_id"]
    result["canonical"] = result["resolved_canonical"]
    result["display_order"] = pd.to_numeric(result["resolved_display_order"], errors="coerce")
    result["is_usable"] = (
        result["status"].ne("failed")
        & result["species_id"].notna()
        & result["hash_matches_manifest"]
        & result["method_matches_manifest"]
        & result["identity_matches_manifest"]
    )
    return result

def _rank_lowest_candidates(
    rows: pd.DataFrame,
    *,
    analysis_id: str,
    required_species: set[str],
    selection_policy: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = rows.loc[rows["species_id"].isin(required_species) & rows["is_usable"]].copy()
    missing = sorted(required_species - set(candidates["species_id"]))
    if missing:
        raise RuntimeError(f"{analysis_id} is missing usable species: {missing}")
    candidates = candidates.sort_values(
        ["species_id", "computed_gibbs_hartree", "manifest_file_name"]
    ).reset_index(drop=True)
    candidates["analysis_id"] = analysis_id
    candidates["selection_policy"] = selection_policy
    candidates["candidate_rank_within_species"] = candidates.groupby("species_id").cumcount() + 1
    candidates["minimum_gibbs_hartree"] = candidates.groupby("species_id")[
        "computed_gibbs_hartree"
    ].transform("min")
    candidates["delta_g_from_min_kcal_mol"] = (
        candidates["computed_gibbs_hartree"] - candidates["minimum_gibbs_hartree"]
    ) * HARTREE_TO_KCAL_MOL
    candidates["candidate_count"] = candidates.groupby("species_id")["local_path"].transform("size")
    candidates["selected_for_analysis"] = candidates["candidate_rank_within_species"].eq(1)
    selected = candidates.loc[candidates["selected_for_analysis"]].copy()
    return candidates, selected

def _select_fixed_inputs(
    rows: pd.DataFrame,
    *,
    analysis_id: str,
    required_species: set[str],
    allowed_roles: set[str],
    selection_policy: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = rows.loc[
        rows["is_usable"]
        & rows["primary_include"]
        & rows["analysis_role"].isin(allowed_roles)
        & rows["species_id"].isin(required_species)
    ].copy()
    counts = candidates.groupby("species_id").size()
    missing = sorted(required_species - set(counts.index))
    duplicated = sorted(counts.loc[counts.ne(1)].index.tolist())
    if missing or duplicated:
        raise RuntimeError(
            f"{analysis_id} fixed-input violation; missing={missing}, nonunique={duplicated}"
        )
    candidates = candidates.sort_values("species_id").reset_index(drop=True)
    candidates["analysis_id"] = analysis_id
    candidates["selection_policy"] = selection_policy
    candidates["candidate_rank_within_species"] = 1
    candidates["minimum_gibbs_hartree"] = candidates["computed_gibbs_hartree"]
    candidates["delta_g_from_min_kcal_mol"] = 0.0
    candidates["candidate_count"] = 1
    candidates["selected_for_analysis"] = True
    return candidates, candidates.copy()

def _analysis_inputs(
    identity: pd.DataFrame,
    analysis_sets: pd.DataFrame,
    registry: pd.DataFrame,
    families: pd.DataFrame,
    panels: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    reaction_ids_by_scope = {
        "all balanced reaction quantities": list(REACTION_ORDER),
        "net nucleotide formation": NET_NUCLEOTIDE_REACTION,
        "nucleoside formation thermochemistry": NUCLEOSIDE_REACTION,
    }
    candidate_parts: list[pd.DataFrame] = []
    selected_parts: list[pd.DataFrame] = []
    reaction_map: dict[str, list[str]] = {}
    for definition in analysis_sets.to_dict("records"):
        analysis_id = definition["analysis_id"]
        method_id = definition["method_id"]
        panel_families = families.loc[
            families["reaction_family_id"].isin(panels[definition["family_panel"]])
        ].copy()
        reaction_ids = reaction_ids_by_scope[definition["reaction_scope"]]
        reaction_map[analysis_id] = reaction_ids
        required = required_species_for_reactions(panel_families, reaction_ids)
        method_rows = identity.loc[identity["method_id"].eq(method_id)]
        if method_id == "1_b3lyp_pcm":
            pool = method_rows.loc[
                method_rows["analysis_role"].eq("conformer_candidate")
                & method_rows["primary_include"]
            ]
            candidates, selected = _rank_lowest_candidates(
                pool,
                analysis_id=analysis_id,
                required_species=required,
                selection_policy=definition["product_policy"],
            )
        else:
            candidates, selected = _select_fixed_inputs(
                method_rows,
                analysis_id=analysis_id,
                required_species=required,
                allowed_roles=_allowed_roles(definition),
                selection_policy=definition["product_policy"],
            )
        for column in ["family_panel", "reaction_scope", "manuscript_status"]:
            candidates[column] = definition[column]
            selected[column] = definition[column]
        candidate_parts.append(candidates)
        selected_parts.append(selected)
    return (
        pd.concat(candidate_parts, ignore_index=True),
        pd.concat(selected_parts, ignore_index=True),
        reaction_map,
    )

def _build_coverage(
    candidates: pd.DataFrame,
    selected: pd.DataFrame,
    analysis_sets: pd.DataFrame,
    families: pd.DataFrame,
    panels: dict[str, list[str]],
    reaction_map: dict[str, list[str]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for definition in analysis_sets.to_dict("records"):
        analysis_id = definition["analysis_id"]
        scoped_families = families.loc[
            families["reaction_family_id"].isin(panels[definition["family_panel"]])
        ].sort_values("display_order")
        candidate_counts = candidates.loc[candidates["analysis_id"].eq(analysis_id)].groupby("species_id").size()
        selected_species = set(selected.loc[selected["analysis_id"].eq(analysis_id), "species_id"])
        required = required_species_for_reactions(scoped_families, reaction_map[analysis_id])
        shared = {
            species: int(candidate_counts.get(species, 0))
            for species in ("ribose", "phosphate", "water", "ribo_phosphate")
        }
        for family in scoped_families.to_dict("records"):
            family_species = {
                "base": family["base_species_id"],
                "nucleoside": family["nucleoside_species_id"],
                "nucleotide": family["nucleotide_species_id"],
            }
            family_required = {
                species for species in family_species.values() if species in required
            }
            family_required.update(species for species in shared if species in required)
            rows.append(
                {
                    "analysis_id": analysis_id,
                    "method_id": definition["method_id"],
                    "method_label": definition["method_label"],
                    "family_panel": definition["family_panel"],
                    "reaction_family_id": family["reaction_family_id"],
                    "family_display_name": family["display_name"],
                    "base_candidate_count": int(candidate_counts.get(family_species["base"], 0)),
                    "nucleoside_candidate_count": int(candidate_counts.get(family_species["nucleoside"], 0)),
                    "nucleotide_candidate_count": int(candidate_counts.get(family_species["nucleotide"], 0)),
                    "ribose_candidate_count": shared["ribose"],
                    "phosphate_candidate_count": shared["phosphate"],
                    "water_candidate_count": shared["water"],
                    "ribo_phosphate_candidate_count": shared["ribo_phosphate"],
                    "required_species_count": len(family_required),
                    "complete_required_coverage": family_required.issubset(selected_species),
                }
            )
    return pd.DataFrame(rows)

def _build_reaction_outputs(
    selected: pd.DataFrame,
    registry: pd.DataFrame,
    families: pd.DataFrame,
    panels: dict[str, list[str]],
    analysis_sets: pd.DataFrame,
    reaction_map: dict[str, list[str]],
    reference_family: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reaction_parts: list[pd.DataFrame] = []
    closure_parts: list[pd.DataFrame] = []
    ranking_parts: list[pd.DataFrame] = []
    for definition in analysis_sets.to_dict("records"):
        analysis_id = definition["analysis_id"]
        chosen = selected.loc[selected["analysis_id"].eq(analysis_id)].copy()
        scoped_families = families.loc[
            families["reaction_family_id"].isin(panels[definition["family_panel"]])
        ].copy()
        reactions, closure = build_reactions(
            chosen,
            registry,
            scoped_families,
            reaction_ids=reaction_map[analysis_id],
        )
        if not reactions["balanced"].all() or reactions["missing_species"].ne("").any():
            raise RuntimeError(f"Reaction construction failed for {analysis_id}")
        metadata = {
            "analysis_id": analysis_id,
            "method_id": definition["method_id"],
            "method_label": definition["method_label"],
            "family_panel": definition["family_panel"],
            "reaction_scope": definition["reaction_scope"],
            "manuscript_status": definition["manuscript_status"],
        }
        for key, value in metadata.items():
            reactions[key] = value
            closure[key] = value
        if not closure.empty and not closure["passes"].all():
            raise RuntimeError(f"Reaction path closure failed for {analysis_id}")
        rankings = add_relative_rankings(
            reactions,
            panel_name=definition["family_panel"],
            reference_family=reference_family,
        )
        reaction_parts.append(reactions)
        if not closure.empty:
            closure_parts.append(closure)
        ranking_parts.append(rankings)
    return (
        pd.concat(reaction_parts, ignore_index=True),
        pd.concat(closure_parts, ignore_index=True) if closure_parts else pd.DataFrame(),
        pd.concat(ranking_parts, ignore_index=True),
    )

def _relevant_identity_rows(
    identity: pd.DataFrame,
    definition: dict[str, str],
    families: pd.DataFrame,
    panels: dict[str, list[str]],
    reaction_ids: list[str],
    full_panel_id: str,
) -> pd.DataFrame:
    scoped = families.loc[families["reaction_family_id"].isin(panels[definition["family_panel"]])]
    required = required_species_for_reactions(scoped, reaction_ids)
    rows = identity.loc[identity["method_id"].eq(definition["method_id"])]
    if definition["method_id"] == "1_b3lyp_pcm":
        include_unresolved_provenance = definition["analysis_id"] == full_panel_id
        return rows.loc[
            rows["analysis_role"].eq("conformer_candidate")
            & rows["primary_include"]
            & (
                rows["species_id"].isin(required)
                | (include_unresolved_provenance & rows["species_id"].isna())
            )
        ]
    return rows.loc[
        rows["primary_include"]
        & rows["analysis_role"].isin(_allowed_roles(definition))
        & rows["species_id"].isin(required)
    ]

def _run_summary(
    identity: pd.DataFrame,
    candidates: pd.DataFrame,
    selected: pd.DataFrame,
    reactions: pd.DataFrame,
    analysis_sets: pd.DataFrame,
    families: pd.DataFrame,
    panels: dict[str, list[str]],
    reaction_map: dict[str, list[str]],
    full_panel_id: str,
) -> pd.DataFrame:
    """One row per analysis set; the full panel also counts records of unresolved identity."""
    rows: list[dict[str, object]] = []
    for definition in analysis_sets.to_dict("records"):
        analysis_id = definition["analysis_id"]
        relevant = _relevant_identity_rows(identity, definition, families, panels, reaction_map[analysis_id], full_panel_id)
        chosen = selected.loc[selected["analysis_id"].eq(analysis_id)]
        analysis_reactions = reactions.loc[reactions["analysis_id"].eq(analysis_id)]
        rows.append(
            {
                **definition,
                "configured_log_count": len(relevant),
                "usable_log_count": int(relevant["is_usable"].sum()),
                "excluded_log_count": int((~relevant["is_usable"]).sum()),
                "candidate_row_count": len(candidates.loc[candidates["analysis_id"].eq(analysis_id)]),
                "selected_species_count": len(chosen),
                "reaction_quantity_count": len(analysis_reactions),
                "run_status": "passed_with_exclusions" if (~relevant["is_usable"]).any() else "passed",
            }
        )
    return pd.DataFrame(rows)
