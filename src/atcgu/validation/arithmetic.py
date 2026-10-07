"""Independent arithmetic audit. Standard library only; no production analysis imports.

Re-derive selection, balanced stoichiometry, reaction energies, references, ranks,
hydration minima and plotted coordinates from the released records and registries,
and compare them with the published results. The SI section recounts selections,
orders and the 64 placement vectors from the published treated energies; the
entropy formulas themselves stay in production code.
"""
from pathlib import Path
import csv
import json
import math
import re
from collections import defaultdict, Counter
from itertools import combinations, product


def audit(root: Path, output: Path) -> Path:
    root, output = root.resolve(), output.resolve()
    if output.is_relative_to(root) and not output.is_relative_to(root / "outputs"):
        raise ValueError(f"Inside the checkout, write outputs under {root / 'outputs'}")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose a fresh audit directory")
    output.mkdir(parents=True, exist_ok=True)
    def read(name):
        with (root / name).open(newline="") as h:
            return list(csv.DictReader(h))
    study = json.loads((root / "config/study.json").read_text())
    factor = study["hartree_to_kcal_mol"]
    records = {}
    for row in read("data/calculations/manifest.csv"):
        r = json.loads((root / row["record_path"]).read_text())
        if r["status"] != "failed":
            r["independent_g"] = r["electronic_energy_hartree"] + r["thermal_gibbs_correction_hartree"]
            if not math.isfinite(r["independent_g"]):
                raise ValueError("Nonfinite energy")
        records[r["calculation_id"]] = r
    manifest = read("config/calculation_manifest.csv")
    families = read("config/reaction_families.csv")
    panels = {"all_21": [r["reaction_family_id"] for r in families], **study["family_panels"]}
    formulas = {r["species_id"]: {e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", r["formula"])} for r in read("config/molecules.csv")}
    charges = {r["species_id"]: int(r["charge"]) for r in read("config/molecules.csv")}
    stoich = defaultdict(dict)
    for row in read("config/reaction_stoichiometry.csv"):
        key = (row["reaction_family_id"], row["reaction_id"])
        if row["species_id"] in stoich[key]:
            raise ValueError("Duplicate stoichiometric species")
        stoich[key][row["species_id"]] = int(row["coefficient"])
    scopes = {"net nucleotide formation": ["P1_total_net_nucleotide_formation"],
              "nucleoside formation thermochemistry": ["P1a_nucleoside_formation"],
              "all balanced reaction quantities": sorted({k[1] for k in stoich})}
    checks = []
    def check(label, actual, expected, tolerance=1e-8):
        numeric = isinstance(actual, (float, int)) and not isinstance(actual, bool)
        error = abs(actual - float(expected)) if numeric else None
        passed = math.isfinite(error) and error <= tolerance if numeric else actual == expected
        checks.append({"check_id": label, "passed": bool(passed), "actual": actual,
                       "expected": expected, "absolute_error": error, "tolerance": tolerance if numeric else ""})
    published = {(r["analysis_id"], r["species_id"]): r for r in read("data/published/continuum/05_selected_species.csv")}
    published_rxn = {(r["analysis_id"], r["reaction_family_id"], r["reaction_id"]): r for r in read("data/published/continuum/07_reaction_energies.csv")}
    published_rank = {(r["analysis_id"], r["reaction_family_id"], r["reaction_id"]): r for r in read("data/published/continuum/09_relative_rankings.csv")}
    selected = {}; deltas = {}; relative = {}
    for definition in read("config/analysis_sets.csv"):
        aid, method = definition["analysis_id"], definition["method_id"]
        keys = [(f, r) for f in panels[definition["family_panel"]] for r in scopes[definition["reaction_scope"]]]
        required = set().union(*(set(stoich[k]) for k in keys))
        choices = {}
        roles = {"method_matched_reference", "required_transferred_nucleotide" if definition["family_panel"] == "targeted_10" else "transferred_nucleoside"}
        for sid in sorted(required):
            candidates = [r for r in manifest if r["method_id"] == method and r["species_id"] == sid
                          and r["primary_include"].lower() == "true" and records[r["calculation_id"]]["status"] != "failed"
                          and r["analysis_role"] in ({"conformer_candidate"} if method == "1_b3lyp_pcm" else roles)]
            if not candidates or (method != "1_b3lyp_pcm" and len(candidates) != 1):
                raise ValueError(f"Independent selection coverage failure: {aid}:{sid}")
            choice = min(candidates, key=lambda r: (records[r["calculation_id"]]["independent_g"], r["file_name"]))
            choices[sid] = records[choice["calculation_id"]]["independent_g"]
            selected[aid, sid] = choice["calculation_id"]
            golden = published[aid, sid]
            check(f"selection:{aid}:{sid}", choice["calculation_id"], golden["calculation_id"])
            check(f"selection_energy:{aid}:{sid}", choices[sid], golden["computed_gibbs_hartree"], 1e-10)
            check(f"candidate_count:{aid}:{sid}", len(candidates), golden["candidate_count"], 0)
        for family, reaction in keys:
            coeffs = stoich[family, reaction]
            residual = Counter()
            for sid, c in coeffs.items():
                for e, n in formulas[sid].items(): residual[e] += c * n
            check(f"balance:{aid}:{family}:{reaction}", all(n == 0 for n in residual.values()) and sum(c * charges[s] for s, c in coeffs.items()) == 0, True)
            delta = math.fsum(c * choices[sid] for sid, c in coeffs.items()) * factor
            deltas[aid, family, reaction] = delta
            check(f"reaction:{aid}:{family}:{reaction}", delta, published_rxn[aid, family, reaction]["delta_g_kcal_mol"])
        for family, reaction in keys:
            value = deltas[aid, family, reaction]
            rel = value - deltas[aid, study["reference_family"], reaction]
            relative[aid, family, reaction] = rel
            rank = 1 + sum(deltas[aid, f, reaction] < value for f in panels[definition["family_panel"]])
            golden = published_rank[aid, family, reaction]
            check(f"relative:{aid}:{family}:{reaction}", rel, golden["relative_g_kcal_mol"])
            check(f"rank:{aid}:{family}:{reaction}", rank, golden["rank_within_panel"], 0)
    if len(selected) != len(published) or len(deltas) != len(published_rxn):
        raise ValueError("Independent cohort size disagrees with the published results")
    refs = {(r["family"], r["role"]): records[r["calculation_id"]]["independent_g"] for r in read("config/microsolvation_references.csv")}
    primary = read("config/microsolvation_candidates.csv")
    scores = {r["candidate_id"]: records[r["calculation_id"]]["independent_g"] - refs[r["family"], "base"] for r in primary}
    best = {family: min([r for r in primary if r["family"] == family], key=lambda r: (scores[r["candidate_id"]], r["candidate_id"])) for family in {r["family"] for r in primary}}
    for golden in read("data/published/figure1/06_primary_selected_candidates.csv"):
        f = golden["family"]; winner = best[f]["candidate_id"]
        check(f"hydration_selection:{f}", winner, golden["candidate_id"])
        rel = (scores[winner] - scores[best["adenine"]["candidate_id"]]) * factor
        check(f"hydration_relative:{f}", rel, golden["relative_g_kcal_mol"])
    conv = study["panel_c_display_convention"]
    offset = -2 * conv["gas_constant_kcal_mol_k"] * conv["temperature_k"] * (math.log(conv["gas_to_solution_factor"]) + math.log(conv["water_reservoir_molarity"]))
    full = "full21_b3lyp_pcm"; water = records[selected[full, "water"]]["independent_g"]
    for row in read("data/published/figure1/figure1_plotted_data.csv"):
        panel, f, series = row["panel"], row["family"], row["series"]
        if panel in {"A", "B"}:
            aid = full if panel == "A" else "method_sensitivity10_" + series
            value = relative[aid, f, "P1_total_net_nucleotide_formation"]
        else:
            f = "tap" if f == "tap_c5link" else f
            # Derive relative dry product-base score directly, as in the display convention.
            dry = refs[f, "continuum_nucleotide"] - refs[f, "base"]
            zero = refs["adenine", "continuum_nucleotide"] - refs["adenine", "base"]
            value = (dry - zero) * factor
            if series != "pcm_only":
                placement = series.rsplit("_", 1)[1]
                row_c = next(r for r in primary if r["family"] == f and r["placement"] == placement)
                cluster = records[row_c["calculation_id"]]["independent_g"]
                value = (cluster - refs[f, "base"] - zero - 2 * water) * conv["hartree_to_kcal_mol"] + offset
        check(f"figure:{panel}:{f}:{series}", value, row["value"], 1e-6)
    figure_points = len(read("data/published/figure1/figure1_plotted_data.csv"))
    # SI: recount from the published treated energies, keyed by the configured identities.
    si = study["si_microsolvation"]
    order, sign = si["family_order"], lambda x: (x > 0) - (x < 0)
    treated = {(r["treatment"], r["calculation_id"]): r for r in read("data/published/si_microsolvation/03_treatment_energies.csv")}
    ref_ids = {(r["family"], r["role"]): r["calculation_id"] for r in read("config/microsolvation_references.csv")}
    for (t, cid), r in treated.items():
        if t == "harmonic":
            check(f"si_harmonic_record:{cid}", float(r["treated_g_hartree"]), records[cid]["independent_g"], 1e-10)
        check(f"si_shift_sum:{t}:{cid}", float(r["harmonic_g_hartree"]) + float(r["entropy_shift_hartree"]), r["treated_g_hartree"], 1e-10)
    published_vectors = {(r["treatment"], r["choice_bits"]): r for r in read("data/published/si_microsolvation/09_placement_combinations.csv")}
    published_placements = {(r["treatment"], r["family"], r["placement"]): r for r in read("data/published/si_microsolvation/05_fixed_reference_placement_scores.csv")}
    published_shifts = {(r["treatment"], r["family"]): r for r in read("data/published/si_microsolvation/08_selected_fixed_reference_shifts.csv")}
    published_occupancy = {(r["treatment"], r["family"], int(r["rank"])): r for r in read("data/published/si_microsolvation/10_rank_occupancy.csv")}
    published_si = {r["treatment"]: r for r in read("data/published/si_microsolvation/thermochemical_summary.csv")}
    vectors = 0
    for t in si["treatments"]:
        g = lambda cid: float(treated[t, cid]["treated_g_hartree"])
        score = {r["candidate_id"]: g(r["calculation_id"]) - g(ref_ids[r["family"], "base"]) for r in primary}
        by_place = {(r["family"], r["placement"]): r["candidate_id"] for r in primary}
        dry = {f: g(ref_ids[f, "continuum_nucleotide"]) - g(ref_ids[f, "base"]) for f in order}
        chosen = {f: min((c for (fam, _), c in by_place.items() if fam == f), key=lambda c: (score[c], c)) for f in order}
        best = {f: score[chosen[f]] for f in order}
        fixed_score = {candidate: (value - zero - 2 * water) * conv["hartree_to_kcal_mol"] + offset
                       for candidate, value in score.items()}
        for row in primary:
            golden = published_placements[t, row["family"], row["placement"]]
            check(f"si_fixed_score:{t}:{row['candidate_id']}", fixed_score[row["candidate_id"]], golden["fixed_reference_score_kcal_mol"], 1e-6)
            check(f"si_fixed_selected:{t}:{row['candidate_id']}", row["candidate_id"] == chosen[row["family"]], golden["selected_for_family"].lower() == "true")
        kept = lambda vals: sum(sign(dry[a] - dry[b]) == sign(vals[a] - vals[b]) for a, b in combinations(order, 2))
        orders, ac = set(), 0
        occupancy = Counter()
        for bits in product([0, 1], repeat=len(order)):
            vals = {f: score[by_place[f, si["placement_bits"][b]]] for f, b in zip(order, bits)}
            ranked = sorted(order, key=vals.get)
            for rank, family in enumerate(ranked, start=1): occupancy[family, rank] += 1
            orders.add(tuple(ranked)); ac += set(ranked[:2]) == {"adenine", "cytosine"}
            golden = published_vectors[t, "".join(map(str, bits))]
            check(f"si_vector_order:{t}:{''.join(map(str, bits))}", " < ".join(ranked), golden["rank_order"])
            check(f"si_vector_c_minus_a:{t}:{''.join(map(str, bits))}", (vals["cytosine"] - vals["adenine"]) * factor, golden["c_minus_a_kcal_mol"], 1e-6)
            check(f"si_vector_pairs:{t}:{''.join(map(str, bits))}", kept(vals), golden["pairwise_retained"], 0)
            vectors += 1
        golden = published_si[t]
        check(f"si_selected_order:{t}", " < ".join(sorted(order, key=lambda f: (best[f], f))), golden["selected_order"])
        check(f"si_c_minus_a:{t}", (best["cytosine"] - best["adenine"]) * factor, golden["c_minus_a"], 1e-9)
        check(f"si_pairwise_retained:{t}", kept(best), golden["pairwise_retained"], 0)
        check(f"si_unique_orders:{t}", len(orders), golden["unique_placement_orders"], 0)
        check(f"si_ac_top_two:{t}", ac, golden["ac_top_two_choices"], 0)
        selected_order = sorted(order, key=lambda f: (best[f], f))
        for family in order:
            selected_rank = selected_order.index(family) + 1
            fixed = fixed_score[chosen[family]]
            shift_golden = published_shifts[t, family]
            harmonic_golden = published_shifts["harmonic", family]
            check(f"si_selected_fixed:{t}:{family}", fixed, shift_golden["selected_fixed_reference_score_kcal_mol"], 1e-6)
            check(f"si_selected_shift:{t}:{family}", fixed - float(harmonic_golden["selected_fixed_reference_score_kcal_mol"]),
                  shift_golden["selected_fixed_reference_shift_from_harmonic_kcal_mol"], 1e-6)
            for rank in range(1, len(order) + 1):
                occ_golden = published_occupancy[t, family, rank]
                check(f"si_rank_occupancy:{t}:{family}:{rank}", occupancy[family, rank], occ_golden["combination_count"], 0)
                check(f"si_selected_rank:{t}:{family}:{rank}", selected_rank, occ_golden["selected_lower_g_rank"], 0)
    if vectors != 64 * len(si["treatments"]) or len(published_vectors) != vectors:
        raise ValueError("Independent SI vector count disagrees with the published results")
    with (output / "independent_arithmetic_checks.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(checks[0])); w.writeheader(); w.writerows(checks)
    summary = {"checks": len(checks), "failed_checks": sum(not r["passed"] for r in checks), "selected_species": len(selected), "reaction_quantities": len(deltas), "figure_points": figure_points, "si_placement_vectors": vectors,
               "maximum_numeric_error": max(r["absolute_error"] or 0 for r in checks),
               "scope": "arithmetic, cohort selection, published Figure 1, fixed-reference SI scores, placement vectors and rank occupancy"}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if summary["failed_checks"]:
        raise ValueError("Independent arithmetic checks failed; inspect the audit CSV")
    return output


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(); p.add_argument("root", type=Path); p.add_argument("output", type=Path)
    args = p.parse_args(); audit(args.root.resolve(), args.output.resolve())
