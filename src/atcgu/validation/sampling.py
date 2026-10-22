"""Replay released preparation candidate sets without generating jobs."""
import csv
import hashlib
import json
import numpy as np
import pandas as pd
from ..sampling import selection as s
from ..constants import HARTREE_TO_KCAL_MOL


def check_run(root, row, study):
    rules = {}
    for rule in study["sampling"]["selection_rules"]:
        pol = {"mode": rule["mode"], "window": rule["energy_window_kcal_mol"]}
        if pol["mode"] == "energy":
            pol["n_keep"] = rule["max_keep"]
        else:
            pol.update(n_bins=rule["rg_bins"], per_bin=rule["per_bin"])
        for moiety in rule["moieties"]:
            rules[moiety] = pol
    cutoff = study["sampling"]["heavy_atom_rmsd_cutoff_angstrom"]
    result = {"sampling_run_id": row["sampling_run_id"]}
    seed, ensemble = root / row["seed_path"], root / row["ensemble_path"]
    lines = seed.read_text().splitlines()
    atoms = [line.split() for line in lines[2:2 + int(lines[0])]]
    if len(atoms) != int(lines[0]) or any(len(a) != 4 for a in atoms):
        raise ValueError("Malformed seed atom count/coordinates")
    symbols = [a[0] for a in atoms]
    xyz = np.array([[float(v) for v in a[1:4]] for a in atoms])
    mask = np.array([a != "H" for a in symbols])
    confs = s.read_multi_xyz(ensemble)
    if not confs or not mask.any() or xyz.shape != (len(symbols), 3) or not np.isfinite(xyz).all():
        raise ValueError("Empty ensemble, missing heavy atoms or invalid seed")
    if any(c[1] != symbols or c[2].shape != xyz.shape or not np.isfinite(c[2]).all() or not np.isfinite(c[0]) for c in confs):
        raise ValueError("Nonfinite energy/coordinates or inconsistent atom order/count")
    if confs[0][0] - min(c[0] for c in confs) > 1e-8:
        raise ValueError("First frame is not the minimum within 1e-8 hartree")
    relative = (np.array([c[0] for c in confs]) - confs[0][0]) * HARTREE_TO_KCAL_MOL
    pol = rules[row["moiety"]]
    def select(rel):
        if pol["mode"] == "energy":
            return s.select_energy(confs, rel, mask, pol, cutoff)
        return s.select_rg_stratified(confs, rel, mask, pol, s.radius_of_gyration(xyz, mask), cutoff)
    keep = select(relative)
    ids = ["dft_seed"] + [f"conf{k+1}" for k in keep]
    # Preserve sampling run ID casing; it differs from species IDs for linkage families.
    mpath = root / "calculation_setup/crest/refinement" / row["sampling_run_id"] / "refinement_manifest.json"
    manifest = json.loads(mpath.read_text())
    saved = [c["candidate"] for c in manifest["candidates"]]
    for candidate in manifest["candidates"]:
        if not np.isfinite(candidate["heavy_atom_rg_angstrom"]):
            raise ValueError("Nonfinite saved Rg")
        energy = candidate["crest_relative_energy_kcal_mol"]
        if candidate["candidate"] != "dft_seed" and (energy is None or not np.isfinite(energy)):
            raise ValueError("Nonfinite saved CREST energy")
    e_errors, rg_errors = [], []
    if ids == saved:
        rg_errors.append(abs(s.radius_of_gyration(xyz, mask) - manifest["candidates"][0]["heavy_atom_rg_angstrom"]))
        for k, c in zip(keep, manifest["candidates"][1:]):
            e_errors.append(abs(relative[k] - c["crest_relative_energy_kcal_mol"]))
            rg_errors.append(abs(s.radius_of_gyration(confs[k][2], mask) - c["heavy_atom_rg_angstrom"]))
    hashes_match = all(hashlib.sha256(path.read_bytes()).hexdigest() == manifest[key] == row[key] for path, key in [(seed, "seed_sha256"), (ensemble, "ensemble_sha256")])
    result.update({"sampling_run_id": row["sampling_run_id"], "policy_matches": pol == manifest["policy"], "hashes_match": hashes_match,
                 "candidates_match": ids == saved, "selected": json.dumps(ids), "published": json.dumps(saved),
                 "max_energy_error_kcal_mol": max(e_errors, default=0),
                 "max_rg_error_angstrom": max(rg_errors, default=0),
                 "inversions_match": sum(confs[i][0] > confs[i+1][0] for i in range(len(confs)-1)) == manifest["adjacent_energy_inversions"],
                 "inversion_count": sum(confs[i][0] > confs[i+1][0] for i in range(len(confs)-1)),
                 "expected_inversion_count": manifest["adjacent_energy_inversions"],
                 "policy": json.dumps(pol, sort_keys=True),
                 "expected_policy": json.dumps(manifest["policy"], sort_keys=True),
                 "seed_energy_null": manifest["candidates"][0]["crest_relative_energy_kcal_mol"] is None,
    })
    failures = [key for key in ["policy_matches", "hashes_match", "candidates_match", "inversions_match", "seed_energy_null"] if not result[key]]
    if result["max_energy_error_kcal_mol"] > 1e-9:
        failures.append("energy differs by more than 1e-9 kcal/mol")
    if result["max_rg_error_angstrom"] > 1e-9:
        failures.append("Rg differs by more than 1e-9 angstrom")
    result.update(passed=not failures, reason="; ".join(failures))
    return result


def validate_sampling(root, output):
    study = json.loads((root / "config/study.json").read_text())
    results = []
    for row in csv.DictReader((root / "config/sampling.csv").open()):
        try:
            result = check_run(root, row, study)
        except (ValueError, KeyError, IndexError, OSError, TypeError) as exc:
            result = {"sampling_run_id": row["sampling_run_id"], "passed": False, "reason": f"{type(exc).__name__}: {exc}"}
        results.append(result)
    pd.DataFrame(results).to_csv(output / "sampling_selection_checks.csv", index=False)
    return [{"item": r["sampling_run_id"], "passed": r["passed"], "status": "sampling_selection", "details": r["reason"]} for r in results]
