"""Prepare external sampling and refinement inputs. No calculation is submitted."""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np

from ..records import read_xyz, write_xyz, SYMBOLS, digest
from ..paths import fresh
from ..building.gaussian import input_text
from ..constants import HARTREE_TO_KCAL_MOL
from .selection import read_multi_xyz, radius_of_gyration, select_energy, select_rg_stratified

WATER_JOB_CPUS = 12
WATER_JOB_MEMORY = "40GB"


def sampling_policy(root: Path) -> dict:
    """Per-moiety conformer-selection rules from config/study.json ``sampling``."""
    sampling = json.loads((root / "config/study.json").read_text())["sampling"]
    rules = {}
    for rule in sampling["selection_rules"]:
        window = float(rule["energy_window_kcal_mol"])
        if rule["mode"] == "rg_stratified":
            pol = {"mode": "rg_stratified", "window": window, "n_bins": int(rule["rg_bins"]), "per_bin": int(rule["per_bin"])}
        elif rule["mode"] == "energy":
            pol = {"mode": "energy", "window": window, "n_keep": int(rule["max_keep"])}
        else:
            raise ValueError(f"Unknown sampling selection mode: {rule['mode']}")
        for moiety in rule["moieties"]:
            if moiety in rules:
                raise ValueError(f"Moiety {moiety} has more than one sampling selection rule")
            rules[moiety] = pol
    return {"rules": rules, "rmsd_cutoff_angstrom": float(sampling["heavy_atom_rmsd_cutoff_angstrom"]),
            "include_dft_seed": bool(sampling["include_dft_seed"])}


def prepare_crest(seed: Path, output: Path, charge: int = 0, multiplicity: int = 1, threads: int = 8) -> Path:
    fresh(output)
    z, xyz = read_xyz(seed)
    write_xyz(output / "seed.xyz", [[int(a), *point.tolist()] for a, point in zip(z, xyz)], "unoptimized sampling seed; angstrom")
    command = f"crest seed.xyz --gfn2 --alpb water --ewin 6.0 --T {threads} --chrg {charge} --uhf {multiplicity-1}"
    (output / "run_crest.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\ncd -- \"$(dirname -- \"$0\")\"\n" + command + "\n")
    (output / "sampling_config.json").write_text(json.dumps({"command": command, "seed_sha256": digest(seed), "charge": charge, "multiplicity": multiplicity, "study_threads": threads}, indent=2) + "\n")
    return output


def prepare_refinement(seed: Path, ensemble: Path, moiety: str, output: Path, policy: dict,
                       charge: int = 0, multiplicity: int = 1) -> Path:
    """Select DFT starts from a CREST ensemble; ``policy`` comes from ``sampling_policy``."""
    if moiety not in policy["rules"]:
        raise ValueError(f"Unknown moiety {moiety!r}; configured: {sorted(policy['rules'])}")
    fresh(output)
    z, xyz = read_xyz(seed)
    confs = read_multi_xyz(ensemble)
    if not confs:
        raise ValueError("Empty CREST ensemble")
    symbols = [SYMBOLS[int(a)] for a in z]
    if any(c[1] != symbols or c[2].shape != xyz.shape or not np.isfinite(c[2]).all() for c in confs):
        raise ValueError("Conformer atom order/composition differs from the seed")
    if not np.isfinite([c[0] for c in confs]).all():
        raise ValueError("Nonfinite CREST energies")
    # Preserve the ensemble's frame order and 1-based conformer numbers. Some
    # ensembles have small adjacent energy inversions; sorting would rename
    # candidates and change which frames the selection rules see.
    inversions = sum(confs[i][0] > confs[i+1][0] for i in range(len(confs)-1))
    if confs[0][0] > min(c[0] for c in confs) + 1e-8:
        raise ValueError("First CREST frame must be the energy minimum")
    pol, rmsd_cutoff = policy["rules"][moiety], policy["rmsd_cutoff_angstrom"]
    rel = (np.array([c[0] for c in confs]) - confs[0][0]) * HARTREE_TO_KCAL_MOL
    mask = z != 1
    if not mask.any():
        raise ValueError("Heavy atoms required for conformer selection")
    if pol["mode"] == "energy":
        keep = select_energy(confs, rel, mask, pol, rmsd_cutoff)
    else:
        keep = select_rg_stratified(confs, rel, mask, pol, radius_of_gyration(xyz, mask), rmsd_cutoff)
    seed_start = [("dft_seed", xyz, None)] if policy["include_dft_seed"] else []
    selected = seed_start + [(f"conf{i+1}", confs[i][2], float(rel[i])) for i in keep]
    opt = "# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)"
    freq = "# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read"
    rows = []
    for name, coords, energy in selected:
        write_xyz(output / f"{name}.xyz", [[int(a), *point.tolist()] for a, point in zip(z, coords)], "DFT refinement start; not an optimized result")
        (output / f"{name}.com").write_text(input_text(name, z, coords, route_opt=opt, route_freq=freq, charge=charge, multiplicity=multiplicity))
        rows.append({"candidate": name, "crest_relative_energy_kcal_mol": energy, "heavy_atom_rg_angstrom": radius_of_gyration(coords, mask), "input_sha256": digest(output / f"{name}.com")})
    (output / "refinement_manifest.json").write_text(json.dumps({"policy": pol, "rmsd_cutoff_angstrom": rmsd_cutoff,
        "ensemble_order_policy": "preserve source order and original 1-based conformer numbers", "adjacent_energy_inversions": inversions,
        "seed_sha256": digest(seed), "ensemble_sha256": digest(ensemble), "candidates": rows}, indent=2) + "\n")
    return output


def prepare_waters(xyz_path: Path, family: str, output: Path, folded: bool = False) -> Path:
    from .water import candidate_pool, choose_two
    fresh(output)
    z, _, metrics, candidates = candidate_pool(xyz_path, family, folded=folded)
    acceptors = [value for value in metrics["base_acceptor_atoms"].split(";") if value]
    selected = choose_two(candidates, len(z), require_different_acceptor=len(acceptors) > 1)
    metadata = []
    opt = "# opt=(calcfc,tight,maxcycles=200,maxstep=10) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) scf=(xqc,maxcycle=512)"
    freq = "# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read"
    for label, candidate in zip(["splitA", "splitB"], selected):
        name = family + "_" + label
        write_xyz(output / f"{name}.xyz", [[int(a), *point.tolist()] for a, point in zip(candidate["z"], candidate["xyz"])], "two-water starting placement; unconstrained optimization required")
        (output / f"{name}.com").write_text(input_text(name, candidate["z"], candidate["xyz"], route_opt=opt, route_freq=freq, cpus=WATER_JOB_CPUS, memory=WATER_JOB_MEMORY))
        clean = {k: v for k, v in candidate.items() if k not in ["z", "xyz"]}
        metadata.append({"placement": label, "input_sha256": digest(output / f"{name}.com"), "construction": clean})
    def default(value):
        if hasattr(value, "item"):
            return value.item()
        if hasattr(value, "tolist"):
            return value.tolist()
        raise TypeError(type(value).__name__)
    (output / "placement_manifest.json").write_text(json.dumps({"source_xyz_sha256": digest(xyz_path), "placements": metadata, "atom_indices": "0-based starting-geometry atom order"}, indent=2, default=default) + "\n")
    return output
