"""Primary two-water split-site starts: one water accepts a phosphate OH, one donates to a base acceptor."""
from pathlib import Path
import numpy as np
import math
from ..records import read_xyz
from ..constants import VDW_RADII as VDW
from ..geometry import bond_graph, hbond_geometry, is_hbond
from ..structure import analyze_structure

def full_water_clearance(z, xyz, solute_count, solute_edges, intended_contacts):
    """Return the closest unintended contact involving an explicit-water atom."""
    excluded = {frozenset(edge) for edge in solute_edges}
    excluded.update({
        frozenset((solute_count, solute_count + 1)),
        frozenset((solute_count, solute_count + 2)),
        frozenset((solute_count + 1, solute_count + 2)),
        frozenset((solute_count + 3, solute_count + 4)),
        frozenset((solute_count + 3, solute_count + 5)),
        frozenset((solute_count + 4, solute_count + 5)),
    })
    excluded.update(frozenset(pair) for pair in intended_contacts)
    contacts = []
    for i in range(len(z)):
        for j in range(i + 1, len(z)):
            if i < solute_count and j < solute_count:
                continue
            if frozenset((i, j)) in excluded:
                continue
            distance = np.linalg.norm(xyz[i] - xyz[j])
            contacts.append((distance / (VDW[int(z[i])] + VDW[int(z[j])]), distance, i, j))
    return min(contacts)

def unit(vector):
    norm = np.linalg.norm(vector)
    if norm < 1e-10:
        raise ValueError("cannot normalize zero vector")
    return vector / norm

def orthogonal_basis(axis):
    axis = unit(axis)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(axis @ trial) > 0.85:
        trial = np.array([0.0, 1.0, 0.0])
    first = unit(trial - (trial @ axis) * axis)
    return first, unit(np.cross(axis, first))

def rotate_point(point, origin, axis, degrees):
    axis = unit(axis)
    angle = math.radians(degrees)
    vector = point - origin
    return origin + (
        vector * math.cos(angle)
        + np.cross(axis, vector) * math.sin(angle)
        + axis * (axis @ vector) * (1 - math.cos(angle))
    )

def one_water_clearance(solute_z, solute_xyz, water, intended_contact):
    full_z = np.concatenate([solute_z, [8, 1, 1]])
    full_xyz = np.vstack([solute_xyz, water])
    n = len(solute_z)
    excluded = {
        frozenset((n, n + 1)),
        frozenset((n, n + 2)),
        frozenset((n + 1, n + 2)),
        frozenset(intended_contact),
    }
    contacts = []
    for i in range(n, n + 3):
        for j in range(i):
            if frozenset((i, j)) in excluded:
                continue
            distance = np.linalg.norm(full_xyz[i] - full_xyz[j])
            contacts.append(distance / (VDW[int(full_z[i])] + VDW[int(full_z[j])]))
    return min(contacts) if contacts else 10.0

def phosphate_water_candidates(z, xyz, donor_o, donor_h, p_atom):
    output = []
    oh_length = np.linalg.norm(xyz[donor_h] - xyz[donor_o])
    bond_axis = xyz[p_atom] - xyz[donor_o]
    for chi in range(0, 360, 10):
        trial = xyz.copy()
        trial[donor_h] = rotate_point(xyz[donor_h], xyz[donor_o], bond_axis, chi)
        direction = unit(trial[donor_h] - trial[donor_o])
        water_o = trial[donor_o] + (oh_length + 1.75) * direction
        away = direction
        u, v = orthogonal_basis(away)
        half = math.radians(104.5 / 2)
        for psi in range(0, 360, 15):
            radial = math.cos(math.radians(psi)) * u + math.sin(math.radians(psi)) * v
            h1 = water_o + 0.96 * (math.cos(half) * away + math.sin(half) * radial)
            h2 = water_o + 0.96 * (math.cos(half) * away - math.sin(half) * radial)
            water = np.array([water_o, h1, h2])
            clearance = one_water_clearance(z, trial, water, (donor_h, len(z)))
            geometry = hbond_geometry(donor_o, donor_h, len(z), np.vstack([trial, water]))
            if clearance >= 0.70 and is_hbond(geometry):
                output.append({
                    "solute_xyz": trial,
                    "water": water,
                    "chi": chi,
                    "psi": psi,
                    "clearance": clearance,
                    "geometry": geometry,
                })
    return sorted(output, key=lambda item: item["clearance"], reverse=True)[:24]

def acceptor_directions(z, xyz, acceptor, adjacency):
    heavy_neighbors = [atom for atom in adjacency[acceptor] if z[atom] != 1]
    ideal = -sum((unit(xyz[atom] - xyz[acceptor]) for atom in heavy_neighbors), np.zeros(3))
    if np.linalg.norm(ideal) < 0.2:
        heavy_points = xyz[[i for i, number in enumerate(z) if number != 1]]
        ideal = xyz[acceptor] - heavy_points.mean(axis=0)
    ideal = unit(ideal)
    u, v = orthogonal_basis(ideal)
    directions = [ideal]
    for cone_angle in (10, 20, 30, 40, 50):
        alpha = math.radians(cone_angle)
        for phi in range(0, 360, 20):
            radial = math.cos(math.radians(phi)) * u + math.sin(math.radians(phi)) * v
            directions.append(unit(math.cos(alpha) * ideal + math.sin(alpha) * radial))
    return directions

def base_water_candidates(z, xyz, acceptor, adjacency):
    output = []
    for direction_id, direction in enumerate(acceptor_directions(z, xyz, acceptor, adjacency)):
        water_o = xyz[acceptor] + 2.80 * direction
        h1_direction = -direction
        h1 = water_o + 0.96 * h1_direction
        u, v = orthogonal_basis(h1_direction)
        theta = math.radians(104.5)
        for psi in range(0, 360, 15):
            radial = math.cos(math.radians(psi)) * u + math.sin(math.radians(psi)) * v
            h2 = water_o + 0.96 * (math.cos(theta) * h1_direction + math.sin(theta) * radial)
            water = np.array([water_o, h1, h2])
            clearance = one_water_clearance(z, xyz, water, (acceptor, len(z) + 1))
            geometry = hbond_geometry(len(z), len(z) + 1, acceptor, np.vstack([xyz, water]))
            if clearance >= 0.70 and is_hbond(geometry):
                output.append({
                    "water": water,
                    "direction_id": direction_id,
                    "psi": psi,
                    "clearance": clearance,
                    "geometry": geometry,
                })
    best_by_direction = {}
    for candidate in output:
        direction_id = candidate["direction_id"]
        if direction_id not in best_by_direction or candidate["clearance"] > best_by_direction[direction_id]["clearance"]:
            best_by_direction[direction_id] = candidate
    return sorted(best_by_direction.values(), key=lambda item: item["clearance"], reverse=True)[:36]

def candidate_pool(xyz_path, family, *, folded: bool) -> tuple[np.ndarray, np.ndarray, dict, list[dict]]:
    """All valid split-site starts, best score first (stable sort).

    Open mode requires an unfolded dry structure and rejects starts that are
    direct- or near-folded. Folded mode requires a direct-folded dry structure
    with at least two base acceptors, keeps the fold, and asserts that only the
    rotated phosphate H moved.
    """
    source = Path(xyz_path)
    z, dry_xyz = read_xyz(source)
    dry_metrics = analyze_structure(z, dry_xyz)
    if dry_metrics["structural_parse"] != "PASS" or dry_metrics["direct_fold"] != folded:
        raise ValueError(f"Expected a {'direct-folded' if folded else 'open'} dry structure: {source}")
    edges, adjacency = bond_graph(z, dry_xyz)
    p_atom = int(np.where(z == 15)[0][0])
    donors = [
        (oxygen, hydrogen)
        for oxygen in adjacency[p_atom] if z[oxygen] == 8
        for hydrogen in adjacency[oxygen] if z[hydrogen] == 1
    ]
    acceptors = [int(value) - 1 for value in dry_metrics["base_acceptor_atoms"].split(";") if value]
    if len(donors) != 2 or len(acceptors) < (2 if folded else 1):
        raise ValueError(f"Unexpected donor/acceptor count in {source}")

    phosphate_sets = {
        (donor_o, donor_h): phosphate_water_candidates(z, dry_xyz, donor_o, donor_h, p_atom)
        for donor_o, donor_h in donors
    }
    base_sets = {
        acceptor: base_water_candidates(z, dry_xyz, acceptor, adjacency)
        for acceptor in acceptors
    }
    candidates: list[dict] = []
    n = len(z)
    full_z = np.concatenate([z, [8, 1, 1, 8, 1, 1]])
    for (donor_o, donor_h), phosphate_options in phosphate_sets.items():
        for acceptor, base_options in base_sets.items():
            for phosphate in phosphate_options:
                for base in base_options:
                    solute_xyz = phosphate["solute_xyz"]
                    water1, water2 = phosphate["water"], base["water"]
                    full_xyz = np.vstack([solute_xyz, water1, water2])
                    p_to_water = hbond_geometry(donor_o, donor_h, n, full_xyz)
                    water_to_base = hbond_geometry(n + 3, n + 4, acceptor, full_xyz)
                    clearance = full_water_clearance(full_z, full_xyz, n, edges, ((donor_h, n), (n + 4, acceptor)))
                    metrics = analyze_structure(full_z, full_xyz)
                    water_oo = float(np.linalg.norm(water1[0] - water2[0]))
                    # Exactly one phosphate->water and one water->base contact pair, no bridging water.
                    if not (
                        metrics["structural_parse"] == "PASS"
                        and (folded or not (metrics["direct_fold"] or metrics["near_fold"]))
                        and metrics["intact_water_count"] == 2
                        and metrics["phosphate_water_contact_count"] == 1
                        and metrics["water_base_contact_count"] == 1
                        and metrics["bridging_water_count"] == 0
                        and is_hbond(p_to_water)
                        and is_hbond(water_to_base)
                        and water_oo >= 2.50
                        and clearance[0] >= 0.70
                    ):
                        continue
                    candidate = {
                        "z": full_z,
                        "xyz": full_xyz,
                        "score": min(phosphate["clearance"], base["clearance"], clearance[0]),
                        "donor_o": donor_o,
                        "donor_h": donor_h,
                        "acceptor": acceptor,
                        "phosphate_chi": phosphate["chi"],
                        "phosphate_water_psi": phosphate["psi"],
                        "base_direction_id": base["direction_id"],
                        "base_water_psi": base["psi"],
                        "p_to_water": p_to_water,
                        "water_to_base": water_to_base,
                        "metrics": metrics,
                        "closest_unintended": clearance,
                        "water_O_O_A": water_oo,
                    }
                    if folded:
                        heavy_displacement = float(
                            np.max(np.linalg.norm(solute_xyz[z != 1] - dry_xyz[z != 1], axis=1))
                        )
                        changed_atoms = [
                            index + 1
                            for index, displacement in enumerate(np.linalg.norm(solute_xyz - dry_xyz, axis=1))
                            if displacement > 1.0e-6
                        ]
                        if heavy_displacement > 1.0e-8 or changed_atoms != [donor_h + 1]:
                            raise RuntimeError("Split-site construction changed atoms beyond the selected phosphate H")
                        candidate["max_solute_heavy_displacement_A"] = heavy_displacement
                        candidate["changed_solute_atoms"] = ";".join(map(str, changed_atoms))
                    candidates.append(candidate)
    if not candidates:
        raise RuntimeError(f"No valid split-site candidates for {family}")
    return z, dry_xyz, dry_metrics, sorted(candidates, key=lambda item: item["score"], reverse=True)


def choose_two(candidates: list[dict], solute_atom_count: int, require_different_acceptor: bool) -> list[dict]:
    """Best start plus the best distinct second start (different phosphate OH, and
    different base acceptor when required), rewarding displacement by 0.03 per angstrom."""
    first = candidates[0]
    alternatives = []
    for candidate in candidates[1:]:
        if candidate["donor_h"] == first["donor_h"]:
            continue
        if require_different_acceptor and candidate["acceptor"] == first["acceptor"]:
            continue
        p_sep = float(
            np.linalg.norm(candidate["xyz"][solute_atom_count] - first["xyz"][solute_atom_count])
        )
        b_sep = float(
            np.linalg.norm(candidate["xyz"][solute_atom_count + 3] - first["xyz"][solute_atom_count + 3])
        )
        if p_sep < 1.25 or b_sep < 1.25:
            continue
        alternatives.append((candidate["score"] + 0.03 * (p_sep + b_sep), candidate))
    if not alternatives:
        raise RuntimeError("Could not select a second, distinct split-site placement")
    return [first, max(alternatives, key=lambda item: item[0])[1]]
