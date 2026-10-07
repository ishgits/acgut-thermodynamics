"""Phosphate/water/base contact classifier shared by water-start preparation and analysis.

Contacts are listed per donor-hydrogen/acceptor pair. Preparation filters on pair
counts (``phosphate_water_contact_count``); analysis reports distinct waters
(``analysis.structure``). Both derive from the same contact lists.
"""
import numpy as np

from .constants import SYMBOLS as PT
from .geometry import bond_graph, components, heterocycle_ring, hbond_geometry as hb_geometry, is_hbond


def analyze_structure(z, xyz):
    """Classify one geometry; returns ``structural_parse`` PASS or FAIL with a ``reason``."""
    edges, adjacency = bond_graph(z, xyz)
    comps = components(adjacency)
    water_components = [
        comp for comp in comps if sorted(z[list(comp)].tolist()) == [1, 1, 8]
    ]
    p_atoms = [i for i, atomic_number in enumerate(z) if atomic_number == 15]
    if len(p_atoms) != 1:
        return {"structural_parse": "FAIL", "P_count": len(p_atoms), "reason": f"Expected one phosphorus atom; found {len(p_atoms)}"}
    p_atom = p_atoms[0]
    solute = next(comp for comp in comps if p_atom in comp)
    try:
        base_ring = heterocycle_ring(z, edges, adjacency, solute)
    except ValueError:
        return {"structural_parse": "FAIL", "P_count": 1, "reason": "no N-containing ring"}
    base_region = set(base_ring)
    for node in base_ring:
        for neighbor in adjacency[node]:
            if neighbor in solute and z[neighbor] in (7, 8):
                base_region.add(neighbor)

    # Chemically filtered geometric acceptor set. Exocyclic NHx and glycosidic
    # ring N atoms with three heavy neighbors are excluded.
    acceptors = []
    for atom in sorted(base_region):
        hydrogens = [neighbor for neighbor in adjacency[atom] if z[neighbor] == 1]
        heavy = [neighbor for neighbor in adjacency[atom] if z[neighbor] != 1]
        if z[atom] == 8 and not hydrogens:
            acceptors.append(atom)
        elif z[atom] == 7 and not hydrogens and len(heavy) <= 2:
            acceptors.append(atom)

    phosphate_donors = []
    for oxygen in adjacency[p_atom]:
        if z[oxygen] != 8:
            continue
        for hydrogen in adjacency[oxygen]:
            if z[hydrogen] == 1:
                phosphate_donors.append((oxygen, hydrogen))

    waters = []
    for water_id, comp in enumerate(water_components, start=1):
        oxygen = next(i for i in comp if z[i] == 8)
        hydrogens = sorted(i for i in comp if z[i] == 1)
        waters.append({"water_id": water_id, "O": oxygen, "H": hydrogens})

    direct_contacts = []
    near_contacts = []
    closest_base = []
    phosphate_to_water = []
    water_to_base = []
    for oxygen, hydrogen in phosphate_donors:
        for acceptor in acceptors:
            values = hb_geometry(oxygen, hydrogen, acceptor, xyz)
            row = {
                "donor_O_atom": oxygen + 1,
                "donor_H_atom": hydrogen + 1,
                "acceptor_atom": acceptor + 1,
                "acceptor_element": PT[z[acceptor]],
                **values,
            }
            if is_hbond(values):
                direct_contacts.append(row)
            if is_hbond(values, near=True):
                near_contacts.append(row)
        if acceptors:
            acceptor = min(acceptors, key=lambda a: np.linalg.norm(xyz[hydrogen] - xyz[a]))
            closest_base.append({
                "donor_O_atom": oxygen + 1,
                "donor_H_atom": hydrogen + 1,
                "acceptor_atom": acceptor + 1,
                "acceptor_element": PT[z[acceptor]],
                **hb_geometry(oxygen, hydrogen, acceptor, xyz),
            })
        for water in waters:
            values = hb_geometry(oxygen, hydrogen, water["O"], xyz)
            if is_hbond(values):
                phosphate_to_water.append({
                    "water_id": water["water_id"],
                    "donor_O_atom": oxygen + 1,
                    "donor_H_atom": hydrogen + 1,
                    "water_O_atom": water["O"] + 1,
                    **values,
                })
    for water in waters:
        for hydrogen in water["H"]:
            for acceptor in acceptors:
                values = hb_geometry(water["O"], hydrogen, acceptor, xyz)
                if is_hbond(values):
                    water_to_base.append({
                        "water_id": water["water_id"],
                        "water_O_atom": water["O"] + 1,
                        "water_H_atom": hydrogen + 1,
                        "acceptor_atom": acceptor + 1,
                        "acceptor_element": PT[z[acceptor]],
                        **values,
                    })
    p_to_water_ids = {row["water_id"] for row in phosphate_to_water}
    water_to_base_ids = {row["water_id"] for row in water_to_base}
    bridging_ids = sorted(p_to_water_ids & water_to_base_ids)
    direct = bool(direct_contacts)
    if direct:
        state = "direct_phosphate_OH_to_base_fold"
    elif bridging_ids:
        state = "open_water_bridge_between_phosphate_and_base"
    elif p_to_water_ids and water_to_base_ids:
        state = "open_water_competes_at_both_sites"
    elif p_to_water_ids:
        state = "open_water_accepts_phosphate_OH_only"
    elif water_to_base_ids:
        state = "open_water_hydrates_base_only"
    else:
        state = "open_no_targeted_water_bridge"

    closest = min(closest_base, key=lambda row: row["H_A_A"]) if closest_base else {}
    ring_centroid = xyz[list(base_ring)].mean(axis=0)
    return {
        "structural_parse": "PASS",
        "natoms": len(z),
        "fragment_sizes": ";".join(str(len(comp)) for comp in sorted(comps, key=len)),
        "intact_water_count": len(waters),
        "phosphate_OH_count": len(phosphate_donors),
        "base_ring_atoms": ";".join(str(i + 1) for i in sorted(base_ring)),
        "base_acceptor_atoms": ";".join(str(i + 1) for i in acceptors),
        "direct_fold": direct,
        "near_fold": bool(near_contacts),
        "fold_state": state,
        "direct_fold_contact_count": len(direct_contacts),
        "phosphate_water_contact_count": len(phosphate_to_water),
        "water_base_contact_count": len(water_to_base),
        "bridging_water_count": len(bridging_ids),
        "bridging_water_ids": ";".join(map(str, bridging_ids)),
        "closest_phosphate_H_to_base_acceptor_A": closest.get("H_A_A", np.nan),
        "closest_contact_angle_deg": closest.get("D_H_A_deg", np.nan),
        "closest_contact_acceptor_atom": closest.get("acceptor_atom", np.nan),
        "closest_contact_acceptor_element": closest.get("acceptor_element", ""),
        "P_to_base_ring_centroid_A": float(np.linalg.norm(xyz[p_atom] - ring_centroid)),
        "P_to_nearest_base_ring_atom_A": float(min(np.linalg.norm(xyz[p_atom] - xyz[i]) for i in base_ring)),
        "solute_heavy_radius_gyration_A": float(
            np.sqrt(np.mean(np.sum((xyz[[i for i in solute if z[i] != 1]] - xyz[[i for i in solute if z[i] != 1]].mean(axis=0)) ** 2, axis=1)))
        ),
        "direct_contacts": direct_contacts,
        "closest_base_contacts": closest_base,
        "phosphate_to_water_contacts": phosphate_to_water,
        "water_to_base_contacts": water_to_base,
    }
