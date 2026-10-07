"""Shared bond-graph, ring, hydrogen-bond and connectivity-fingerprint helpers.

Two distance conventions are kept deliberately: the bond graph uses
0.1 < d < 1.25 * sum(r) and the fingerprint 0.1 < d <= 1.25 * sum(r).
Unifying them would change stored fingerprints.
"""
from __future__ import annotations
import hashlib
import math
import re

import numpy as np

from .constants import (BOND_SCALE, COVALENT_RADII, HB_ANGLE, HB_DA, HB_HA, MIN_BOND_DISTANCE,
                        NEAR_ANGLE, NEAR_DA, NEAR_HA, SYMBOLS)

ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")


def formula_counts(formula: str | None) -> dict[str, int]:
    """Element counts of a Hill-style formula such as C5H5N5 or Gaussian's H2O(1+).

    A trailing parenthesized charge is ignored and repeated elements are summed.
    Returns {} for an empty formula or one with characters outside element/count tokens.
    """
    if not formula:
        return {}
    formula = re.sub(r"\([^)]*\)$", "", formula)
    parts = ELEMENT.findall(formula)
    if not parts or "".join(element + count for element, count in parts) != formula:
        return {}
    counts: dict[str, int] = {}
    for element, count in parts:
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


def formula_key(formula: str | None) -> str:
    """Order-independent formula identity, e.g. 'C:5|H:5|N:5'."""
    return "|".join(f"{element}:{count}" for element, count in sorted(formula_counts(formula).items()))


def bond_graph(z: np.ndarray, xyz: np.ndarray) -> tuple[set[tuple[int, int]], list[set[int]]]:
    distances = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=-1)
    edges = {
        (i, j) for i in range(len(z)) for j in range(i + 1, len(z))
        if z[i] in COVALENT_RADII and z[j] in COVALENT_RADII
        and MIN_BOND_DISTANCE < distances[i, j] < BOND_SCALE * (COVALENT_RADII[z[i]] + COVALENT_RADII[z[j]])
    }
    adjacency = [set() for _ in z]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    return edges, adjacency


def components(adjacency: list[set[int]]) -> list[set[int]]:
    unseen = set(range(len(adjacency)))
    output = []
    while unseen:
        first = unseen.pop()
        component, stack = {first}, [first]
        while stack:
            for neighbor in adjacency[stack.pop()] - component:
                component.add(neighbor)
                unseen.discard(neighbor)
                stack.append(neighbor)
        output.append(component)
    return output


def alternative_path(adjacency: list[set[int]], start: int, target: int, edge: tuple[int, int]) -> bool:
    stack, visited = [start], {start}
    while stack:
        node = stack.pop()
        for neighbor in adjacency[node]:
            if {node, neighbor} == set(edge):
                continue
            if neighbor == target:
                return True
            if neighbor not in visited:
                visited.add(neighbor)
                stack.append(neighbor)
    return False


def ring_systems(z: np.ndarray, edges: set[tuple[int, int]], adjacency: list[set[int]], solute: set[int]) -> list[set[int]]:
    """Connected heavy-atom ring systems within one solute component."""
    cycle_edges = [
        (i, j) for i, j in edges
        if i in solute and j in solute and z[i] != 1 and z[j] != 1
        and alternative_path(adjacency, i, j, (i, j))
    ]
    cycle_adjacency = [set() for _ in z]
    for i, j in cycle_edges:
        cycle_adjacency[i].add(j)
        cycle_adjacency[j].add(i)
    remaining = {node for edge in cycle_edges for node in edge}
    rings = []
    while remaining:
        first = remaining.pop()
        ring, stack = {first}, [first]
        while stack:
            for neighbor in cycle_adjacency[stack.pop()] - ring:
                ring.add(neighbor)
                remaining.discard(neighbor)
                stack.append(neighbor)
        rings.append(ring)
    return rings


def heterocycle_ring(z: np.ndarray, edges: set[tuple[int, int]], adjacency: list[set[int]], solute: set[int]) -> set[int]:
    """Largest nitrogen-containing ring system: the nucleobase."""
    rings = [ring for ring in ring_systems(z, edges, adjacency, solute) if any(z[i] == 7 for i in ring)]
    if not rings:
        raise ValueError("No nitrogen-containing ring found")
    return max(rings, key=len)


def hbond_geometry(donor: int, hydrogen: int, acceptor: int, xyz: np.ndarray) -> dict[str, float]:
    left, right = xyz[donor] - xyz[hydrogen], xyz[acceptor] - xyz[hydrogen]
    cosine = np.clip(left @ right / np.linalg.norm(left) / np.linalg.norm(right), -1.0, 1.0)
    return {
        "H_A_A": float(np.linalg.norm(xyz[hydrogen] - xyz[acceptor])),
        "D_A_A": float(np.linalg.norm(xyz[donor] - xyz[acceptor])),
        "D_H_A_deg": float(np.degrees(np.arccos(cosine))),
    }


def is_hbond(values: dict[str, float], near: bool = False) -> bool:
    h_a, d_a, angle = (NEAR_HA, NEAR_DA, NEAR_ANGLE) if near else (HB_HA, HB_DA, HB_ANGLE)
    return values["H_A_A"] <= h_a and values["D_A_A"] <= d_a and values["D_H_A_deg"] >= angle


def fingerprint(geometry: list) -> str | None:
    """Element-labelled connectivity hash of [[Z, x, y, z], ...] in angstrom."""
    if not geometry:
        return None
    atoms = [(SYMBOLS[int(z)], int(z), tuple(xyz)) for z, *xyz in geometry]
    neighbors: list[list[int]] = [[] for _ in atoms]
    for i, (_, left, point) in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            _, right, other = atoms[j]
            if MIN_BOND_DISTANCE < math.dist(point, other) <= BOND_SCALE * (COVALENT_RADII[left] + COVALENT_RADII[right]):
                neighbors[i].append(j)
                neighbors[j].append(i)
    labels = [symbol for symbol, _, _ in atoms]
    for _ in range(6):
        labels = [hashlib.sha256((labels[i] + "|" + "|".join(sorted(labels[j] for j in linked))).encode()).hexdigest()[:20]
                  for i, linked in enumerate(neighbors)]
    payload = "|".join(sorted(labels)) + ";" + "|".join(str(len(items)) for items in sorted(neighbors, key=len))
    return hashlib.sha256(payload.encode()).hexdigest()
