"""Explicit neutral condensation with audited atom choices and reproducible seeds.

RDKit is optional; these are starting structures, not quantum-optimized results.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

from ..paths import fresh


def load_sdf(path: Path):
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(path), removeHs=False)
    if mol is None:
        raise ValueError(f"Cannot read SDF: {path}")
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        raise ValueError("Construction indices refer to the bundled explicit-H SDF; do not silently add/reorder atoms")
    return mol


def condense(base, donor, *, base_atom: int, donor_carbon: int, donor_oxygen: int,
             seed: int = 61453, expected_formula: str | None = None):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    base, donor = Chem.Mol(base), Chem.Mol(donor)
    if Chem.GetFormalCharge(base) != 0 or Chem.GetFormalCharge(donor) != 0:
        raise ValueError("This study builder supports neutral condensation only; define charged chemistry explicitly")
    if donor.GetAtomWithIdx(donor_carbon).GetSymbol() != "C":
        raise ValueError("Donor attachment must be the explicitly selected carbon")
    if donor.GetAtomWithIdx(donor_oxygen).GetSymbol() != "O" or not donor.GetBondBetweenAtoms(donor_carbon, donor_oxygen):
        raise ValueError("Selected donor OH is not attached to the selected carbon")
    def hydrogens(mol, index):
        return sorted(a.GetIdx() for a in mol.GetAtomWithIdx(index).GetNeighbors() if a.GetAtomicNum() == 1)
    base_h, oh_h = hydrogens(base, base_atom), hydrogens(donor, donor_oxygen)
    if not base_h or len(oh_h) != 1:
        raise ValueError("Condensation requires removable base-H and donor-OH")
    nbase = base.GetNumAtoms()
    base.GetAtomWithIdx(base_atom).SetAtomMapNum(9001)
    donor.GetAtomWithIdx(donor_carbon).SetAtomMapNum(9002)
    combined = Chem.RWMol(Chem.CombineMols(base, donor))
    combined.AddBond(base_atom, nbase + donor_carbon, Chem.BondType.SINGLE)
    for i in sorted([base_h[0], nbase + donor_oxygen, nbase + oh_h[0]], reverse=True):
        combined.RemoveAtom(i)
    product = combined.GetMol()
    Chem.SanitizeMol(product)
    formula = rdMolDescriptors.CalcMolFormula(product)
    if expected_formula and formula != expected_formula:
        raise ValueError(f"Product formula {formula} differs from configured {expected_formula}")
    if Chem.GetFormalCharge(product) != 0 or len(Chem.GetMolFrags(product)) != 1:
        raise ValueError("Expected one neutral product")
    product.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    status = AllChem.EmbedMolecule(product, params)
    if status != 0:
        status = AllChem.EmbedMolecule(product, randomSeed=int(seed), useRandomCoords=True)
    if status != 0:
        raise ValueError("RDKit embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(product):
        ff = "MMFF94"
        converged = AllChem.MMFFOptimizeMolecule(product, maxIters=2000) == 0
    else:
        ff = "UFF"
        converged = AllChem.UFFOptimizeMolecule(product, maxIters=2000) == 0
    metadata = {"formula": formula, "charge": 0, "seed": int(seed), "force_field": ff,
                "force_field_converged": bool(converged), "atom_indices": "0-based explicit-H source SDF",
                "base_attachment_atom": base_atom, "donor_carbon": donor_carbon, "donor_oxygen": donor_oxygen,
                "removed_base_hydrogen": base_h[0], "removed_donor_hydrogen": oh_h[0],
                "validation_scope": "neutral formula/linkage construction; review absolute stereochemistry before calculations"}
    return product, metadata


def build_from_registry(root: Path, registry: Path, output: Path) -> Path:
    from rdkit import Chem, rdBase
    from ..records import digest, within
    fresh(output)
    rows = list(csv.DictReader(registry.open()))
    summary = []
    for row in rows:
        base_path, donor_path = within(root, row["base_sdf"]), within(root, row["donor_sdf"])
        product, metadata = condense(load_sdf(base_path), load_sdf(donor_path),
            base_atom=int(row["base_attachment_atom"]), donor_carbon=int(row["donor_carbon"]),
            donor_oxygen=int(row["donor_oxygen"]), seed=int(row["seed"]), expected_formula=row["expected_formula"])
        name = row["species_id"]
        if not name or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in name):
            raise ValueError("Use a safe, unique species ID")
        if (output / f"{name}.sdf").exists():
            raise ValueError("Construction registry contains duplicate species IDs")
        for atom in product.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.MolToMolFile(product, str(output / f"{name}.sdf"))
        (output / f"{name}.xyz").write_text(Chem.MolToXYZBlock(product))
        metadata.update(species_id=name, rdkit_version=rdBase.rdkitVersion,
                        base_sdf_sha256=digest(base_path), donor_sdf_sha256=digest(donor_path),
                        sdf_sha256=digest(output / f"{name}.sdf"), xyz_sha256=digest(output / f"{name}.xyz"))
        (output / f"{name}.json").write_text(json.dumps(metadata, indent=2) + "\n")
        summary.append(metadata)
    (output / "construction_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    return output
