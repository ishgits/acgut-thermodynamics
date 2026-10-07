"""Gaussian opt→frequency inputs, adapted from the cited PubChem pipeline."""
from __future__ import annotations
import json
import re
from pathlib import Path

from ..records import SYMBOLS, read_xyz, read_bundle, digest
from ..paths import fresh
import numpy as np


def method_routes(method: dict) -> tuple[str, str]:
    """Build default-Opt rerun routes from a configured scientific method."""
    functional, basis = method["functional"].lower(), method["basis_set"].lower()
    model = {"PCM": "iefpcm", "SMD": "smd"}[method["solvent_model"].upper()]
    level = f"{functional}/{basis} scrf=({model},solvent={method['solvent'].lower()})"
    return (
        f"# opt {level}",
        f"# freq {level} temperature={method['temperature_k']:g} geom=allchk guess=read",
    )


def input_text(name: str, z, xyz, *, route_opt: str, route_freq: str,
               charge: int = 0, multiplicity: int = 1, cpus: int = 8, memory: str = "6GB") -> str:
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or not re.fullmatch(r"\d+(?:MB|GB)", memory.upper()):
        raise ValueError("Invalid input basename or Gaussian memory")
    if (cpus < 1 or multiplicity < 1
            or not re.search(r"(?:^|\s)opt(?=\s|=|$)", route_opt, re.I)
            or not re.search(r"(?:^|\s)freq(?=\s|=|$)", route_freq, re.I)):
        raise ValueError("Invalid resources or opt/frequency routes")
    if not re.search(r"(?:^|\s)geom\s*=\s*allchk(?=\s|$)", route_freq, re.I):
        raise ValueError("Frequency job must read the optimized geometry with Geom=AllChk")
    if len(z) != len(xyz):
        raise ValueError("Coordinate and element counts differ")
    if not np.isfinite(xyz).all() or not len(z) or any(int(a) not in SYMBOLS for a in z):
        raise ValueError("Input requires finite CHNOP coordinates")
    coordinates = "\n".join(f"{SYMBOLS[int(a)]} {x:.6f} {y:.6f} {v:.6f}" for a, (x, y, v) in zip(z, xyz))
    header = f"%nprocshared={cpus}\n%mem={memory}\n%chk={name}.chk\n"
    return f"{header}{route_opt}\n\n{name}; regenerated study input\n\n{charge} {multiplicity}\n{coordinates}\n\n--Link1--\n{header}{route_freq}\n\n"


def prepare_study_inputs(root: Path, output: Path, cpus: int = 8, memory: str = "6GB") -> Path:
    import csv
    import yaml
    fresh(output)
    rows = list(csv.DictReader((root / "data/calculations/manifest.csv").open()))
    methods = yaml.safe_load((root / "config/methods.yaml").read_text())["methods"]
    method_fields = ("functional", "basis_set", "solvent_model", "solvent")
    prepared = []
    for row in rows:
        record = read_bundle(root, row["record_path"], row["record_sha256"])
        if record["status"] == "failed" or "optimized.xyz" not in record["files"]:
            prepared.append({"calculation_id": row["calculation_id"], "status": "not_prepared", "reason": "failed source or no optimized geometry"})
            continue
        matches = [key for key, method in methods.items()
                   if all(str(record[field]).casefold() == str(method[field]).casefold()
                          for field in method_fields)]
        if len(matches) != 1:
            raise ValueError(f"Expected one configured method for {row['calculation_id']}: {matches}")
        method_id = matches[0]
        method = methods[method_id]
        if record["temperature_k"] != method["temperature_k"] or record["pressure_atm"] != method["pressure_atm"]:
            raise ValueError(f"Thermochemistry state disagrees with configured method: {row['calculation_id']}")
        geometry = (root / row["record_path"]).parent / "optimized.xyz"
        z, xyz = read_xyz(geometry)
        name = row["calculation_id"]
        route_opt, route_freq = method_routes(method)
        text = input_text(name, z, xyz, route_opt=route_opt, route_freq=route_freq,
                          charge=record["charge"], multiplicity=record["multiplicity"], cpus=cpus, memory=memory)
        path = output / f"{name}.com"
        path.write_text(text)
        prepared.append({"calculation_id": name, "status": "prepared", "input_sha256": digest(path),
                         "coordinate_source": str(geometry.relative_to(root)), "geometry_sha256": digest(geometry),
                         "method_id": method_id, "optimization": "Gaussian default Opt",
                         "temperature_k": method["temperature_k"], "cpus": cpus, "memory": memory})
    (output / "input_manifest.json").write_text(json.dumps(prepared, indent=2) + "\n")
    return output


def prepare_input(root: Path, geometry: Path, method_id: str, name: str, output: Path,
                  charge: int = 0, multiplicity: int = 1, cpus: int = 8, memory: str = "6GB") -> Path:
    import yaml
    fresh(output)
    method = yaml.safe_load((root / "config/methods.yaml").read_text())["methods"][method_id]
    opt, freq = method_routes(method)
    z, xyz = read_xyz(geometry)
    text = input_text(name, z, xyz, route_opt=opt, route_freq=freq,
                      charge=charge, multiplicity=multiplicity, cpus=cpus, memory=memory)
    (output / f"{name}.com").write_text(text)
    (output / "input_manifest.json").write_text(json.dumps({"method_id": method_id, "method": method,
        "charge": charge, "multiplicity": multiplicity, "geometry_sha256": digest(geometry),
        "input_sha256": digest(output / f"{name}.com"), "cpus": cpus, "memory": memory}, indent=2) + "\n")
    return output
