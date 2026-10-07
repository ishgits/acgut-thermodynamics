"""Portable calculation bundles: explicit units, stable IDs, and checked file IO."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .constants import SYMBOLS

NUMBERS = {v: k for k, v in SYMBOLS.items()}


def digest(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def within(root: Path, relative: str) -> Path:
    """Prevent a manifest from reading outside its repository."""
    candidate = (root / relative).resolve()
    if not candidate.is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes repository: {relative}")
    return candidate


def read_xyz(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    lines = Path(path).read_text().splitlines()
    while lines and not lines[-1].strip():
        lines.pop()
    n = int(lines[0])
    if n < 1 or len(lines) != n + 2:
        raise ValueError(f"Expected one complete XYZ geometry: {path}")
    rows = [line.split() for line in lines[2:]]
    if any(len(row) != 4 or row[0] not in NUMBERS for row in rows):
        raise ValueError(f"Unsupported or malformed XYZ atoms: {path}")
    z = np.array([NUMBERS[row[0]] for row in rows], dtype=int)
    xyz = np.array([[float(v) for v in row[1:]] for row in rows])
    if not np.isfinite(xyz).all():
        raise ValueError(f"Nonfinite XYZ coordinates: {path}")
    return z, xyz


def write_xyz(path: str | Path, geometry: list[list[float]], comment: str) -> None:
    text = [str(len(geometry)), comment]
    text.extend(f"{SYMBOLS[int(z)]} {x:.6f} {y:.6f} {v:.6f}" for z, x, y, v in geometry)
    Path(path).write_text("\n".join(text) + "\n")


def read_bundle(root: Path, relative: str, expected_sha256: str | None = None) -> dict:
    path = within(root, relative)
    if expected_sha256 and digest(path) != expected_sha256:
        raise ValueError(f"Calculation record checksum mismatch: {relative}")
    record = json.loads(path.read_text())
    if record.get("schema_version") != "1.0":
        raise ValueError(f"Unsupported record schema: {relative}")
    if path.parent.name != record["calculation_id"]:
        raise ValueError(f"Calculation ID and directory disagree: {relative}")
    for name, sha in record["files"].items():
        target = within(path.parent, name)
        if digest(target) != sha:
            raise ValueError(f"Calculation artifact checksum mismatch: {target.name}")
    if "frequencies.csv" not in record["files"]:
        raise ValueError("Frequency artifact must have a declared checksum")
    freq_path = within(path.parent, "frequencies.csv")
    with freq_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    modes = [int(row["mode_index"]) for row in rows]
    if modes != list(range(1, len(rows) + 1)):
        raise ValueError(f"Missing, reordered, or duplicated frequency modes: {relative}")
    frequencies = [float(row["frequency_cm1"]) for row in rows]
    if not np.isfinite(frequencies).all():
        raise ValueError(f"Nonfinite frequencies: {relative}")
    record["frequencies_cm1"] = frequencies
    geometry_name = record.get("geometry_file")
    if geometry_name:
        if geometry_name not in record["files"]:
            raise ValueError("Geometry artifact must have a declared checksum")
        z, xyz = read_xyz(within(path.parent, geometry_name))
        record["geometry"] = [[int(a), *point.tolist()] for a, point in zip(z, xyz)]
    else:
        record["geometry"] = []
    return record


def write_bundle(root: Path, record: dict) -> str:
    relative = f"data/calculations/{record['calculation_id']}/record.json"
    path = within(root, relative)
    if path.exists():
        raise FileExistsError(f"Record already exists: {relative}")
    path.parent.mkdir(parents=True, exist_ok=True)
    record = dict(record)
    for field in ("route", "route_optimization", "route_frequency"):
        record.pop(field, None)
    frequencies = record.pop("frequencies_cm1", [])
    geometry = record.pop("geometry", [])
    initial = record.pop("initial_geometry", [])
    record["schema_version"] = "1.0"
    record["files"] = {}
    freq_path = path.parent / "frequencies.csv"
    with freq_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mode_index", "frequency_cm1"])
        writer.writerows((i, value) for i, value in enumerate(frequencies, 1))
    record["files"][freq_path.name] = digest(freq_path)
    if geometry:
        name = "optimized.xyz" if record["status"] != "failed" else "incomplete.xyz"
        write_xyz(path.parent / name, geometry,
                  f"{record['calculation_id']}; angstrom; atom order from source job; {record['status']}")
        record["geometry_file"] = name
        record["files"][name] = digest(path.parent / name)
    else:
        record["geometry_file"] = None
    if initial:
        name = "starting.xyz"
        write_xyz(path.parent / name, initial,
                  f"{record['calculation_id']}; angstrom; first printed source-job orientation; unoptimized start")
        record["files"][name] = digest(path.parent / name)
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    return relative
