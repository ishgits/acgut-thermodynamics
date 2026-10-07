"""Content-based Gaussian parsing for optimized frequency calculations."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

from ..constants import COVALENT_RADII, SYMBOLS, TEMPERATURE_K, TEMPERATURE_TOLERANCE_K
from ..geometry import fingerprint, formula_counts, formula_key


FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][+-]?\d+)?"
ROUTE_START = re.compile(r"^\s*#\s*.+$", re.MULTILINE)
SEPARATOR = re.compile(r"^\s*-{10,}\s*$")
METHOD = re.compile(r"([A-Za-z0-9ωΩ-]+)\s*/\s*([A-Za-z0-9+*(),-]+)", re.IGNORECASE)
SOLVENT = re.compile(r"solvent\s*=\s*([A-Za-z0-9_-]+)", re.IGNORECASE)
STOICHIOMETRY = re.compile(r"Stoichiometry\s+([A-Za-z0-9()+-]+)")
CHARGE_MULT = re.compile(r"Charge\s*=\s*([+-]?\d+)\s+Multiplicity\s*=\s*(\d+)")
FREQUENCY = re.compile(r"Frequencies\s+--\s+(.+)")
TEMPERATURE_PRESSURE = re.compile(
    r"Temperature\s+(" + FLOAT + r")\s+Kelvin\.\s+Pressure\s+(" + FLOAT + r")\s+Atm\."
)
SCF = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+(" + FLOAT + r")")
THERMAL_G = re.compile(r"Thermal correction to Gibbs Free Energy=\s+(" + FLOAT + r")")
PRINTED_G = re.compile(r"Sum of electronic and thermal Free Energies=\s+(" + FLOAT + r")")
STANDARD_ORIENTATION = re.compile(r"^\s*Standard orientation:\s*$", re.MULTILINE)
ORIENTATION = re.compile(
    r"(?:Standard|Input) orientation:\s*\n\s*-+\n[^\n]*\n[^\n]*\n\s*-+\n(.*?)\n\s*-+",
    re.S,
)
OPT_COMPLETE = re.compile(
    r"Optimization completed(?:\.|\s+on the basis of negligible forces\.)",
    re.IGNORECASE,
)
NORMAL_TERMINATION = re.compile(r"Normal termination of Gaussian")
ERROR_TERMINATION = re.compile(r"Error termination", re.IGNORECASE)

ATOMIC_SYMBOLS = SYMBOLS
SUPPORTED_ELEMENTS = frozenset(SYMBOLS[z] for z in COVALENT_RADII)


def _number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _last(pattern: re.Pattern[str], text: str) -> float | None:
    matches = pattern.findall(text)
    return _number(matches[-1]) if matches else None


def _route_text(text: str, start: int) -> str:
    lines = text[start:].splitlines()
    result: list[str] = []
    found_separator = False
    for line in lines:
        if result and SEPARATOR.match(line):
            found_separator = True
            break
        # Gaussian prefixes echoed lines with one space. Other whitespace is
        # route content: a wrap can occur between keywords or inside a keyword.
        result.append(line[1:] if line.startswith(" ") else line)
        if len(result) > 20:
            break
    return "".join(result).strip() if found_separator else ""


def _route_parts(route: str) -> tuple[str | None, str | None, str | None, str | None]:
    method = METHOD.search(route)
    compact = re.sub(r"\s+", "", route).casefold()
    solvent = SOLVENT.search(route)
    functional = method.group(1).replace("ω", "w").replace("Ω", "w").upper() if method else None
    basis = method.group(2).upper() if method else None
    model = "SMD" if "smd" in compact else ("PCM" if "pcm" in compact or "iefpcm" in compact else None)
    solvent_name = solvent.group(1).upper() if solvent else None
    return functional, basis, model, solvent_name


def _connectivity_fingerprint(text: str) -> str | None:
    blocks = orientations(text)
    if not blocks or not blocks[-1][1]:
        return None
    return fingerprint(blocks[-1][1])


def orientations(text: str) -> list[tuple[int, list[list[float]]]]:
    """Read Cartesian blocks and their character offsets without inventing coordinates."""
    output = []
    for match in ORIENTATION.finditer(text):
        rows = []
        for line in match.group(1).splitlines():
            fields = line.split()
            if len(fields) != 6:
                raise ValueError("Malformed Gaussian orientation row")
            number = int(fields[1])
            if number not in ATOMIC_SYMBOLS:
                raise ValueError(f"Unsupported atomic number: {number}")
            rows.append([number, *[_number(v) for v in fields[3:6]]])
        output.append((match.start(), rows))
    return output


def parse_log(path: str | Path) -> dict[str, object]:
    """Parse the final frequency job and retain file-wide provenance checks."""
    path = Path(path).resolve()
    raw = path.read_bytes()
    text = raw.decode(errors="ignore")
    routes = list(ROUTE_START.finditer(text))
    route_records = [(_route_text(text, match.start()), match.start()) for match in routes]
    route_records = [(route, start) for route, start in route_records if route]
    frequency_routes = [(route, start) for route, start in route_records if "freq" in route.casefold()]
    thermo_route, thermo_start = frequency_routes[-1] if frequency_routes else (route_records[-1] if route_records else ("", 0))
    later_starts = [start for _, start in route_records if start > thermo_start]
    thermo_end = min(later_starts) if later_starts else len(text)
    thermo_segment = text[thermo_start:thermo_end]
    functional, basis, solvent_model, solvent = _route_parts(thermo_route)

    formula_matches = list(STOICHIOMETRY.finditer(thermo_segment)) or list(STOICHIOMETRY.finditer(text))
    formula = formula_matches[-1].group(1) if formula_matches else None
    counts = formula_counts(formula)
    charge_matches = list(CHARGE_MULT.finditer(thermo_segment)) or list(CHARGE_MULT.finditer(text))
    charge = int(charge_matches[-1].group(1)) if charge_matches else None
    multiplicity = int(charge_matches[-1].group(2)) if charge_matches else None
    frequencies = [_number(token) for match in FREQUENCY.findall(thermo_segment) for token in match.split()]
    tp_matches = TEMPERATURE_PRESSURE.findall(thermo_segment)
    temperature = _number(tp_matches[-1][0]) if tp_matches else None
    pressure = _number(tp_matches[-1][1]) if tp_matches else None
    electronic = _last(SCF, thermo_segment)
    correction = _last(THERMAL_G, thermo_segment)
    printed = _last(PRINTED_G, thermo_segment)
    computed = electronic + correction if electronic is not None and correction is not None else None
    residual = computed - printed if computed is not None and printed is not None else None
    normal_after_thermo = bool(NORMAL_TERMINATION.search(thermo_segment))
    geometry_blocks = orientations(thermo_segment)
    geometry_source = "frequency_job" if frequency_routes else "unfinished_job_no_frequency"
    geometry_offset = thermo_start
    if not geometry_blocks:
        previous_start = max([start for _, start in route_records if start < thermo_start], default=0)
        geometry_blocks = orientations(text[previous_start:thermo_start])
        geometry_source = "preceding_job_before_frequency"
        geometry_offset = previous_start
    geometry = geometry_blocks[-1][1] if geometry_blocks else []
    geometry_line = (
        text[:geometry_offset + geometry_blocks[-1][0]].count("\n") + 1
        if geometry_blocks else None
    )
    initial_blocks = orientations(text[:thermo_start] if frequency_routes and thermo_start else text)
    initial_geometry = initial_blocks[0][1] if initial_blocks else []

    failures: list[str] = []
    warnings: list[str] = []
    if not frequency_routes:
        failures.append("missing_frequency_job")
    if not normal_after_thermo:
        failures.append("frequency_job_not_normally_terminated")
    if ERROR_TERMINATION.search(text[thermo_start:]):
        failures.append("error_termination_after_frequency_job")
    if not OPT_COMPLETE.search(text[:thermo_start]):
        failures.append("missing_completed_optimization_before_frequency_job")
    if None in (functional, basis, solvent_model, solvent):
        failures.append("incomplete_method_fingerprint")
    if not counts or not set(counts).issubset(SUPPORTED_ELEMENTS):
        failures.append("missing_or_unsupported_stoichiometry")
    if charge is None or multiplicity is None:
        failures.append("missing_charge_or_multiplicity")
    if not frequencies:
        failures.append("missing_frequencies")
    if any(value < 0 for value in frequencies):
        failures.append("imaginary_frequency_present")
    # Mode counts differ for linear molecules and monatomic species.
    linear = False
    if len(geometry) > 1:
        import numpy as np
        xyz = np.array([row[1:] for row in geometry])
        singular = np.linalg.svd(xyz - xyz.mean(axis=0), compute_uv=False)
        linear = bool(len(geometry) == 2 or singular[1] < 1e-6)
    natoms = sum(counts.values()) if counts else None
    expected_modes = (0 if natoms == 1 else 3 * natoms - (5 if linear else 6)) if natoms else None
    if expected_modes is not None and len(frequencies) != expected_modes:
        failures.append("frequency_count_not_expected_for_geometry")
    if temperature is None or abs(temperature - TEMPERATURE_K) > TEMPERATURE_TOLERANCE_K:
        failures.append("temperature_is_not_298K")
    if pressure is None:
        warnings.append("pressure_not_found")
    if None in (electronic, correction, printed, computed):
        failures.append("missing_gibbs_thermochemistry")
    elif residual is not None and abs(residual) > 2e-5:
        failures.append("computed_and_printed_gibbs_disagree")
    connectivity = _connectivity_fingerprint(thermo_segment)
    if connectivity is None and geometry_source == "preceding_job_before_frequency":
        connectivity = _connectivity_fingerprint(text[previous_start:thermo_start])
    if connectivity is None:
        failures.append("missing_connectivity_fingerprint")

    return {
        "log_path": str(path),
        "file_name": path.name,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "status": "failed" if failures else ("warning" if warnings else "valid"),
        "failures": "; ".join(failures),
        "warnings": "; ".join(warnings),
        "functional": functional,
        "basis_set": basis,
        "solvent_model": solvent_model,
        "solvent": solvent,
        "route_fingerprint": "|".join(str(v or "UNKNOWN") for v in (functional, basis, solvent_model, solvent)),
        "formula": formula,
        "formula_key": formula_key(formula),
        "connectivity_fingerprint": connectivity,
        "charge": charge,
        "multiplicity": multiplicity,
        "atom_count": sum(counts.values()) if counts else None,
        "normal_termination_count": len(NORMAL_TERMINATION.findall(text)),
        "optimization_completed": bool(OPT_COMPLETE.search(text[:thermo_start])),
        "frequency_count": len(frequencies),
        "expected_frequency_count": expected_modes,
        "imaginary_frequency_count": sum(value < 0 for value in frequencies),
        "minimum_frequency_cm1": min(frequencies) if frequencies else None,
        "temperature_k": temperature,
        "pressure_atm": pressure,
        "electronic_energy_hartree": electronic,
        "thermal_gibbs_correction_hartree": correction,
        "computed_gibbs_hartree": computed,
        "printed_gibbs_hartree": printed,
        "gibbs_residual_hartree": residual,
        "frequencies_cm1": frequencies,
        "geometry": geometry,
        "initial_geometry": initial_geometry,
        "geometry_source": geometry_source,
        "geometry_source_line": geometry_line,
        "geometry_precision_angstrom": 1e-6,
        "atom_indices": "1-based source orientation order",
        "linear_geometry": linear,
        "thermo_segment_first_line": text[:thermo_start].count("\n") + 1,
        "thermo_segment_last_line": text[:thermo_end].count("\n") + 1,
        "gaussian_version": next(iter(re.findall(r"Gaussian\s+16[^\n]*Revision[^\n]*", text)), None),
    }
