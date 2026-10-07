"""Independent line-based extraction used for preparation checks.

Public logs are not bundled; this utility can audit a supplied new log. It
deliberately shares no parsing code with ``atcgu.extraction.gaussian``.
"""
import hashlib
import re

NUMBER = r"[-+]?\d*\.?\d+(?:[DdEe][-+]?\d+)?"
SEPARATOR = r"\s*-{10,}\s*"
THERMO_LABELS = {
    "Zero-point correction=": "zpe_correction_hartree",
    "Thermal correction to Energy=": "thermal_energy_correction_hartree",
    "Thermal correction to Enthalpy=": "thermal_enthalpy_correction_hartree",
    "Thermal correction to Gibbs Free Energy=": "thermal_gibbs_correction_hartree",
    "Sum of electronic and zero-point Energies=": "electronic_plus_zpe_hartree",
    "Sum of electronic and thermal Free Energies=": "printed_gibbs_hartree",
    "Sum of electronic and thermal Enthalpies=": "H_hartree",
}


def num(text):
    return float(text.replace("D", "E").replace("d", "e"))


def parse(path):
    raw = path.read_bytes()
    lines = raw.decode(errors="replace").splitlines()
    result = {"sha256": hashlib.sha256(raw).hexdigest()}
    source_lines = {}
    frequencies = []
    routes = []
    route_starts = []
    orientation_lines = []
    for index, line in enumerate(lines):
        if line.strip().startswith("#"):
            route_parts = []
            for route_line in lines[index:index + 22]:
                if re.fullmatch(SEPARATOR, route_line):
                    break
                route_parts.append(route_line.strip())
            else:
                # No separator within 22 lines: not a complete route echo.
                continue
            routes.append("".join(route_parts))
            route_starts.append(index)
        if "orientation:" in line and ("Standard" in line or "Input" in line):
            orientation_lines.append(index)

    # Last thermochemistry block is authoritative; bound it to its enclosing route.
    thermo_lines = [index for index, line in enumerate(lines)
                    if line.strip().startswith("Thermal correction to Gibbs Free Energy=")]
    pivot = thermo_lines[-1] if thermo_lines else len(lines) - 1
    begin = max([start for start in route_starts if start <= pivot], default=0)
    end = min([start for start in route_starts if start > pivot], default=len(lines))
    segment = lines[begin:end]
    for line_number, line in enumerate(segment, begin + 1):
        stripped = line.strip()
        for label, key in THERMO_LABELS.items():
            if stripped.startswith(label):
                values = re.findall(NUMBER, stripped[len(label):])
                if values:
                    result[key] = num(values[0])
                    source_lines[key] = line_number
        if "SCF Done:" in line:
            value = line.split("=", 1)[1].split()[0]
            result["electronic_energy_hartree"] = num(value)
            source_lines["electronic_energy_hartree"] = line_number
        if stripped.startswith("Frequencies --") or re.match(r"Frequencies\s+--", stripped):
            frequencies.extend(num(value) for value in stripped.split("--")[1].split())
        if stripped.startswith("Stoichiometry "):
            result["formula"] = stripped.split()[1]
            source_lines["formula"] = line_number
        match = re.search(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", line)
        if match:
            result.update(charge=int(match[1]), multiplicity=int(match[2]))
        match = re.search(r"Temperature\s+(" + NUMBER + r")\s+Kelvin\.\s+Pressure\s+(" + NUMBER + r")", line)
        if match:
            result.update(temperature_k=num(match[1]), pressure_atm=num(match[2]))
        if re.fullmatch(r"Total\s+" + NUMBER + r"\s+" + NUMBER + r"\s+" + NUMBER, stripped):
            result["entropy_cal_mol_k"] = num(stripped.split()[-1])
            source_lines["entropy_cal_mol_k"] = line_number

    result.update(
        frequency_count=len(frequencies),
        imaginary_frequency_count=sum(value < 0 for value in frequencies),
        normal_termination_count=sum("Normal termination of Gaussian" in line for line in lines),
        error_termination_count=sum("Error termination" in line for line in lines),
        error_via_count=sum("Error termination via" in line for line in lines),
    )
    if frequencies:
        result["minimum_frequency_cm1"] = min(frequencies)
    if "electronic_energy_hartree" in result and "thermal_gibbs_correction_hartree" in result:
        result["computed_gibbs_hartree"] = result["electronic_energy_hartree"] + result["thermal_gibbs_correction_hartree"]
    if "formula" in result:
        formula = re.sub(r"\([^)]*\)$", "", result["formula"])
        result["atom_count"] = sum(int(count or 1) for _, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula))
    result["thermo_segment_first_line"] = begin + 1
    result["thermo_segment_last_line"] = end
    result["normal_in_thermo_segment"] = any("Normal termination of Gaussian" in line for line in segment)
    result["optimization_completed"] = any("Optimization completed" in line for line in lines[:begin])

    relevant = [index for index in orientation_lines if begin <= index < end]
    if not relevant:
        previous = max([start for start in route_starts if start < begin], default=0)
        relevant = [index for index in orientation_lines if previous <= index < begin]
    if relevant:
        atoms = []
        for line in lines[relevant[-1] + 5:]:
            if re.fullmatch(SEPARATOR, line):
                break
            fields = line.split()
            if len(fields) == 6:
                try:
                    atoms.append([int(fields[1])] + list(map(float, fields[3:])))
                except ValueError:
                    pass
        result["geometry"] = atoms
    result["frequencies"] = frequencies
    result["source_lines"] = source_lines
    return result
