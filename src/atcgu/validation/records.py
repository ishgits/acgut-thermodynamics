"""Scientific checks derived from the released geometry, frequencies and metadata."""
from __future__ import annotations
from collections import Counter
import math

import numpy as np

from ..constants import SYMBOLS as ATOMIC_SYMBOLS, TEMPERATURE_TOLERANCE_K
from ..geometry import fingerprint, formula_counts


def check_record(record: dict, expected: dict | None = None) -> list[str]:
    """Return all violations. Missing facts are failures for usable records."""
    failures = []
    if record.get("status") == "failed":
        return ["source_calculation_failed"]
    frequencies = record["frequencies_cm1"]
    geometry = record["geometry"]
    n = len(geometry)
    if not n:
        return ["missing_geometry"]
    xyz = np.array([row[1:] for row in geometry])
    singular = np.linalg.svd(xyz - xyz.mean(axis=0), compute_uv=False)
    linear = n == 2 or (n > 2 and singular[1] < 1e-6)
    expected_modes = 0 if n == 1 else 3 * n - (5 if linear else 6)
    if len(frequencies) != expected_modes or record.get("frequency_count") != len(frequencies):
        failures.append("mode_count_mismatch")
    if any(f <= 0 for f in frequencies):
        failures.append("nonpositive_frequency")
    if record.get("imaginary_frequency_count") != sum(f < 0 for f in frequencies):
        failures.append("imaginary_count_mismatch")
    if record.get("minimum_frequency_cm1") != (min(frequencies) if frequencies else None):
        failures.append("minimum_frequency_mismatch")
    counts = Counter(ATOMIC_SYMBOLS[int(row[0])] for row in geometry)
    if counts != formula_counts(record.get("formula")) or n != record.get("atom_count"):
        failures.append("geometry_formula_mismatch")
    if fingerprint(geometry) != record.get("connectivity_fingerprint"):
        failures.append("connectivity_fingerprint_mismatch")
    fields = ["electronic_energy_hartree", "thermal_gibbs_correction_hartree", "printed_gibbs_hartree"]
    values = [record.get(key) for key in fields]
    if any(v is None or not math.isfinite(v) for v in values):
        failures.append("missing_or_nonfinite_energy")
    else:
        computed = values[0] + values[1]
        reported = record.get("computed_gibbs_hartree")
        if reported is None or not math.isfinite(reported) or abs(computed - reported) > 1e-12:
            failures.append("gibbs_sum_mismatch")
        if abs(computed - values[2]) > 2e-5:
            failures.append("printed_gibbs_mismatch")
    if not record.get("optimization_completed") or not record.get("frequency_job_normally_terminated"):
        failures.append("unverified_optimization_or_frequency_termination")
    if expected:
        for key in ["functional", "basis_set", "solvent_model", "solvent", "charge", "multiplicity"]:
            if str(record.get(key)).upper() != str(expected[key]).upper():
                failures.append(f"unexpected_{key}")
        for key, tolerance in [("temperature_k", TEMPERATURE_TOLERANCE_K), ("pressure_atm", 1e-6)]:
            value = record.get(key)
            if value is None or not math.isfinite(value) or abs(value - float(expected[key])) > tolerance:
                failures.append(f"unexpected_{key}")
    return failures
