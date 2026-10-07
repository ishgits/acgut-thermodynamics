"""Entropy-only thermochemical treatments applied to released frequency records.

Only the vibrational entropy changes; E, ZPE and H keep their harmonic values.
Frequencies are unscaled, matching Gaussian's harmonic G.
"""
from __future__ import annotations
import math

from ..constants import HARTREE_TO_KCAL_MOL, TEMPERATURE_K, TEMPERATURE_TOLERANCE_K

R_J_MOL_K = 8.31446261815324
PLANCK_J_S = 6.62607015e-34
BOLTZMANN_J_K = 1.380649e-23
LIGHT_CM_S = 2.99792458e10
J_MOL_PER_HARTREE = 4.184 * 1000.0 * HARTREE_TO_KCAL_MOL
GRIMME_BAV_KG_M2 = 1.0e-44

TREATMENTS = ("harmonic", "floor_50", "floor_100", "grimme_qrrho_100")
TREATMENT_LABELS = {
    "harmonic": "Harmonic",
    "floor_50": "50 cm⁻¹ floor",
    "floor_100": "100 cm⁻¹ floor",
    "grimme_qrrho_100": "Grimme qRRHO (100 cm⁻¹)",
}


def _harmonic_entropy_per_mode(frequency_cm1: float, temperature_k: float) -> float:
    x = PLANCK_J_S * LIGHT_CM_S * frequency_cm1 / (BOLTZMANN_J_K * temperature_k)
    return R_J_MOL_K * (x / math.expm1(x) - math.log(-math.expm1(-x)))


def entropy_floor_shift_hartree(
    frequencies_cm1: list[float], temperature_k: float, cutoff_cm1: float
) -> float:
    """Return the Cramer/Truhlar frequency-raising entropy correction."""
    if cutoff_cm1 <= 0:
        return 0.0
    entropy_removed = sum(
        _harmonic_entropy_per_mode(freq, temperature_k)
        - _harmonic_entropy_per_mode(max(freq, cutoff_cm1), temperature_k)
        for freq in frequencies_cm1 if freq > 0
    )
    return temperature_k * entropy_removed / J_MOL_PER_HARTREE


def grimme_qrrho_shift_hartree(
    frequencies_cm1: list[float],
    temperature_k: float,
    cutoff_cm1: float = 100.0,
    bav_kg_m2: float = GRIMME_BAV_KG_M2,
) -> float:
    """Return Grimme's entropy-only qRRHO correction to harmonic G.

    The implementation follows the free-rotor interpolation used by GoodVibes,
    with an exponent of four and the global average moment of inertia from the
    original Grimme treatment. Frequencies are intentionally unscaled so this
    remains a thermochemical-treatment sensitivity check against Gaussian G.
    """
    correction_j_mol = 0.0
    for frequency in frequencies_cm1:
        if frequency <= 0:
            continue
        s_rrho = _harmonic_entropy_per_mode(frequency, temperature_k)
        mu = PLANCK_J_S / (8.0 * math.pi**2 * frequency * LIGHT_CM_S)
        mu_prime = mu * bav_kg_m2 / (mu + bav_kg_m2)
        factor = 8.0 * math.pi**3 * mu_prime * BOLTZMANN_J_K * temperature_k / PLANCK_J_S**2
        s_free_rotor = R_J_MOL_K * (0.5 + 0.5 * math.log(factor))
        damping = 1.0 / (1.0 + (cutoff_cm1 / frequency) ** 4)
        s_qrrho = damping * s_rrho + (1.0 - damping) * s_free_rotor
        correction_j_mol += temperature_k * (s_rrho - s_qrrho)
    return correction_j_mol / J_MOL_PER_HARTREE


def treated_energies(record: dict) -> dict[str, float]:
    """Harmonic G plus each entropy-only shift, from a ``read_bundle`` record (hartree)."""
    frequencies = record["frequencies_cm1"]
    temperature = record["temperature_k"]
    harmonic = record["computed_gibbs_hartree"]
    name = record.get("calculation_id", "record")
    if len(frequencies) != record["frequency_count"]:
        raise ValueError(f"{name}: frequency file and frequency_count disagree")
    if any(f <= 0 for f in frequencies):
        raise ValueError(f"{name}: entropy treatments require all-real positive modes")
    if temperature is None or abs(temperature - TEMPERATURE_K) > TEMPERATURE_TOLERANCE_K:
        raise ValueError(f"{name}: expected {TEMPERATURE_K} K thermochemistry, found {temperature}")
    if harmonic is None or not math.isfinite(harmonic):
        raise ValueError(f"{name}: missing harmonic Gibbs energy")
    return {
        "harmonic": harmonic,
        "floor_50": harmonic + entropy_floor_shift_hartree(frequencies, temperature, 50.0),
        "floor_100": harmonic + entropy_floor_shift_hartree(frequencies, temperature, 100.0),
        "grimme_qrrho_100": harmonic + grimme_qrrho_shift_hartree(frequencies, temperature, 100.0),
    }
