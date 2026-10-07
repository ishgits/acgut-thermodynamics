from pathlib import Path

import pandas as pd
import pytest

from atcgu.analysis.microsolvation import run_primary
from atcgu.analysis.microsolvation_si import run_si
from atcgu.analysis.thermochemistry_sensitivity import (entropy_floor_shift_hartree, grimme_qrrho_shift_hartree,
                                                        treated_energies)
from atcgu.analysis.workflow import run_continuum
from atcgu.records import read_bundle


@pytest.fixture(scope="module")
def si_outputs(tmp_path_factory):
    root = Path(__file__).resolve().parents[1]
    tables = tmp_path_factory.mktemp("si_tables")
    continuum = run_continuum(root, tables)
    primary = run_primary(root, tables, continuum["05_selected_species"])
    return run_si(root, tables, primary)


def test_floor_shift_zero_when_all_modes_reach_cutoff():
    assert entropy_floor_shift_hartree([100.0, 250.0, 3000.0], 298.0, 100.0) == 0.0
    assert entropy_floor_shift_hartree([20.0, 250.0], 298.0, 100.0) > 0.0


def test_qrrho_shift_vanishes_for_high_frequencies():
    assert abs(grimme_qrrho_shift_hartree([3000.0, 3500.0], 298.0)) < 1e-9
    assert grimme_qrrho_shift_hartree([15.0], 298.0) > grimme_qrrho_shift_hartree([3000.0], 298.0)


def test_treated_energies_guard_record_contract(root):
    record = read_bundle(root, "data/calculations/calc_557d5df31f5d3868/record.json")
    energies = treated_energies(record)
    assert energies["harmonic"] == record["computed_gibbs_hartree"]
    for broken in [{"temperature_k": 310.0}, {"frequency_count": record["frequency_count"] + 1}]:
        with pytest.raises(ValueError):
            treated_energies({**record, **broken})
    with pytest.raises(ValueError, match="positive"):
        treated_energies({**record, "frequencies_cm1": [-5.0, *record["frequencies_cm1"][1:]]})


def test_si_row_counts(si_outputs):
    counts = {name: len(frame) for name, frame in si_outputs.items()}
    assert counts == {"si_03_treatment_energies": 96, "si_05_fixed_reference_placement_scores": 48,
                      "si_06_selected_candidates": 24, "si_07_continuum_rankings": 24,
                      "si_08_selected_fixed_reference_shifts": 24, "si_09_placement_combinations": 256,
                      "si_10_rank_occupancy": 144, "si_thermochemical_summary": 4}
    energies = si_outputs["si_03_treatment_energies"]
    for _, group in energies.groupby("treatment"):
        assert group.record_type.value_counts().to_dict() == {"candidate": 12, "reference": 12}
    vectors = si_outputs["si_09_placement_combinations"]
    assert vectors.groupby("treatment").choice_bits.nunique().eq(64).all()
    occupancy = si_outputs["si_10_rank_occupancy"]
    assert occupancy.groupby(["treatment", "family"]).combination_count.sum().eq(64).all()


def test_fixed_reference_endpoints_and_selections(si_outputs):
    scores = si_outputs["si_05_fixed_reference_placement_scores"]
    harmonic = scores.loc[scores.treatment.eq("harmonic")].set_index(["family", "placement"])
    assert harmonic.fixed_reference_score_kcal_mol["adenine", "splitB"] == pytest.approx(0.668710, abs=1e-6)
    assert harmonic.fixed_reference_score_kcal_mol["adenine", "splitA"] == pytest.approx(4.485969, abs=1e-6)
    assert harmonic.fixed_reference_score_kcal_mol["cytosine", "splitB"] == pytest.approx(0.226373, abs=1e-6)
    assert harmonic.fixed_reference_score_kcal_mol["cytosine", "splitA"] == pytest.approx(2.405896, abs=1e-6)
    expected = {
        "harmonic": {"tap": "splitA"},
        "floor_50": {},
        "floor_100": {},
        "grimme_qrrho_100": {"tap": "splitA"},
    }
    selected = scores.loc[scores.selected_for_family].set_index(["treatment", "family"])
    for treatment, overrides in expected.items():
        for family in scores.family.unique():
            assert selected.placement[treatment, family] == overrides.get(family, "splitB")
    pairs = scores.sort_values(["treatment", "family", "score_hartree", "candidate_id"])
    assert pairs.groupby(["treatment", "family"]).head(1).selected_for_family.all()


def test_fixed_reference_shifts_include_placement_switches(si_outputs):
    selected = si_outputs["si_06_selected_candidates"].set_index(["treatment", "family"])
    shifts = si_outputs["si_08_selected_fixed_reference_shifts"].set_index(["treatment", "family"])
    for (treatment, family), row in shifts.iterrows():
        expected = (selected.fixed_reference_score_kcal_mol[treatment, family]
                    - selected.fixed_reference_score_kcal_mol["harmonic", family])
        assert row.selected_fixed_reference_shift_from_harmonic_kcal_mol == pytest.approx(expected, abs=1e-12)
        assert row.placement_switched_from_harmonic == (row.selected_placement != row.harmonic_selected_placement)


def test_guanine_splitb_is_selected_harmonic_guanine(root, si_outputs):
    selected = si_outputs["si_06_selected_candidates"]
    guanine = selected.loc[selected.treatment.eq("harmonic") & selected.family.eq("guanine")].iloc[0]
    assert guanine.placement == "splitB"
    published = pd.read_csv(root / "data/published/figure1/06_primary_selected_candidates.csv")
    assert guanine.candidate_id == published.set_index("family").candidate_id["guanine"]
