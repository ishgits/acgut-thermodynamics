import json
from pathlib import Path
import pytest
from atcgu.analysis.thermochemistry_sensitivity import (
    _harmonic_entropy_per_mode, entropy_floor_shift_hartree, grimme_qrrho_shift_hartree,
)

FIXTURE = json.loads((Path(__file__).parent / 'fixtures/entropy.json').read_text())


@pytest.mark.parametrize('row', FIXTURE['modes'])
def test_independent_mode(row):
    v = float(row['frequency_cm1'])
    assert _harmonic_entropy_per_mode(v, 298) == pytest.approx(float(row['harmonic_entropy_j_mol_k']), abs=1e-12, rel=0)
    for cutoff in (50, 100):
        assert entropy_floor_shift_hartree([v], 298, cutoff) == pytest.approx(float(row[f'floor_{cutoff}']), abs=1e-15, rel=0)
    assert grimme_qrrho_shift_hartree([v], 298) == pytest.approx(float(row['grimme_qrrho_100']), abs=1e-15, rel=0)


def test_independent_spectrum():
    modes = list(map(float, FIXTURE['spectrum_cm1']))
    expected = FIXTURE['spectrum_shifts_hartree']
    for cutoff in (50, 100):
        assert entropy_floor_shift_hartree(modes, 298, cutoff) == pytest.approx(float(expected[f'floor_{cutoff}']), abs=1e-15, rel=0)
    assert grimme_qrrho_shift_hartree(modes, 298) == pytest.approx(float(expected['grimme_qrrho_100']), abs=1e-15, rel=0)
