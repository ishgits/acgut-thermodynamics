import csv
import json
import shutil
import numpy as np
import pytest
from atcgu.sampling import selection as s
from atcgu.validation.sampling import check_run, validate_sampling


def frames(radii):
    return [(0, ['C', 'C'], np.array([[-r, 0., 0.], [r, 0., 0.]])) for r in radii]


def test_first_above_window_stops():
    assert s.select_energy(frames([1, 2, 3]), np.array([0, 4, 1]), np.array([True, True]), {'window': 3, 'n_keep': 5}, .5) == [0]


def test_rmsd_boundary():
    f = frames([1, 1.5, 1.499999])
    mask = np.array([True, True])
    assert s.kabsch_rmsd(f[0][2], f[1][2], mask) == .5
    assert s.select_energy(f, [0, 1, 2], mask, {'window': 3, 'n_keep': 5}, .5) == [0, 1]


def test_rg_final_edge_and_across_bin_dedup():
    mask = np.array([True, True])
    policy = {'window': 6, 'n_bins': 3, 'per_bin': 2}
    assert s.select_rg_stratified(frames([1, 2, 3]), np.array([0, 1, 6]), mask, policy, 3, .5) == [0, 1, 2]
    assert s.select_rg_stratified(frames([1, 1.3, 1.6]), np.array([0, 1, 2]), mask, policy, 1.6, .5) == [0, 2]


def test_disposable_manifest_mismatch(root, tmp_path):
    for folder in ['config', 'data', 'src', 'calculation_setup/crest/refinement']:
        shutil.copytree(root / folder, tmp_path / folder)
    row = next(csv.DictReader((tmp_path / 'config/sampling.csv').open()))
    path = tmp_path / 'calculation_setup/crest/refinement' / row['sampling_run_id'] / 'refinement_manifest.json'
    manifest = json.loads(path.read_text()); manifest['candidates'][1]['candidate'] = 'conf_wrong'
    path.write_text(json.dumps(manifest))
    checks = validate_sampling(tmp_path, tmp_path)
    assert not checks[0]['passed'] and 'candidates_match' in checks[0]['details']
    assert all(c['passed'] for c in checks[1:])
    from atcgu.cli import main
    with pytest.raises(SystemExit) as error:
        main(['--root', str(tmp_path), 'validate', '--out', str(tmp_path / 'outputs/failed-validation')])
    assert error.value.code == 1
    summary = json.loads((tmp_path / 'outputs/failed-validation/summary.json').read_text())
    assert summary['failed_checks'] == 1


@pytest.mark.parametrize('bad', ['nan', 'order', 'empty'])
def test_malformed_ensemble(root, tmp_path, bad):
    row = next(csv.DictReader((root / 'config/sampling.csv').open()))
    row['seed_path'] = str(root / row['seed_path'])
    ensemble = tmp_path / 'ensemble.xyz'
    ensemble.write_text('' if bad == 'empty' else '1\n' + ('nan' if bad == 'nan' else '0') + '\nC 0 0 0\n')
    row['ensemble_path'] = str(ensemble)
    with pytest.raises(ValueError):
        check_run(root, row, json.loads((root / 'config/study.json').read_text()))
