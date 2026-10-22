import json
import shutil
import pandas as pd
import pytest
from atcgu.validation.regression import compare_table
from atcgu.cli import main


@pytest.mark.parametrize('kind', ['value_mismatch', 'missing_key', 'extra_key', 'duplicate_key', 'nonfinite', 'missing_column'])
def test_structured_comparison(tmp_path, kind):
    expected = pd.DataFrame({'id': ['a', 'b'], 'value': [1., 2.]})
    path = tmp_path / 'expected.csv'; expected.to_csv(path, index=False)
    actual = expected.copy()
    if kind == 'value_mismatch': actual.loc[0, 'value'] = 9
    if kind == 'missing_key': actual = actual.iloc[1:]
    if kind == 'extra_key': actual.loc[2] = ['c', 3.]
    if kind == 'duplicate_key': actual.loc[2] = actual.loc[0]
    if kind == 'nonfinite': actual.loc[0, 'value'] = float('inf')
    if kind == 'missing_column': actual = actual.drop(columns='value')
    result = compare_table(actual, path, ['id'], ['value'])
    assert not result['passed']
    assert kind in {r['issue_type'] for r in result['mismatches']}


def test_cli_failure_preserves_all_reports(root, tmp_path):
    fixture = tmp_path / 'fixture'
    for folder in ['config', 'data', 'src']:
        shutil.copytree(root / folder, fixture / folder)
    path = fixture / 'data/published/continuum/07_reaction_energies.csv'
    frame = pd.read_csv(path); frame.loc[0, 'delta_g_kcal_mol'] += 1
    frame.loc[1, 'reaction_family_id'] = 'missing_family'
    frame.to_csv(path, index=False)
    out = tmp_path / 'failed'
    with pytest.raises(SystemExit) as error:
        main(['--root', str(fixture), 'reproduce', '--no-figures', '--out', str(out)])
    assert error.value.code == 1
    checks = pd.read_csv(out / 'reproduction_checks.csv')
    assert len(checks) == 25 and checks.iloc[-1].passed
    mismatch = pd.read_csv(out / 'reproduction_mismatches.csv')
    assert {'missing_key', 'extra_key', 'value_mismatch'} <= set(mismatch.issue_type)
    assert 'delta_g_kcal_mol' in set(mismatch.column)
    assert not json.loads((out / 'run_metadata.json').read_text())['all_regressions_passed']
    assert 'FAILED' in (out / 'VALIDATION_REPORT.md').read_text()


def test_partial_failure(root, tmp_path, monkeypatch):
    from atcgu import runner
    def broken(*args): raise RuntimeError('deliberate later failure')
    monkeypatch.setattr(runner, 'run_primary', broken)
    out = tmp_path / 'partial'
    with pytest.raises(RuntimeError, match='deliberate'):
        runner.run(root, out, figures=False)
    assert (out / 'tables/09_relative_rankings.csv').exists()
    assert json.loads((out / 'run_metadata.json').read_text())['stage'] == 'primary microsolvation'
    assert (out / 'VALIDATION_REPORT.md').exists()
