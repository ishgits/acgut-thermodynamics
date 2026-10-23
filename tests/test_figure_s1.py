import pandas as pd
import pytest
import yaml
from atcgu.plotting.decomposition import load_rankings, render_decomposition
from atcgu.validation.regression import compare_table


def test_s1_points(root, tmp_path):
    config = yaml.safe_load((root / 'config/figures.yaml').read_text())
    render_decomposition(root / 'data/published/continuum', tmp_path, 80, config)
    actual = pd.read_csv(tmp_path / 'figure_s1_plotted_data.csv')
    assert len(actual) == 63
    result = compare_table(actual, root / 'data/published/supplement/figure_s1_plotted_data.csv',
                           ['panel', 'analysis_id', 'reaction_id', 'reaction_family_id'], ['display_order', 'canonical', 'value_kcal_mol'], 1e-6)
    assert result['passed']
    actual.loc[0, 'display_order'] = 9
    assert not compare_table(actual, root / 'data/published/supplement/figure_s1_plotted_data.csv', ['panel', 'reaction_family_id'], ['display_order'])['passed']


@pytest.mark.parametrize('mutation', ['duplicate', 'panel', 'closure'])
def test_s1_guards(root, tmp_path, mutation):
    data = pd.read_csv(root / 'data/published/continuum/09_relative_rankings.csv')
    if mutation == 'duplicate': data = pd.concat([data, data.iloc[[0]]])
    if mutation == 'panel': data = data.loc[data.reaction_id.ne('P1a_nucleoside_formation')]
    if mutation == 'closure': data.loc[0, 'relative_g_kcal_mol'] += 1
    path = tmp_path / 'rankings.csv'; data.to_csv(path, index=False)
    with pytest.raises(ValueError): load_rankings(path)
