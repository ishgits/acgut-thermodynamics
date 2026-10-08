import shutil

import pytest
import yaml

from atcgu.records import read_bundle, within
from atcgu.validation.records import check_record
from atcgu.analysis.reactions import build_reactions
from atcgu.analysis.microsolvation import run_primary


def record(root):
    return read_bundle(root, "data/calculations/calc_557d5df31f5d3868/record.json")


@pytest.mark.parametrize("field,value", [("charge", 1), ("multiplicity", 3), ("basis_set", "STO-3G"), ("solvent", "METHANOL"), ("solvent_model", "SMD"), ("functional", "WB97XD"), ("pressure_atm", 2.0), ("temperature_k", float("nan"))])
def test_wrong_method_or_state_rejected(root, field, value):
    r = record(root); r[field] = value
    method = yaml.safe_load((root / "config/methods.yaml").read_text())["methods"]["1_b3lyp_pcm"]
    assert "unexpected_" + field in check_record(r, method)


def test_missing_and_nonpositive_modes_rejected(root):
    r = record(root); r["frequencies_cm1"] = r["frequencies_cm1"][:-1]
    assert "mode_count_mismatch" in check_record(r)
    r = record(root); r["frequencies_cm1"][0] = 0.0
    assert "nonpositive_frequency" in check_record(r)


def test_path_escape_and_geometry_checksum_guard(root, tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        within(root, "../outside.xyz")
    source = root / "data/calculations/calc_557d5df31f5d3868"
    target = tmp_path / "data/calculations" / source.name
    shutil.copytree(source, target)
    with (target / "optimized.xyz").open("a") as h: h.write("changed\n")
    with pytest.raises(ValueError, match="checksum"):
        read_bundle(tmp_path, str((target / "record.json").relative_to(tmp_path)))


def test_duplicate_species_and_charge_balance(root):
    import pandas as pd
    selected = pd.read_csv(root / "data/published/continuum/05_selected_species.csv")
    selected = selected.loc[selected.analysis_id.eq("full21_b3lyp_pcm")]
    registry = pd.read_csv(root / "config/molecules.csv")
    families = pd.read_csv(root / "config/reaction_families.csv")
    with pytest.raises(ValueError, match="duplicate"):
        build_reactions(pd.concat([selected, selected.iloc[:1]]), registry, families)
    registry.loc[registry.species_id.eq("water"), "charge"] = 1
    rxns, _ = build_reactions(selected, registry, families)
    assert not rxns.balanced.any()


def test_missing_primary_placement_stops_comparison(root, tmp_path):
    import pandas as pd
    config = tmp_path / "config"; config.mkdir()
    df = pd.read_csv(root / "config/microsolvation_candidates.csv").iloc[:-1]
    df.to_csv(config / "microsolvation_candidates.csv", index=False)
    with pytest.raises(ValueError, match="twelve"):
        run_primary(tmp_path, tmp_path, pd.DataFrame())
