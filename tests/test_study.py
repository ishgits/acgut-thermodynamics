import csv
import json

import pytest

from atcgu.runner import run
from atcgu.validation.arithmetic import audit


def test_complete_numerical_reproduction(root, tmp_path):
    out = run(root, tmp_path / "reproduction", figures=False)
    checks = list(csv.DictReader((out / "reproduction_checks.csv").open()))
    assert len(checks) == 25 and all(r["passed"] == "True" for r in checks)
    with pytest.raises(FileExistsError):
        run(root, out, figures=False)


def test_independent_standard_library_audit(root, tmp_path):
    out = audit(root, tmp_path / "audit")
    summary = json.loads((out / "summary.json").read_text())
    assert summary["failed_checks"] == 0
    assert summary["reaction_quantities"] == 186
    assert summary["figure_points"] == 79
    assert summary["si_placement_vectors"] == 256



def test_neutral_linkage_builder_reproduces_defined_chemistry(root, tmp_path):
    pytest.importorskip("rdkit")
    import csv
    from atcgu.building.molecules import build_from_registry
    from atcgu.records import read_xyz
    from atcgu.validation.records import fingerprint
    output = build_from_registry(root, root / "config/construction.csv", tmp_path / "build")
    infos = json.loads((output / "construction_manifest.json").read_text())
    expected = {r["species_id"]: r["expected_formula"] for r in csv.DictReader((root / "config/construction.csv").open())}
    assert len(infos) == len(expected)
    for info in infos:
        assert info["charge"] == 0 and info["formula"] == expected[info["species_id"]]
        z, xyz = read_xyz(output / (info["species_id"] + ".xyz"))
        geometry = [[int(a), *point.tolist()] for a, point in zip(z, xyz)]
        registry = {r["species_id"]: r for r in csv.DictReader((root / "config/molecules.csv").open())}
        assert fingerprint(geometry) == registry[info["species_id"]]["connectivity_fingerprint"]



def test_microsolvation_298_k_common_shift(root, tmp_path):
    import math
    import pandas as pd

    out = run(root, tmp_path / "water_correction", figures=False)
    convention = json.loads((root / "config/study.json").read_text())["panel_c_display_convention"]
    assert convention["temperature_k"] == 298.0
    assert convention["gas_to_solution_factor"] == pytest.approx(24.4654 * 298 / 298.15)
    expected_offset = -2 * 0.0019872041 * 298 * math.log((24.4654 * 298 / 298.15) * 55.34)
    energy = pd.read_csv(out / "tables/micro_05_primary_candidate_energies.csv", float_precision="round_trip")
    refs = pd.read_csv(out / "tables/micro_02_reference_manifest_and_validation.csv", float_precision="round_trip")
    dry = refs.loc[refs.role.eq("continuum_nucleotide")].set_index("family").computed_gibbs_hartree
    selected = pd.read_csv(out / "tables/05_selected_species.csv", float_precision="round_trip")
    water = float(selected.loc[selected.analysis_id.eq("full21_b3lyp_pcm") & selected.species_id.eq("water"), "computed_gibbs_hartree"].iloc[0])
    raw = (energy.g_hartree - energy.family.map(dry) - 2 * water) * 627.509474
    assert (energy.association_shift_kcal_mol - raw).tolist() == pytest.approx([expected_offset] * 12, abs=1e-10)
    adenine = energy.loc[energy.family.eq("adenine") & energy.selected_for_family, "association_shift_kcal_mol"].iloc[0]
    assert adenine == pytest.approx(0.6687099313612619, abs=1e-9)
