import csv
import json

import pytest

from atcgu.runner import run
from atcgu.validation.arithmetic import audit
from atcgu.sampling.prepare import prepare_refinement, sampling_policy
from atcgu.building.gaussian import prepare_study_inputs


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


def test_refinement_preserves_ensemble_conformer_numbers(root, tmp_path):
    seed = root / "data/sampling/adenine_ribose_phosphate/seed.xyz"
    ensemble = seed.with_name("crest_conformers.xyz")
    output = prepare_refinement(seed, ensemble, "nucleotide", tmp_path / "refine", sampling_policy(root))
    info = json.loads((output / "refinement_manifest.json").read_text())
    assert info["adjacent_energy_inversions"] == 1
    assert {c["candidate"] for c in info["candidates"]} == {"dft_seed", "conf1", "conf13", "conf22", "conf52", "conf75", "conf113"}


def test_neutral_linkage_builder_reproduces_defined_chemistry(root, tmp_path):
    pytest.importorskip("rdkit")
    from atcgu.building.molecules import build_from_registry
    from atcgu.records import read_xyz
    from atcgu.validation.records import fingerprint
    output = build_from_registry(root, root / "examples/add_molecule/construction.csv", tmp_path / "build")
    info = json.loads((output / "construction_manifest.json").read_text())[0]
    assert info["charge"] == 0 and info["formula"] == "C9H13N2O10P"
    z, xyz = read_xyz(output / (info["species_id"] + ".xyz"))
    geometry = [[int(a), *point.tolist()] for a, point in zip(z, xyz)]
    registry = {r["species_id"]: r for r in csv.DictReader((root / "config/molecules.csv").open())}
    assert fingerprint(geometry) == registry[info["species_id"]]["connectivity_fingerprint"]


def test_study_reruns_use_optimized_geometries_and_plain_opt(root, tmp_path):
    import numpy as np
    from atcgu.records import read_xyz, SYMBOLS, digest

    out = prepare_study_inputs(root, tmp_path / "gaussian")
    manifest = json.loads((out / "input_manifest.json").read_text())
    assert len(manifest) == 420
    assert sum(row["status"] == "not_prepared" for row in manifest) == 3
    for row in manifest:
        if row["status"] != "prepared":
            continue
        path = out / (row["calculation_id"] + ".com")
        assert path.read_bytes() == (root / "calculation_setup/gaussian/jobs" / path.name).read_bytes()
        text = path.read_text()
        routes = [line for line in text.splitlines() if line.startswith("#")]
        assert routes[0].startswith("# opt ") and "opt=" not in routes[0]
        assert "temperature=298 geom=allchk guess=read" in routes[1]
        geometry = root / row["coordinate_source"]
        assert geometry.name == "optimized.xyz" and digest(geometry) == row["geometry_sha256"]
        z, xyz = read_xyz(geometry)
        supplied = text.split("\n\n")[2].splitlines()[1:]
        assert [line.split()[0] for line in supplied] == [SYMBOLS[int(a)] for a in z]
        actual = np.array([[float(v) for v in line.split()[1:]] for line in supplied])
        np.testing.assert_allclose(actual, xyz, atol=5e-7, rtol=0)


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
