import pytest

from atcgu.extraction.gaussian import parse_log, _route_text
from atcgu.extraction.workflow import extract
from atcgu.validation.independent_extraction import parse as independent
from atcgu.records import read_bundle


def test_complete_scoped_job(water_log):
    r = parse_log(water_log)
    assert r["status"] == "valid"
    assert r["frequency_count"] == 3
    assert abs(r["computed_gibbs_hartree"] + 76.395) < 1e-12
    assert r["geometry"] == independent(water_log)["geometry"]


def test_trailing_single_point_does_not_replace_frequency_energy(water_log):
    with water_log.open("a") as h:
        h.write("# sp wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water)\n --------------------\n SCF Done: E(RWB97XD) = -999.0 A.U.\n Normal termination of Gaussian\n")
    r = parse_log(water_log)
    assert r["electronic_energy_hartree"] == -76.4
    assert r["functional"] == "B3LYP"
    assert r["solvent_model"] == "PCM"
    assert independent(water_log)["electronic_energy_hartree"] == -76.4


def test_unfinished_final_frequency_job_is_not_replaced_by_previous_success(water_log):
    with water_log.open("a") as h:
        h.write("# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)\n --------------------\n Error termination\n")
    r = parse_log(water_log)
    assert r["status"] == "failed"
    assert r["computed_gibbs_hartree"] is None


def test_complete_frequency_list_keeps_imaginary_modes(water_log):
    water_log.write_text(water_log.read_text().replace("1600.0", "-100.0"))
    r = parse_log(water_log)
    assert r["status"] == "failed"
    assert r["frequencies_cm1"] == [-100.0, 3500.0, 3600.0]


def test_input_orientation_supported(water_log):
    water_log.write_text(water_log.read_text().replace("Standard orientation", "Input orientation"))
    assert parse_log(water_log)["status"] == "valid"


def test_new_log_extracts_to_portable_bundle_without_copying_source(water_log, tmp_path):
    output = tmp_path / "records"
    extract([water_log], output)
    path = next((output / "data/calculations").glob("*/record.json"))
    r = read_bundle(output, str(path.relative_to(output)))
    assert r["frequencies_cm1"] == [1600.0, 3500.0, 3600.0]
    assert len(r["geometry"]) == 3
    assert not list(output.rglob("*.log"))
    assert not {"route", "route_optimization", "route_frequency"}.intersection(r)


@pytest.mark.parametrize("boundary", ["between", "inside", "leading_space"])
def test_wrapped_route_preserves_keyword_boundaries(water_log, boundary):
    route = "# freq b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 Geom=AllChk Guess=Read"
    split = route.index("Geom")
    if boundary == "inside":
        split += 5
    elif boundary == "leading_space":
        split -= 1
    echo = " " + route[:split] + "\n " + route[split:] + "\n --------------------\n"
    assert _route_text(echo, 0) == route
    text = water_log.read_text()
    original = text[text.index("# freq"):text.index(" --------------------", text.index("# freq"))]
    water_log.write_text(text.replace(original, echo.split(" --------------------")[0]))
    record = parse_log(water_log)
    assert record["status"] == "valid"
    assert record["functional"] == "B3LYP"
    assert record["solvent_model"] == "SMD"
    assert record["temperature_k"] == 298.0


def test_independent_disagreement_writes_no_bundles(water_log, tmp_path, monkeypatch):
    import atcgu.extraction.workflow as workflow

    def disagreeing(path):
        values = independent(path)
        values["electronic_energy_hartree"] += 1.0
        return values

    monkeypatch.setattr(workflow, "independent_parse", disagreeing)
    output = tmp_path / "records"
    with pytest.raises(ValueError, match="no bundles were written"):
        extract([water_log], output)
    assert not (output / "data/calculations").exists()
    assert (output / "independent_extraction_checks.csv").is_file()
    assert (output / "extraction_inventory.csv").is_file()
    assert not (output / "manifest_rows.csv").exists()
