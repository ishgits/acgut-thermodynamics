"""Panel, schema, extension-output and output-path contracts."""
import csv
import json
import shutil

import pandas as pd
import pytest

from atcgu.analysis.workflow import run_continuum
from atcgu.extraction.workflow import extract
from atcgu.paths import check_output
from atcgu.runner import run


def _copy_config(root, tmp_path, drop_family):
    shutil.copytree(root / "config", tmp_path / "config")
    (tmp_path / "data").symlink_to(root / "data", target_is_directory=True)
    path = tmp_path / "config/reaction_families.csv"
    families = pd.read_csv(path, dtype=str, keep_default_na=False)
    families.loc[families.reaction_family_id.ne(drop_family)].to_csv(path, index=False)
    return tmp_path


@pytest.mark.parametrize("family,message", [("hypoxanthine", "expects 21 families"), ("xanthine", "absent from config/reaction_families.csv")])
def test_shrunken_family_panel_fails_before_tables(root, tmp_path, family, message):
    checkout = _copy_config(root, tmp_path / "checkout", family)
    tables = tmp_path / "tables"
    tables.mkdir()
    with pytest.raises(ValueError, match=message):
        run_continuum(checkout, tables)
    assert not any(tables.iterdir())


def test_all_records_and_fresh_extraction_match_schema(root, water_log, tmp_path):
    jsonschema = pytest.importorskip("jsonschema")
    validator = jsonschema.Draft202012Validator(json.loads((root / "data/schemas/calculation_record.schema.json").read_text()))
    rows = list(csv.DictReader((root / "data/calculations/manifest.csv").open()))
    assert len(rows) == 420
    for row in rows:
        errors = [e.message for e in validator.iter_errors(json.loads((root / row["record_path"]).read_text()))]
        assert not errors, (row["calculation_id"], errors)
    output = extract([water_log], tmp_path / "extracted")
    record = json.loads(next(output.rglob("record.json")).read_text())
    assert not list(validator.iter_errors(record))


def test_extraction_writes_manifest_rows_in_release_shape(root, water_log, tmp_path):
    output = extract([water_log], tmp_path / "extracted")
    with (root / "data/calculations/manifest.csv").open() as h:
        expected = next(csv.reader(h))
    with (output / "manifest_rows.csv").open() as h:
        rows = list(csv.DictReader(h))
    assert list(rows[0]) == expected and len(rows) == 1
    assert (output / rows[0]["record_path"]).is_file()


def test_output_path_guard(root, tmp_path):
    for inside in ["config", "data/published", "src", "results", ".", "docs"]:
        with pytest.raises(ValueError, match="outputs"):
            check_output(root, root / inside / "run")
    assert check_output(root, root / "outputs/review") == (root / "outputs/review").resolve()
    assert check_output(root, tmp_path / "anywhere") == (tmp_path / "anywhere").resolve()
    with pytest.raises(ValueError, match="outputs"):
        run(root, root / "results/reproduction", figures=False)


def test_water_counts_are_distinct_waters_not_contacts(root):
    from atcgu.analysis.structure import analyze_microsolvation_structure
    from atcgu.records import read_xyz
    from atcgu.structure import analyze_structure
    row = next(r for r in csv.DictReader((root / "config/microsolvation_candidates.csv").open())
               if r["family"] == "guanine" and r["placement"] == "splitA")
    path = root / row["record_path"].replace("record.json", "optimized.xyz")
    assert analyze_microsolvation_structure(path)["phosphate_to_water_hb_count"] == 1
    assert analyze_structure(*read_xyz(path))["phosphate_water_contact_count"] == 2
