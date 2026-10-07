"""Complete released-bundle integrity and scientific record audit."""
from pathlib import Path
import csv
import json

from ..paths import check_output, fresh
from ..records import digest, read_bundle, within
from .records import check_record


def validate(root: Path, output: Path) -> Path:
    fresh(check_output(root, output))
    checks, record_status = [], {}
    for row in csv.DictReader((root / "data/calculations/manifest.csv").open()):
        record_status[row["calculation_id"]] = row["status"]
        record = read_bundle(root, row["record_path"], row["record_sha256"])
        errors = check_record(record)
        expected_failure = row["status"] == "failed" and errors == ["source_calculation_failed"]
        checks.append({"item": row["calculation_id"], "passed": not errors or expected_failure,
                       "status": "documented_source_failure" if expected_failure else "usable",
                       "details": ";".join(errors)})
    for row in csv.DictReader((root / "data/source/manifest.csv").open()):
        checks.append({"item": row["path"], "passed": digest(within(root, row["path"])) == row["sha256"], "status": "source_integrity", "details": ""})
    for row in csv.DictReader((root / "config/sampling.csv").open()):
        for field in ["seed", "ensemble"]:
            if row[field + "_path"]:
                checks.append({"item": row[field + "_path"], "passed": digest(within(root, row[field + "_path"])) == row[field + "_sha256"], "status": "sampling_integrity", "details": ""})
    # Re-validate method, state, identity and cohort membership in the continuum workflow;
    # only records that failed at the source may be unusable.
    from ..analysis.workflow import load_inventory
    identity, sets, registry, families, panels = load_inventory(root)
    unexpected = identity.loc[~identity.is_usable & identity.calculation_id.map(record_status).ne("failed")]
    checks.append({"item": "continuum_method_state_identity", "passed": unexpected.empty, "status": "scientific_identity", "details": ";".join(unexpected.calculation_id)})
    from ..analysis.microsolvation import run_primary
    from ..analysis.selection import _analysis_inputs
    _, selected, _ = _analysis_inputs(identity, sets, registry, families, panels)
    (output / "primary_tables").mkdir()
    try:
        run_primary(root, output / "primary_tables", selected)
        passed, details = True, "12 primary candidates; two intact waters; both contacts and dry-reference identity checked"
    except (ValueError, KeyError) as exc:
        passed, details = False, str(exc)
    checks.append({"item": "primary_water_structure_and_references", "passed": passed, "status": "scientific_structure", "details": details})
    with (output / "record_checks.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(checks[0])); w.writeheader(); w.writerows(checks)
    summary = {"checks": len(checks), "failed_checks": sum(not r["passed"] for r in checks),
               "documented_failed_source_calculations": sum(r["status"] == "documented_source_failure" for r in checks)}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if summary["failed_checks"]:
        raise ValueError("Record validation failed; inspect record_checks.csv")
    return output
