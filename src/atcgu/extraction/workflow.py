"""Extract supplied Gaussian logs into portable bundles."""
from pathlib import Path
import csv
import json
import math

from .gaussian import parse_log
from ..records import write_bundle, digest
from ..paths import fresh
from ..validation.independent_extraction import parse as independent_parse

COMPARED_FIELDS = ["electronic_energy_hartree", "thermal_gibbs_correction_hartree", "printed_gibbs_hartree",
                   "charge", "multiplicity", "formula", "temperature_k", "pressure_atm"]


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def extract(paths: list[Path], output: Path) -> Path:
    """Parse and cross-check every log first; write bundles only when every check passes."""
    fresh(output)
    records, checks, seen = [], [], set()
    for path in paths:
        record = parse_log(path)
        record.pop("log_path", None)
        sha = record["sha256"]
        if sha in seen:
            raise ValueError("Duplicate source bytes in extraction inputs")
        seen.add(sha)
        record.update(calculation_id="calc_" + sha[:16], source_log_sha256=sha,
                      source_log_name=path.name,
                      frequency_job_normally_terminated="frequency_job_not_normally_terminated" not in record["failures"],
                      units={"coordinates": "angstrom", "frequencies": "cm^-1", "energies": "hartree per particle",
                             "temperature": "K", "pressure": "atm", "entropy": "cal mol^-1 K^-1"})
        independent = independent_parse(path)
        for key in COMPARED_FIELDS:
            a, b = record.get(key), independent.get(key)
            passed = a == b if a is None or b is None or isinstance(a, str) else math.isclose(a, b, abs_tol=1e-10, rel_tol=0)
            checks.append({"calculation_id": record["calculation_id"], "field": key, "passed": passed})
        checks.append({"calculation_id": record["calculation_id"], "field": "full_frequencies",
                       "passed": record["frequencies_cm1"] == independent["frequencies"]})
        for key in ["zpe_correction_hartree", "thermal_energy_correction_hartree", "thermal_enthalpy_correction_hartree",
                    "electronic_plus_zpe_hartree", "H_hartree", "entropy_cal_mol_k", "source_lines"]:
            record[key] = independent.get(key)
        records.append(record)
    if not records:
        raise ValueError("No input logs")
    _write_csv(output / "independent_extraction_checks.csv", checks)
    if not all(r["passed"] for r in checks):
        _write_csv(output / "extraction_inventory.csv",
                   [{"calculation_id": r["calculation_id"], "record_path": "", "record_sha256": "",
                     "source_log_sha256": r["sha256"], "status": r["status"], "failures": r["failures"]} for r in records])
        raise ValueError("Independent extraction disagreement; no bundles were written. "
                         "Inspect independent_extraction_checks.csv")
    inventory, manifest = [], []
    for record in records:
        relative = write_bundle(output, record)
        written = json.loads((output / relative).read_text())
        sha = record["source_log_sha256"]
        inventory.append({"calculation_id": record["calculation_id"], "record_path": relative,
                          "record_sha256": digest(output / relative), "source_log_sha256": sha,
                          "status": record["status"], "failures": record["failures"]})
        manifest.append({"calculation_id": written["calculation_id"], "record_path": relative,
                         "record_sha256": digest(output / relative), "source_log_sha256": sha,
                         "source_log_name": written["source_log_name"], "status": written["status"],
                         "frequency_count": written["frequency_count"], "atom_count": written["atom_count"],
                         "geometry_file": written["geometry_file"] or ""})
    # manifest_rows.csv has exactly the columns of data/calculations/manifest.csv, in order.
    _write_csv(output / "extraction_inventory.csv", inventory)
    _write_csv(output / "manifest_rows.csv", manifest)
    (output / "README.md").write_text("# New extracted records\n\nInspect `extraction_inventory.csv` and `independent_extraction_checks.csv` before use. "
                                      "`manifest_rows.csv` can be appended to `data/calculations/manifest.csv` after review; "
                                      "a valid row does not approve the chemical identity. See docs/EXTENDING.md.\n")
    return output
