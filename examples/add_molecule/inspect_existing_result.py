"""Independent, minimal energy calculation on a known linkage example."""
from pathlib import Path
import csv
import json

root = Path(__file__).resolve().parents[2]
selected = [r for r in csv.DictReader((root / "data/published/continuum/05_selected_species.csv").open()) if r["analysis_id"] == "full21_b3lyp_pcm"]
inventory = {r["calculation_id"]: r["record_path"] for r in csv.DictReader((root / "data/calculations/manifest.csv").open())}
g = {}
for row in selected:
    record = json.loads((root / inventory[row["calculation_id"]]).read_text())
    g[row["species_id"]] = record["electronic_energy_hartree"] + record["thermal_gibbs_correction_hartree"]

def net(base, nucleotide):
    return (g[nucleotide] + 2 * g["water"] - g[base] - g["ribose"] - g["phosphate"]) * 627.5094740631

example = net("barbituric_acid", "barbituric_acid_clink_ribose_phosphate")
adenine = net("adenine", "adenine_ribose_phosphate")
print(f"C-linked barbituric net reaction: {example:.9f} kcal/mol")
print(f"Relative to adenine: {example - adenine:.9f} kcal/mol")
