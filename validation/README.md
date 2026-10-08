# Validation evidence

These files record checks that cannot be rerun from the public release, because the Gaussian and CREST logs are not distributed. Everything else (`atcgu reproduce`, `atcgu audit`, `atcgu validate`, `pytest`) regenerates its own report in a fresh `outputs/` directory.

| File | Coverage | Result |
|---|---:|---|
| `extraction_comparisons.csv` | 420 source logs; 4,620 field comparisons between the production and the independent Gaussian parsers | All agree |
| `extraction_summary.json` | Counts behind the comparison above | No disagreements; 3 source calculations failed and are excluded |
| `construction_checks.csv` | 20 neutral construction recipes | Formula and bond graph match; force-field convergence recorded |
| `refinement_preparation.csv` | 66 CREST ensembles | DFT refinement starts prepared from each |
| `water_preparation.csv` | 12 primary two-water placements | Atom order and alignment against the first source orientation recorded |
| `preparation_coverage.csv` | 66 species | Which starting structures, recipes, seeds, ensembles and records are included |

`../calculation_setup/microsolvation/original_input_manifest.csv` indexes the 12 water-start geometries that were submitted.

See the [workflow](../docs/WORKFLOW.md) to reproduce results and trace values through the released records.
