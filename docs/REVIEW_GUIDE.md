# Reviewing the results

Begin with a new environment and fresh output directory. Run `reproduce`, `audit`, `validate` and the tests.

## Trace a plotted value

For a Figure 1A/B point, locate panel/family/series in `figure1_plotted_data.csv`, identify the corresponding analysis/reaction in `09_relative_rankings.csv`, inspect the balanced coefficients and absolute value in `07_reaction_energies.csv`, then follow each species through `05_selected_species.csv` to its calculation bundle. Recalculate E + thermal G correction and sum the signed coefficients. Subtract that method's reference-family (adenine) reaction. Inspect all alternatives in `04_analysis_candidates.csv` and confirm the declared selection policy.

For Figure 1C, follow the selected **and other** placement through `micro_05_primary_candidate_energies.csv`. Inspect the full XYZ including both waters, complete frequency list, solute identity, fragment composition and contacts. Recompute the displayed fixed-reference value using `micro_03_fixed_reference_convention.csv` and the equation in [Methods](METHODS.md), including harmonic continuum water and the −8.539633 kcal/mol convention offset. Do not confuse the within-hydrated-condition ranking table with the actual common-zero plotted coordinates.

For the SI, verify all 48 rows of `05_fixed_reference_placement_scores.csv`: only each cluster and corresponding base change treatment. Check that the filled endpoint is the lower score in every pair and recompute `08_selected_fixed_reference_shifts.csv` from independently selected minima. Review `09_placement_combinations.csv` separately from the placement-range figure, then confirm every family/treatment row in `10_rank_occupancy.csv` sums to 64. The occupancy counts enumerate choices and do not estimate physical probabilities.

## Review the chemistry and scope

- Compare base/product SDF bonds, linkage indices, intended tautomers and absolute stereochemistry with the optimized XYZs. The graph fingerprint does not prove stereochemistry.
- Confirm the neutral-singlet choices, functional/basis/solvent, thermochemistry temperature/pressure and Gaussian revisions in each selected record.
- Confirm full primary-pool selection versus fixed transferred products at alternate methods. Diagnostic-only alternative products (`diagnostic_only_legacy_product`) are never selected.
- Inspect all three failed continuum candidates and their exclusion reasons. Missing data are explicit, not zero-filled.
- Review every primary water arrangement, unconstrained final geometry and contact classification. Two starts are a sensitivity experiment, not a convergence claim.
- Original water-start XYZ geometries are indexed in `calculation_setup/microsolvation/original_input_manifest.csv`.

## Preparation evidence

`validation/` keeps only the evidence that needs the unreleased logs: `extraction_comparisons.csv` (4,620 checks from two separate extraction implementations for 420 records), `extraction_summary.json`, `construction_checks.csv` (all 20 rebuilt neutral products against the registry's formula and graph, with convergence), `refinement_preparation.csv`, `water_preparation.csv` and `preparation_coverage.csv`. See [its README](../validation/README.md).

Regenerate the result checks yourself: `atcgu reproduce` writes `reproduction_checks.csv`, `atcgu audit` writes `independent_arithmetic_checks.csv`, and `atcgu validate` writes `record_checks.csv`.
