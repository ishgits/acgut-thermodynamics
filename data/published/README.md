# Published results

Each table here is a published result. `atcgu reproduce` regenerates it and compares the fresh values by scientific key (`src/atcgu/validation/regression.py`). Hartree quantities normally use a 1e-9 tolerance; displayed kcal/mol coordinates and shifts use 1e-6. Categorical fields (IDs, choices, orders, ranks and flags) must match exactly. `atcgu audit` independently re-derives the selections, reaction quantities, Figure 1 coordinates, fixed-reference SI scores, 64-vector rankings and occupancy counts with the standard library.

Table numbers follow the analysis pipeline; gaps are intentional.

## `continuum/`: continuum study

Generator: `atcgu.analysis.workflow.run_continuum` (selection in `analysis/selection.py`, reactions in `analysis/reactions.py`).

| Table | Rows | Compared by | Columns compared | Role |
|---|---|---|---|---|
| `01_input_inventory.csv` | 408 | calculation_id | method_id, analysis_role, primary_include, manifest_species_id, local_path, expected/parsed SHA, hash_matches_manifest | Every configured continuum record |
| `02_log_validation.csv` | 408 | calculation_id | status, is_usable; `failures` on non-failed rows only (a fresh run appends `source_calculation_failed` to the three failed rows) | Validation outcome per record |
| `03_thermochemistry.csv` | 408 | calculation_id | E, thermal G correction, computed G, printed G, residual; usable rows only (energies are blank for the three failed rows) | Harmonic energies |
| `04_analysis_candidates.csv` | 503 | analysis_id, calculation_id | candidate rank, selected flag, computed G | All candidates per analysis |
| `05_selected_species.csv` | 206 | analysis_id, species_id | calculation_id, computed G, candidate count | Selected structure per species |
| `06_sampling_coverage.csv` | 81 | analysis_id, reaction_family_id | all candidate counts, required count, complete_required_coverage | Coverage per family |
| `07_reaction_energies.csv` | 186 | analysis, family, reaction | delta_g_kcal_mol | Balanced reaction quantities |
| `08_path_closure.csv` | 21 | analysis_id, reaction_family_id | three residual columns, passes | Two-path closure |
| `09_relative_rankings.csv` | 186 | analysis, family, reaction | relative_g_kcal_mol, rank_within_panel | Adenine-relative values (Figure 1A/B) |
| `10_legacy_vs_transferred_audit.csv` | 45 | method_id, species_id | transferred − legacy, selected file | 45-row matched-start policy audit (diagnostic; never used for selection) |

The canonical-five nucleoside rows in `04`–`09` (`canonical_nucleoside5_*`) are supporting data outside Figure 1.

## `figure1/`: primary two-water comparison (Figure 1C)

Generator: `atcgu.analysis.microsolvation.run_primary`; plotted points from `atcgu.plotting.figure1.render_figure1`.

| Table | Rows | Compared by | Columns compared | Role |
|---|---|---|---|---|
| `02_reference_manifest_and_validation.csv` | 12 | family, role | calculation_id, species_id, expected SHA, computed G | Base and continuum-nucleotide references |
| `03_fixed_reference_convention.csv` | 1 | reference_family | Harmonic continuum-adenine, base, water, conversion and correction terms | Explicit Figure 1C/SI display reference |
| `04_structure_validation.csv` | 12 | candidate_id | structural, classification and displayed-score columns | Water-structure classification |
| `05_primary_candidate_energies.csv` | 12 | candidate_id | G, fixed-reference score, fold_state, intact waters, both H-bond counts | All 12 primary clusters |
| `06_primary_selected_candidates.csv` | 6 | family | candidate_id, relative G, fixed-reference score, rank | Lower-G placement per family |
| `16_panel_c_plotted_data.csv` | 6 | family | continuum and hydrated relative G, selected candidate | Within-condition panel C values |
| `figure1_plotted_data.csv` | 79 | panel, family, series | value (full run only) | Every Figure 1 point |

## `si_microsolvation/`: SI sensitivity

Generator: `atcgu.analysis.microsolvation_si.run_si` (treatments in `analysis/thermochemistry_sensitivity.py`). Both SI plotted-data sets are published under `results/supplement/` and compared in a full run.

| Table | Rows | Compared by | Columns compared | Role |
|---|---|---|---|---|
| `03_treatment_energies.csv` | 96 | treatment, calculation_id | treated G, entropy shift | 24 records × 4 treatments |
| `05_fixed_reference_placement_scores.csv` | 48 | treatment, family, placement | candidate ID, treated cluster-minus-base score, fixed-reference score, selected flag | Four connected placement pairs per family |
| `06_selected_candidates.csv` | 24 | treatment, family | candidate_id, relative G, fixed-reference score, rank | Selection per treatment |
| `07_continuum_rankings.csv` | 24 | treatment, family | relative G, rank | Continuum-only ranking per treatment |
| `08_selected_fixed_reference_shifts.csv` | 24 | treatment, family | selected/harmonic candidates and fixed-reference score difference | Minima-to-minima entropy shift, including switches |
| `09_placement_combinations.csv` | 256 | treatment, choice_bits | rank_order, C − A, pairs retained, A/C top two | 64 placement vectors × 4 treatments |
| `10_rank_occupancy.csv` | 144 | treatment, family, rank | count, selected rank | Counts across enumerated choices |
| `thermochemical_summary.csv` | 4 | treatment | selected order, C − A, pairs retained, unique orders, A/C top-two count | Summary quoted in the SI |
