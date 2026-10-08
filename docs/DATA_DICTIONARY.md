# Data dictionary

CSV files are UTF-8 with headers. Empty cells mean missing/not applicable; they are never replaced with zero. Paths are relative to the repository root.

## IDs and configuration

| Field/location | Meaning |
|---|---|
| `species_id` | Stable modeled chemical species; includes linkage identity for constructed products |
| `reaction_family_id` | Free-base/product family; two families can share a free base |
| `method_id` | Functional/continuum treatment in `methods.yaml` |
| `analysis_id` | Candidate cohort, species selection and reaction scope in `analysis_sets.csv` |
| `calculation_id` | `calc_` plus the first 16 hex characters of the source-log SHA-256; full source SHA is retained and IDs are unique in the release |
| `config/calculation_manifest.csv` | All 405 continuum records, source hashes, roles, identity claims, method and primary inclusion flags |
| `config/microsolvation_candidates.csv` | The 12 primary cluster candidates and expected fragment/contact criteria |
| `config/microsolvation_references.csv` | Base and selected continuum-nucleotide reference for each hydrated family; `folded_start` marks the nucleotides whose water starts keep the direct fold |
| `config/construction.csv` | Explicit-H source SDFs, **zero-based** attachment/OH atom indices, seed and expected product formula |
| `config/sampling.csv` | Sampling-run ID, mapped species, moiety, seed/ensemble paths and hashes; sampling-run IDs need not be unique species IDs |
| `config/reaction_stoichiometry.csv` | Integer signed coefficients: negative reactants, positive products |
| `config/study.json`, `figures.yaml` | Scope, reference family, family panels, sampling selection rules, numerical constants and figure conventions |

`primary_include` describes cohort eligibility; `is_usable` in the output tables is computed from record validation. `analysis_role=diagnostic_only_legacy_product` marks a diagnostic-only alternative structure that is never selected. `required_transferred_nucleotide` and `transferred_nucleoside` indicate fixed alternate-method structures; `method_matched_reference` identifies their reference calculations. `file_name` identifies the source log by name and `local_path` is the public record path; no log bytes are bundled.

## Calculation bundles

`data/calculations/manifest.csv` is the shared index. Each record has `record.json`, `frequencies.csv` and its available geometry. Every valid record includes `optimized.xyz`; failed records use `incomplete.xyz` when coordinates were available. `starting.xyz` is the first printed source orientation, not a claimed optimized result. The manifest and `record.files` provide SHA-256 checksums.

| Record field | Units/definition |
|---|---|
| `schema_version` | Current contract `1.0` |
| `status`, `failures`, `warnings` | Explicit extraction status and reasons; failed candidates cannot be selected |
| `functional`, `basis_set`, `solvent_model`, `solvent` | Parsed from the selected frequency job's route |
| `charge`, `multiplicity` | Gaussian molecular state |
| `temperature_k`, `pressure_atm` | Thermochemistry state in K and atm |
| `electronic_energy_hartree` | Final SCF energy in the selected frequency job |
| `thermal_gibbs_correction_hartree` | Printed thermal correction to G |
| `computed_gibbs_hartree` | Electronic energy plus correction |
| `printed_gibbs_hartree`, `gibbs_residual_hartree` | Gaussian printed total and computed-minus-printed difference |
| `zpe_correction_hartree`, `thermal_energy_correction_hartree`, `thermal_enthalpy_correction_hartree` | Complete printed scalar thermochemical corrections |
| `electronic_plus_zpe_hartree`, `H_hartree` | Printed energy plus ZPE and energy plus thermal enthalpy |
| `entropy_cal_mol_k` | Printed total entropy, cal mol^-1 K^-1 |
| `frequency_count`, `expected_frequency_count` | Full mode count; 3N−6 nonlinear, 3N−5 linear, zero monatomic |
| `imaginary_frequency_count`, `minimum_frequency_cm1` | Counts of negative modes and minimum printed mode |
| `source_log_sha256`, `source_log_name` | Original byte hash and source filename |
| `source_lines`, `thermo_segment_first_line`, `thermo_segment_last_line` | One-based source evidence locations/job boundaries |
| `gaussian_version` | Recorded Gaussian software revision |
| `atom_indices`, `geometry_fragments`, `distance_inferred_edges` | **One-based** source Cartesian ordering/derived connectivity; no independent chemical bond-order assignment |
| `geometry_precision_angstrom` | Source output precision, 10^-6 Å |

`frequencies.csv` holds **every printed vibrational mode**, one row per mode, with contiguous one-based `mode_index` and `frequency_cm1` (cm^-1). Negative frequencies are not removed. Failed sources may have incomplete or empty spectra; the empty file still has a header and the record explains the failure. XYZ coordinates are Å in source atom order. Water atoms are retained in the full cluster geometry.

Each `record.json` is the authoritative record; the loader recomputes facts from its linked geometry and complete spectrum. Water preparation manifests use **zero-based** starting-geometry atom indices, as do construction SDF choices; record fragments/edges use the explicitly documented one-based source ordering.

Original water-start XYZ coordinates are recorded to 10^-9 Å.

## Numerical results

Energies named `_hartree` are per particle; reaction/relative quantities named `_kcal_mol` are kcal/mol. `computed_gibbs_hartree` always uses E + thermal correction. `delta_g_kcal_mol` is products minus reactants. `relative_g_kcal_mol` subtracts the reference-family (adenine) reaction in the same scope. Rank 1 is lowest, with minimum ranks for exact ties (SI tables break exact ties by row order, `rank(method="first")`). Candidate ranks use energy and original source filename.

Figure 1C has two valid reference frames. `micro_16_panel_c_plotted_data.csv` retains **within-condition** adenine-relative hydrated values for ranking. The actual plotted common-zero values are in `figures/figure1_plotted_data.csv` and the candidate table's explicit `fixed_reference_score_kcal_mol`. The legacy-compatible `common_reference_value_kcal_mol` is retained as an equivalent diagnostic. These differ from the within-condition values by the hydrated adenine displacement. `micro_03_fixed_reference_convention.csv` records every fixed term used in the displayed score.

Structural contact cutoffs are H···A ≤ 2.20 Å, D···A ≤ 3.20 Å, angle D–H···A ≥ 130°. Near contacts use 2.50 Å, 3.50 Å and 120°. `fold_state` and water-fragment/contact counts are derived classifications, not annotations copied from a filename. `phosphate_to_water_hb_count` is the number of distinct waters with ≥1 qualifying H-bond from a phosphate OH; `water_to_base_hb_count` is the number of distinct waters with ≥1 qualifying H-bond to a base acceptor. The placement manifests under `calculation_setup/microsolvation/` instead count individual contacts (`phosphate_water_contact_count`, `water_base_contact_count`).

JSON schema: [calculation_record.schema.json](../data/schemas/calculation_record.schema.json). Schema documents the contract; runtime loaders additionally check actual files, mode order/count, integrity, geometry, method and state.

## SI microsolvation tables

`data/published/si_microsolvation/` holds the published SI tables; `acgut reproduce` writes the same files to `tables/si/`. `treatment` is one of `harmonic`, `floor_50`, `floor_100` or `grimme_qrrho_100`; see the manuscript Methods for the computational details.

| Table | Key | Fields |
|---|---|---|
| `03_treatment_energies.csv` | treatment, calculation_id | `record_type` (candidate/reference), `role` (`two_water_cluster`, `base`, `continuum_nucleotide`), `harmonic_g_hartree`, `entropy_shift_hartree` (treated − harmonic), `treated_g_hartree`, record-relative `record_path` |
| `05_fixed_reference_placement_scores.csv` | treatment, family, placement | All 48 fixed-reference scores; treated cluster-minus-base score, fixed harmonic reference/water terms, selected flag and within-pair rank |
| `06_selected_candidates.csv` | treatment, family | Lower-score placement; `score_hartree` = G(cluster) − G(base); `fixed_reference_score_kcal_mol`; diagnostic `relative_g_kcal_mol` vs that treatment's selected adenine; `rank` |
| `07_continuum_rankings.csv` | treatment, family | `score_hartree` = G(continuum nucleotide) − G(base), referenced and ranked the same way |
| `08_selected_fixed_reference_shifts.csv` | treatment, family | Treatment-specific and harmonic selections, fixed-reference scores, their difference and whether the placement switched |
| `09_placement_combinations.csv` | treatment, choice_bits | `choice_bits`: six digits for adenine, cytosine, imidazole, melamine, TAP, guanine (0 = splitA, 1 = splitB); chosen candidate per family; `rank_order` (ties keep that family order); `c_minus_a_kcal_mol`; `pairwise_retained` (of 15 pairs, same sign as continuum); `ac_top_two` |
| `10_rank_occupancy.csv` | treatment, family, rank | Count among 64 enumerated placement choices; selected-lower-G rank and indicator |
| `thermochemical_summary.csv` | treatment | `selected_order`; `c_minus_a` (selected cytosine, kcal/mol); `pairwise_retained` (continuum vs selected); `unique_placement_orders` and `ac_top_two_choices` over the 64 vectors |

A catalogue of every published table is in [data/published/README.md](../data/published/README.md).
