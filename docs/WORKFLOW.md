# Workflow

Start in the repository folder. Reproducing the included results needs Python 3.12 and the bundled data; it does not require running new chemistry calculations.

## Reproduce the results

1. **Set up Python** (once):

   ```bash
   python3.12 -m venv .venv
   source .venv/bin/activate
   python -m pip install -r environments/reproduction.txt
   python -m pip install --no-deps -e .
   ```

2. **Run the analysis:**

   ```bash
   atcgu reproduce --out outputs/first_run
   ```

3. **Check the run.** Open `outputs/first_run/VALIDATION_REPORT.md`. Its first result should say **28 of 28 published-result comparisons passed**. If the command stops or a check fails, see [Troubleshooting](TROUBLESHOOTING.md). Use a new name under `outputs/` for each run; existing runs are not overwritten.

4. **Explore the results.** Open `outputs/first_run/figures/` for Figure 1, the SI fixed-reference placement figure and the separate rank-occupancy figure. In `outputs/first_run/tables/`, start with `09_relative_rankings.csv` for the main comparison, `micro_06_primary_selected_candidates.csv` for the water-placement choices, `si/05_fixed_reference_placement_scores.csv` for all displayed endpoints and `si/thermochemical_summary.csv` for the sensitivity summary. The report links to the other useful tables.

5. **Follow one result back to its input.** For a guided walkthrough, install Jupyter with `python -m pip install '.[notebooks]'` if needed, then run `jupyter lab notebooks/01_reproduce_and_review.ipynb` and work through the cells.

## Go further

- For extra checks, run `atcgu audit --out outputs/audit_01` and `atcgu validate --out outputs/records_01`. Each command needs a new output folder.
- To redraw the figures from an existing run's tables, run `atcgu plot --tables outputs/first_run/tables --out outputs/plots_01`.
- To add new calculations, follow [Adding a molecule](EXTENDING.md): prepare and run the external jobs, extract and check their results, add the checked records to the repository, then run `atcgu analyze --out outputs/new_analysis`. [Calculation setup](CALCULATIONS.md) gives the job-preparation steps.
