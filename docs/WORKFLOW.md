# Workflow

This is the whole user workflow: reproduce, verify, then learn by following the derivation. Everything starts from a checkout of this repository.

## 1. Reproduce the results

Python 3.12, bundled data only — no chemistry software needed.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r environments/reproduction.txt
python -m pip install --no-deps -e .
atcgu reproduce --out outputs/first_run
```

## 2. Check the run

Open `outputs/first_run/VALIDATION_REPORT.md`. It should say **28 of 28 published-result comparisons passed** — every recomputed table matched against `data/published/`. That report is the verification; if a check fails, see [Troubleshooting](TROUBLESHOOTING.md). Use a fresh directory under `outputs/` for each run.

For an independent second opinion: `atcgu audit --out outputs/audit_01` re-derives the selections, energies, ranks, and plotted points using only the Python standard library (no project code). `atcgu validate --out outputs/records_01` checks the calculation records, source structures, and sampling checksums.

## 3. Follow the derivation

Install Jupyter (`python -m pip install '.[notebooks]'`) and open `notebooks/01_reproduce_and_review.ipynb`. It walks the analysis the way the study was actually done:

1. **The data** — what the 420 records contain and what each field means
2. **Baseline reaction energies** — the net nucleotide-formation reaction across all 21 families
3. **Method/solvation sensitivity** — ten families under four computational treatments
4. **Water placement** — the two-water test for the six leading families
5. **Decomposition** — the net reaction split into nucleoside formation and phosphorylation

The derivations are written out in the notebook cells; plumbing (record parsing, plotting) lives in `src/atcgu/` and is called from the notebook, not duplicated in it. To trace any single number — say, the 0.36 kcal/mol adenine–cytosine separation — start from `outputs/first_run/tables/09_relative_rankings.csv` and follow the notebook back to the records it came from. The [data dictionary](DATA_DICTIONARY.md) defines every table and field.

## 4. See how the inputs were made

`calculation_setup/` holds the exact inputs that were run (Gaussian `.com` files, CREST jobs, water-placement starts) — a static record. The two programmatic construction steps are scripted and readable: `src/atcgu/building/molecules.py` (+ `config/construction.csv`) built the 20 non-PubChem starting structures, and `src/atcgu/sampling/water.py` built the two-water starting arrangements. [Calculation setup](CALCULATIONS.md) describes the layout.

That's the entire workflow. There is deliberately no machinery here for extending the study to new molecules — this repo documents and verifies the published analysis; extending it is future work.
