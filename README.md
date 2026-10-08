# Thermodynamics of nucleotide formation

**Purpose.** This repository is the reproducible companion to our study of nucleotide-formation free energies across 21 canonical and noncanonical base families. It lets readers trace reported values to the released calculation records and recompute them in minutes. The manuscript and preprint are not public yet; a link will be added here when available.

The study compares 21 base families at B3LYP/PCM, repeats ten of them with three other functional/solvent combinations, and tests six nucleotide products with two explicit waters (manuscript Figure 1); the SI varies the entropy treatment and the water placements.

## How the repo is organized

```
data/calculations/           src/atcgu/                 data/published/
420 calculation records ──▶ analysis and checks ──▶ published tables
         │                         │                         │
         └─────────────────────────┴─────────────────────────┘
                                   │
                       notebooks/ (guided review)

data/source/ and calculation_setup/ document starting structures and inputs.
```

Three ideas organize the repo: **evidence** (the records), **derivation** (`src/atcgu/` code that turns records into results), and **claims** (the published tables). `atcgu reproduce` recomputes the tables and checks them against the published values. The notebook is a guided review of selected steps and results.

## Start here (10 minutes)

You need Python 3.12. No Gaussian, CREST, or quantum chemistry software — the records already contain everything.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r environments/reproduction.txt
python -m pip install --no-deps -e .
atcgu reproduce --out outputs/my_review
```

Then open `outputs/my_review/VALIDATION_REPORT.md`: 28 comparisons between the recomputed tables and the published ones, all of which must pass. That report *is* the verification — there is deliberately no separate integrity-checking machinery in this repo; git tracks the files, and the reproduce report checks the numbers.

Next, open `notebooks/01_reproduce_and_review.ipynb`. To run it interactively, install Jupyter with `python -m pip install '.[notebooks]'`. It walks through the evidence and headline results:

1. **The data** — what the 420 records contain and what each field means
2. **Baseline reaction energies** — the net nucleotide-formation reaction across all 21 families
3. **Method/solvation sensitivity** — the same reaction for ten families under four computational treatments
4. **Water placement** — the two-water microsolvation test for the six leading families
5. **Decomposition** — splitting the net reaction into nucleoside formation and phosphorylation

The notebook computes baseline reaction energies from the published selected-species table and uses other published intermediate tables to explain later results. For the full derivation from calculation records, run `atcgu reproduce`.

## The commands

| Command | What it does |
|---|---|
| `atcgu reproduce --out <dir>` | Recompute every table, render Figure 1 + SI figures, compare against `data/published/` (28 checks) |
| `atcgu audit --out <dir>` | Independent re-derivation of selections, energies, ranks and plotted points using only the standard library |
| `atcgu validate --out <dir>` | Check the calculation records, source structures, and sampling checksums |
| `atcgu plot --tables <dir>/tables --out <dir>` | Re-render the figures from an existing reproduction |

Each run needs a fresh output directory. Add `--no-figures` to `reproduce` to skip rendering.

## What's where

| Location | Contents |
|---|---|
| `notebooks/` | Guided review of the records, selected calculations, and headline results |
| `data/calculations/` | 420 calculation records (408 continuum, 12 two-water clusters): `record.json` + `optimized.xyz` + `frequencies.csv` + `starting.xyz` per record |
| `data/source/` | 66 starting structures and a source manifest (`pubchem` CIDs are embedded in the SDF files; `constructed` recipes are in the construction script) |
| `data/published/` | The published tables: `continuum/`, `figure1/`, `si_microsolvation/` |
| `data/sampling/` | CREST sampling inputs and ensembles behind the conformer selection |
| `calculation_setup/` | CREST and water-placement inputs, plus Gaussian rerun inputs generated from optimized geometries |
| `config/` | Species, methods, reaction families, selection rules, figure conventions |
| `src/atcgu/` | The analysis, validation, and plotting code |
| `src/atcgu/building/molecules.py` | How the 20 programmatically constructed starting structures were built (with `config/construction.csv`) |
| `src/atcgu/sampling/water.py` | How the two-water microsolvation starting structures were constructed |
| `results/` | Rendered Figure 1, SI figures, captions, and plotted data |
| `docs/` | [Workflow](docs/WORKFLOW.md), [data dictionary](docs/DATA_DICTIONARY.md), [calculation setup](docs/CALCULATIONS.md), [troubleshooting](docs/TROUBLESHOOTING.md) |

## Notes

- Structure acquisition and Gaussian input generation were adapted from Ishaan Madan's [PubChem Gaussian Pipeline, DOI 10.5281/zenodo.21365151](https://zenodo.org/records/21365151); see [NOTICE.md](NOTICE.md). The shipped study inputs and the adaptations here define this study's workflow — don't assume upstream defaults match.
- The manuscript's Methods section will be the scientific reference for *what* was calculated. Its link will be added when the manuscript or preprint is public.
- `.github/workflows/reproduce.yml` reruns `atcgu reproduce` on every push, so the verification stays live.
