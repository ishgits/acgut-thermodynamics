# Nucleotide-formation thermochemistry: an open, verifiable study

**Purpose.** This repository is the complete, reproducible companion to our study of nucleotide-formation free energies across 21 canonical and noncanonical base families. It exists so that anyone — an undergrad learning computational chemistry, a reviewer checking a number, or a future collaborator — can see exactly where every reported value comes from and recompute it in minutes. If you want the science, read the manuscript first; this repo shows the work behind it.

The study compares 21 base families at B3LYP/PCM, repeats ten of them with three other functional/solvent combinations, and tests six nucleotide products with two explicit waters (manuscript Figure 1); the SI varies the entropy treatment and the water placements.

## How the repo is organized

```
                        ┌─────────────────────────────┐
                        │        THE DERIVATION         │
                        │  notebooks/  ← start here    │
                        │  (the analysis, step by step) │
                        └──────┬──────────────┬─────────┘
                               │              │
              ┌────────────────┘              └────────────────┐
              ▼                                                ▼
┌──────────────────────────┐                    ┌──────────────────────────┐
│       THE EVIDENCE         │                    │        THE CLAIMS        │
│  data/calculations/        │                    │  data/published/         │
│  420 calculation records:  │─── recompute ───▶  │  the manuscript's        │
│  geometry, frequencies,    │                    │  tables and plotted      │
│  energies, method, state   │                    │  data, as published      │
└──────────────────────────┘                    └──────────────────────────┘
┌──────────────────────────┐
│    HOW IT WAS COMPUTED     │
│  data/source/ — 66 starting structures (PubChem or constructed; see manifest) │
│  calculation_setup/ — the exact Gaussian/CREST inputs that were run (static)  │
│  src/atcgu/ — the analysis and plotting code the notebook calls               │
└──────────────────────────┘
```

Three ideas, that's the whole repo: **evidence** (the records), **derivation** (the notebook + code that turns records into results), **claims** (the published tables the derivation must reproduce). `atcgu reproduce` checks the third against the second, automatically.

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

Next, open `notebooks/01_reproduce_and_review.ipynb`. It walks the analysis the way the study was actually done:

1. **The data** — what the 420 records contain and what each field means
2. **Baseline reaction energies** — the net nucleotide-formation reaction across all 21 families
3. **Method/solvation sensitivity** — the same reaction for ten families under four computational treatments
4. **Water placement** — the two-water microsolvation test for the six leading families
5. **Decomposition** — splitting the net reaction into nucleoside formation and phosphorylation

The derivations (reaction energies, rankings, the numbers in the Results) are written out in the notebook; plumbing (record parsing, plotting) lives in `src/atcgu/` and is called, not hidden.

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
| `notebooks/` | The analysis, step by step — the place to learn and to work |
| `data/calculations/` | 420 calculation records (408 continuum, 12 two-water clusters): `record.json` + `optimized.xyz` + `frequencies.csv` + `starting.xyz` per record |
| `data/source/` | 66 starting structures and a manifest saying, truthfully, where each came from (`pubchem` → see manuscript Table 1; `constructed` → see the construction script) |
| `data/published/` | The published tables: `continuum/`, `figure1/`, `si_microsolvation/` |
| `data/sampling/` | CREST sampling inputs and ensembles behind the conformer selection |
| `calculation_setup/` | The exact inputs that were run (Gaussian `.com`, CREST jobs, water-placement starts) — a static record, not a generator |
| `config/` | Species, methods, reaction families, selection rules, figure conventions |
| `src/atcgu/` | The analysis, validation, and plotting code |
| `src/atcgu/building/molecules.py` | How the 20 programmatically constructed starting structures were built (with `config/construction.csv`) |
| `src/atcgu/sampling/water.py` | How the two-water microsolvation starting structures were constructed |
| `results/` | Rendered Figure 1, SI figures, captions, and plotted data |
| `docs/` | [Workflow](docs/WORKFLOW.md), [data dictionary](docs/DATA_DICTIONARY.md), [calculation setup](docs/CALCULATIONS.md), [troubleshooting](docs/TROUBLESHOOTING.md) |

## Notes

- Structure acquisition and Gaussian input generation were adapted from Ishaan Madan's [PubChem Gaussian Pipeline, DOI 10.5281/zenodo.21365151](https://zenodo.org/records/21365151); see [NOTICE.md](NOTICE.md). The shipped study inputs and the adaptations here define this study's workflow — don't assume upstream defaults match.
- The manuscript's Methods section is the scientific reference for *what* was calculated; this repo shows *how the numbers were derived from the records*. Nothing here repeats the paper.
- `.github/workflows/reproduce.yml` reruns `atcgu reproduce` on every push, so the verification stays live.
