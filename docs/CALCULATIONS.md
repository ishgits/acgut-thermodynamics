# Preparing and rerunning external calculations

Rerunning the calculations requires Gaussian 16 and CREST/GFN2-xTB for conformer searching. Gaussian revisions are recorded per calculation.

The study used CREST 3.0.2 (commit 6914a25) with GFN2-xTB through the tblite backend. Software settings are indexed in `provenance/sampling_software.csv`. The general preparation command defaults to 8 threads; change this with `--threads`.

## Ready-made study inputs

`calculation_setup/gaussian/jobs/input_manifest.json` indexes all 420 records: 417 usable jobs have `.com` inputs, and three failed sources are skipped. Rerun inputs use each record's `optimized.xyz`, plain `Opt` at its configured method, and frequency analysis at 298 K. CPU/memory defaults are 8 cores/6 GB, configurable with `--cpus` and `--memory`. Refinement and water-placement inputs (`prepare-refinement`, `prepare-waters`) reproduce the submitted tight-optimization routes (`opt=(tight,calcfc,maxcycles=200)` and `opt=(calcfc,tight,maxcycles=200,maxstep=10)`), whereas `prepare-gaussian` reruns use default `Opt`.

The twelve original water-start XYZ geometries are retained in `calculation_setup/microsolvation/original_inputs/`, indexed by `original_input_manifest.csv`.

## New or transferred Gaussian input

```bash
atcgu prepare-input data/calculations/calc_565406d4de58dfbc/optimized.xyz \
  --method 2_b3lyp_smd --name adenine_nucleotide_smd \
  --cpus 8 --memory 6GB --out outputs/transferred_input
```

This takes the optimized B3LYP/PCM nucleotide and prepares plain `Opt` followed by frequencies at 298 K at the requested configured method. Set charge/multiplicity for different chemical states; the study uses neutral singlets.

After loading Gaussian, run the prepared input:

```bash
bash calculation_setup/scheduler_templates/run_gaussian.sh outputs/transferred_input/adenine_nucleotide_smd.com
```

The wrapper refuses to overwrite an existing log. For scheduled jobs, customize the modules, account, queue, memory and wall time in `gaussian.slurm`.

## CREST and refinement

```bash
atcgu prepare-crest data/sampling/adenine_ribose_phosphate/seed.xyz \
  --threads 8 --out outputs/adenine_crest
# With CREST/xTB installed:
bash outputs/adenine_crest/run_crest.sh
```

The prepared command uses `--gfn2 --alpb water --ewin 6.0`, charge and unpaired-electron count.

To prepare refinements from the released ensemble:

```bash
atcgu prepare-refinement data/sampling/adenine_ribose_phosphate/seed.xyz \
  data/sampling/adenine_ribose_phosphate/crest_conformers.xyz \
  --moiety nucleotide --out outputs/adenine_refinement
```

The manifest records the selection policy, conformer numbers, relative energies and radius of gyration. Supply `--moiety ribo_phosphate` for that reference and `--moiety reference` for ribose/phosphoric acid; each moiety's rule is read from `sampling` in `config/study.json`. Ready-made jobs are in `calculation_setup/crest/`.

## Primary water placements

```bash
atcgu prepare-waters data/calculations/calc_0fc80346cb0cfb89/optimized.xyz \
  --family imidazole --out outputs/imidazole_water_starts
atcgu prepare-waters data/calculations/calc_565406d4de58dfbc/optimized.xyz \
  --family adenine --folded --out outputs/adenine_water_starts
```

The open construction rule was used for G/TAP/M/I; `--folded` invokes the A/C rule (`folded_start` in `config/microsolvation_references.csv`). Construction scans phosphate OH torsions and water directions/orientations, rejects unintended contacts, checks intended H-bond geometry and separates two starts. It preserves solute heavy atoms initially, while allowing the selected OH hydrogen to rotate. The optimizations then have no geometry constraints. Both waters are present simultaneously. A/B labels identify starts, not energetic order or final folding states.

The prepared input uses B3LYP/6-311++G(2df,2p) with IEFPCM water. Standard template resources are 12 cores/40 GB. Placement parameters are in each `placement_manifest.json`.

`workflows/prepare_study.py` regenerates `calculation_setup/gaussian/`, `crest/jobs/` and `crest/refinement/` byte for byte. On Apple silicon the guanine splitA and melamine splitB water starts differ from the shipped inputs in the last digits; the shipped inputs are authoritative.

## Extract and inspect new results

```bash
atcgu extract outputs/transferred_input/adenine_nucleotide_smd.log \
  --out outputs/transferred_records
```

The extractor selects the final frequency job and retains complete frequencies and job boundaries. Supported jobs require optimization evidence, normal frequency-job termination and CHNOP coordinates. Every record is compared field by field with an independent parser before anything is written; if any comparison disagrees, the command writes only `independent_extraction_checks.csv` and `extraction_inventory.csv` and exits with an error, so no bundles are produced. Inspect failure reasons and extraction comparisons before adding records to a study.
