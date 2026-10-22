# Calculation setup

This page documents the released calculation inputs. `calculation_setup/` is a static record; nothing here reruns the chemistry. Its Gaussian `.com` files are reconstructed rerun inputs made from the released optimized geometries, not copies of the original submitted inputs.

## What ran

- **Gaussian 16** for DFT optimizations and frequencies (revision recorded per calculation in its record). Rerun-style inputs use each record's `optimized.xyz`, plain `Opt` at the configured method, frequency analysis at 298 K.
- **CREST 3.0.2** (commit 6914a25) with **GFN2-xTB** through the tblite backend for conformer searching (`--gfn2 --alpb water --ewin 6.0`). Software settings are indexed in `provenance/sampling_software.csv`.

## What's in `calculation_setup/`

| Location | Contents |
|---|---|
| `gaussian/jobs/` | Reconstructed `.com` rerun inputs for all 417 records |
| `crest/jobs/` | CREST search inputs: seed geometries, sampling configs, run scripts |
| `crest/refinement/` | DFT refinement starts selected from the CREST ensembles, with `refinement_manifest.json` recording the selection policy, conformer numbers, and relative energies |
| `microsolvation/jobs/` | The two-water starting arrangements per family |
| `microsolvation/original_inputs/` | The twelve original water-start XYZ geometries, indexed by `original_input_manifest.csv` |
| `scheduler_templates/` | Example Slurm/batch wrappers (customize modules, account, queue for your system) |

Refinement and water-placement inputs used tight optimization routes (`opt=(tight,calcfc,maxcycles=200)` and `opt=(calcfc,tight,maxcycles=200,maxstep=10)`). On Apple silicon, regenerating the guanine splitA and melamine splitB water starts can differ from the shipped inputs in the last digits — the shipped inputs are authoritative.

## Programmatic construction

Two construction steps were scripted rather than hand-built:

- `src/atcgu/building/molecules.py` (+ `config/construction.csv`) built the 20 non-PubChem starting structures by explicit neutral condensation with audited atom choices and fixed random seeds.
- `src/atcgu/sampling/water.py` built the two-water starting arrangements: one water accepts from a phosphate OH, the other donates to a base acceptor, scanning phosphate OH torsions and water orientations, rejecting unintended contacts. A/B labels identify starts, not energetic order.

## Records

Every calculation's inputs and outputs are bundled in `data/calculations/<id>/`: `record.json` (energies, method, state, checksums), `optimized.xyz`, `frequencies.csv`, and `starting.xyz`. The notebook starts from these bundles — no Gaussian or CREST needed.

## Released candidate verification

```bash
acgut validate --out outputs/candidate_review
python - <<'PYCODE'
import pandas as pd
r = pd.read_csv("outputs/candidate_review/sampling_selection_checks.csv")
print(r[r.sampling_run_id.isin(["adenine", "adenine_ribose_phosphate"])].to_string(index=False))
PYCODE
```

The pure selector was recovered from Git commit `4d2957ac3023deea7c9bc320fc53f7acca583111`; the preparation adapter was deliberately not restored. `config/sampling.csv` supplies case-sensitive sampling run IDs (including `tap_C5link`), seed/ensemble paths and moieties. `config/study.json` supplies policy. The verifier writes one row per run, then aggregates failures without changing any manifest.

Heavy-atom RMSD excludes H, centers corresponding atoms and uses proper-rotation Kabsch alignment without symmetry permutations. Rg is unweighted. Duplicate means RMSD strictly below 0.5 Å; equality survives. Energy selection keeps file order, stops at the first frame above 3 kcal/mol, and keeps at most five distinct frames. Rg selection admits energies ≤6 kcal/mol, uses three equal-width bins spanning eligible minimum Rg to max(eligible maximum Rg, seed Rg), adds 1e-9 Å to the final edge, includes lower/excludes upper edges, and retains at most two per bin with deduplication across all previous bins. The recovered `np.argsort` energy tie behavior depends on NumPy; use `environments/reproduction.txt` (NumPy 2.3.5), rather than inventing a new tie rule. Final selected frames are energy-ordered with the recovered stable Python sort.

Released preparation candidate sets contain `dft_seed` plus selected CREST frames and **no seed-similarity filter**. All 66 manifests reproduce using 627.5094740631 kcal/mol per hartree. This establishes candidate IDs/order, metrics, policy, hashes and adjacent-energy inversions, not executed jobs or byte-identical Gaussian inputs. Actual calculation evidence is the calculation-record inventory.

Historical campaign dispatch used 627.509 and an additional step:

```text
survivors = selected CREST frames with RMSD(frame, seed) >= 0.5 Å
if survivors is empty: record "existing"; create no new jobs
else: campaign starts = [dft_seed] + survivors
```

That historical replay matched 11 pilot and 12 follow-up species. Applying it counterfactually to all 66 released ensembles changes 31 candidate lists; this does not establish incorrect executed calculations. Conversion alone changes no selected IDs but can change relative energies by 4.43e-6 kcal/mol. The 420 historical logs included three failures; 417 valid calculation records are released.
