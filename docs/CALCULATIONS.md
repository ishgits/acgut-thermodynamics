# Calculation setup

This page documents the released calculation inputs. `calculation_setup/` is a static record; nothing here reruns the chemistry. Its Gaussian `.com` files are reconstructed rerun inputs made from the released optimized geometries, not copies of the original submitted inputs.

## What ran

- **Gaussian 16** for DFT optimizations and frequencies (revision recorded per calculation in its record). Rerun-style inputs use each record's `optimized.xyz`, plain `Opt` at the configured method, frequency analysis at 298 K.
- **CREST 3.0.2** (commit 6914a25) with **GFN2-xTB** through the tblite backend for conformer searching (`--gfn2 --alpb water --ewin 6.0`). Software settings are indexed in `provenance/sampling_software.csv`.

## What's in `calculation_setup/`

| Location | Contents |
|---|---|
| `gaussian/jobs/` | Reconstructed `.com` rerun inputs for 417 usable records; three failed sources were skipped and are documented in the records |
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
