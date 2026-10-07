# Adding a molecule

The analysis consumes defined species, calculation records and explicit candidate cohorts. A new name or a PubChem hit alone cannot determine the intended protonation, tautomer, linkage or stereochemistry.

The supplied scientific parser, graph checks and coordinate/input helpers support CHNOP. Additional elements require explicit element/radius support and appropriate tests before using the same identification/validation rules.

## Runnable preparation example

`examples/add_molecule/construction.csv` reconstructs the C-linked barbituric nucleotide from its audited base/donor atom indices. It uses an existing study species so its formula and graph can be compared against bundled evidence. Run:

```bash
python -m pip install -r environments/structures.txt
atcgu build --registry examples/add_molecule/construction.csv --out outputs/linkage_example
atcgu prepare-input outputs/linkage_example/barbituric_acid_clink_ribose_phosphate.xyz \
  --method 1_b3lyp_pcm --name linkage_example --out outputs/linkage_example_gaussian
```

These commands work offline and create a construction manifest plus a Gaussian input. The example's [README](../examples/add_molecule/README.md) explains the attachment/removal choices and downstream reaction calculation.

## New chemical identity and source

Create a separate checkout for an extension and keep this release's published results unchanged. Choose a new stable species ID; define formula, charge, multiplicity, tautomer, protonation, linkage and intended stereochemistry. For PubChem acquisition put explicit CID choices in a CSV with `species_id,cid` and use `atcgu acquire --registry your_cids.csv --out outputs/acquired`. Network responses can change, so save response bytes and hashes as the new dataset's evidence.

For a constructed product, copy the example registry and choose the zero-based attachment atom in the explicit-H base SDF, donor carbon and donor OH oxygen. The neutral builder removes the attached base H plus donor OH/OH-H, adds the linkage and embeds/relaxes the product. Review formula, formal charge, fragment count, bond order and stereochemistry. Do not reuse atom indices after SDF reordering or implicit-H conversion. Charged/cationic alkylation requires a separately defined builder and analysis policy; it is not silently treated as neutral water-loss condensation.

## Calculate and sample

Prepare and run an initial opt/frequency job using `prepare-input`. Extract and validate it. Use its optimized geometry as the CREST seed, then select and prepare refinements with `atcgu prepare-refinement --moiety ...`, which applies that moiety's rule from `sampling` in `config/study.json`. Keep the initial DFT seed in the candidate set. Validate every output, including unsuccessful candidates. Record any additional tautomers or changes to the sampling design.

For alternate-method comparisons, begin with the selected primary-method geometry and use `prepare-input --method ...`. Calculate appropriate method-matched references. A transferred reoptimization is a distinct sampling design from a full alternate-method search.

## Register reviewed records and analyze

Register reviewed records in this order, in the extension checkout:

1. **Review.** Inspect `extraction_inventory.csv` and `independent_extraction_checks.csv` from `atcgu extract`, then each extracted `record.json`, geometry and frequency list: identity, absolute stereochemistry, charge/multiplicity, method and job boundaries.
2. **Copy bundles.** Copy each reviewed `data/calculations/<calculation_id>/` directory from the extraction output into the checkout's `data/calculations/`.
3. **Append manifest rows.** Append the matching rows of the extraction's `manifest_rows.csv` to `data/calculations/manifest.csv`. It already has the release manifest's nine columns in order (calculation ID, record path and SHA, source-log SHA and name, status, mode count, atom count, geometry file).
4. **Add configuration rows.** Add the species, family and candidate rows described below.
5. **Analyze.** Run `atcgu analyze --out outputs/extended_continuum`.

A valid manifest row means only that the bundle is complete and checksummed; it does not mean the chemical identity has been approved.

Add the species to `config/molecules.csv`. Its distance-inferred graph fingerprint is in the extracted record; review it against the intended structure first. Add a family to `config/reaction_families.csv` with its base/nucleoside/nucleotide IDs. Shared ribose, phosphate, ribose phosphate and water references must be present wherever required.

Append each candidate to `config/calculation_manifest.csv` using the existing schema. Fields required for selection include method ID, expected functional/solvent model, species and class, role, inclusion, original source filename, full SHA, public `local_path`, calculation ID and record hash. Role `conformer_candidate` is for the primary sampled pool. Alternate-method fixed products use `required_transferred_nucleotide` or `transferred_nucleoside`; references use `method_matched_reference`. Diagnostics must remain excluded from primary selection. Assign a new analysis ID and update declared family counts/labels in `analysis_sets.csv`.

The `all_21` family-panel key names the full panel and reads all rows in the family registry. You can extend that full continuum panel and change its analysis ID/declared size; if you rename the full panel's analysis ID, update `panel_a_analysis_id` in `config/figures.yaml` too. Targeted and canonical lists come from `family_panels` in `study.json`; edit them explicitly when defining different subsets. `load_inventory` returns these panels, checks each analysis set's declared family count against them, and passes them to every selection step. The primary water panel and Figure 1 remain fixed to this released study. A larger study needs its own figure/panel definitions and layout.

Run `atcgu analyze --out outputs/extended_continuum` to compute the new continuum tables without demanding equality to original results. Recheck coverage, method/state, all balanced reactions, path closure and reference families before interpreting the extension. Add its explicit stoichiometry to `reaction_stoichiometry.csv` for an independent auditor; the shipped `audit` intentionally compares only with this study's published results. Version the new dataset, record checksums for its results and add independent checks for its scientific claims.
