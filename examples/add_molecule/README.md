# Worked neutral-linkage preparation example

The example reconstructs the study's **C-linked barbituric nucleotide**, with offline preparation steps and the matching study calculation for comparison.

The registry selects explicit-H barbituric-base atom 5, ribose-phosphate donor carbon 12 and donor OH oxygen 4, all **zero-based SDF indices**. Condensation removes the base H plus donor OH/OH-H and adds a C–C linkage, giving neutral C9H13N2O10P. N-linked barbituric products are chemically distinct; the same molecular formula does not justify merging them.

From the repository root:

```bash
atcgu build --registry examples/add_molecule/construction.csv --out outputs/example_build
atcgu prepare-input outputs/example_build/barbituric_acid_clink_ribose_phosphate.xyz \
  --method 1_b3lyp_pcm --name barbituric_clink_example --out outputs/example_input
python examples/add_molecule/inspect_existing_result.py
```

The final command calculates the example's balanced net reaction and adenine-relative score from the published selected records. For a new species, run the Gaussian and CREST/refinement jobs, then extract and validate the outputs. [Extending](../../docs/EXTENDING.md) explains registration.
