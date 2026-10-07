# Computational definitions and scope

The modeled chemical states are neutral singlets. A family consists of a free base and the explicitly defined nucleoside/nucleotide linkage products. C-linked and N-linked barbituric families share the free-base species but have distinct product identities. These comparisons concern the modeled states and sampled structures, not all tautomers, protonation states or stereoisomers in bulk water.

## Continuum study

All four methods use 6-311++G(2df,2p) with aqueous IEFPCM or SMD: B3LYP/PCM, B3LYP/SMD, wB97X-D/PCM and wB97X-D/SMD. Harmonic thermochemistry is printed at 298.0 K and 1 atm. The analysis computes G = final frequency-job electronic energy + thermal Gibbs correction in hartree per particle. Printed total G provides a rounding check. The conversion is 627.5094740631 kcal mol^-1 per hartree.

The full 21-family analysis selects the lowest usable harmonic G for every required species from the full B3LYP/PCM pool. Ties use the original source filename. The ten-family method comparison selects the fixed transferred nucleotide at each alternate method, together with the specified method-matched reactants/references. This is reoptimization of selected B3LYP/PCM structures, not independent conformer sampling at every method. For 45 transferred products a second, diagnostic-only optimization also exists; `10_legacy_vs_transferred_audit.csv` compares the two, and the diagnostic structure is never selected.

The two reaction paths are:

1. Base + ribose → nucleoside + water; nucleoside + H3PO4 → nucleotide + water.
2. Ribose + H3PO4 → ribose phosphate + water; base + ribose phosphate → nucleotide + water.

Both sum to base + ribose + H3PO4 → nucleotide + 2 water. Elements, charge, required species and path closure are checked. `config/reaction_stoichiometry.csv` provides signed coefficients independently of the production reaction builder. All relative values subtract the reference-family reaction (adenine, `reference_family` in `config/study.json`) within the same method and reaction scope. Ranks within a panel use `rank(method="min")`, so tied energies would share a rank; the SI tables sort deterministically and rank with `method="first"`. The data contain no ties. Lower scores indicate more favorable modeled formation. They are not kinetic predictions.

## Structure sampling

CREST uses GFN2-xTB, ALPB water, a 6 kcal/mol retained window and its default iMTD-GC workflow. The seed coordinates and the CREST ensembles used are included. Energies are relative to the first frame, which must be the ensemble minimum.

- Nucleotides and ribose phosphate: frames within 6 kcal/mol are divided into three heavy-atom radius-of-gyration bins, and the lowest two per bin are kept. The upper bin spans at least the seed's radius.
- Bases, nucleosides and the references ribose and phosphoric acid (moiety `reference`): frames are taken in file order until the first one above 3 kcal/mol, keeping up to five.
- A frame is skipped when its heavy-atom Kabsch RMSD to any already kept frame is below 0.5 Å. Every selection also carries the DFT seed.

The rules, windows and RMSD cutoff are read from `sampling` in `config/study.json`; `config/sampling.csv` assigns each run its moiety. The selector preserves CREST file order and conformer IDs.

## Primary Figure 1C

Six nucleotide products each have two primary split-site starts. One water accepts from phosphate OH; a second donates to a base acceptor. A/B use distinct phosphate OH groups. They also target distinct base acceptors except imidazole, which has the same ring nitrogen approached differently. Each start must have exactly one phosphate→water and one water→base contact (counted as donor-H/acceptor pairs), no bridging water, and an O···O distance of at least 2.50 Å. Open families start from an unfolded dry nucleotide and reject folded starts; adenine and cytosine (`folded_start` in `config/microsolvation_references.csv`) start from their direct-folded dry minimum, keep the fold, and move only the rotated phosphate H. All optimizations are unconstrained. Two intact water fragments and the specified contacts/folding classification are checked from the final XYZ, where contacts are counted as distinct waters.

Within-condition hydrated scores use G(cluster) − G(base), referenced to the selected adenine cluster. The figure instead keeps continuum adenine as a common zero. For family i it displays:

`x_dry(i) = ([G(nucleotide_i) − G(base_i)] − [G(nucleotide_A) − G(base_A)]) × conversion`

`x_hydrated(i) = x_dry(i) + [G(cluster_i) − G(nucleotide_i) − 2 G(water)] × 627.509474 + C`

`C = −2 R T ln(24.453091396947848) − 2 R T ln(55.34) = −8.539633166577238 kcal/mol`

The correction and Gaussian thermochemistry both use 298.0 K. R = 0.0019872041 kcal mol^-1 K^-1; the 1 atm to 1 M factor is evaluated at that temperature. Constants are defined in `config/figures.yaml` and `config/study.json`.

Reactants retain their continuum treatment. Placement spans show differences between available optimized starts, not statistical error bars. Two starts do not establish equilibrium hydration-shell or conformational convergence. Figure 1C uses harmonic thermochemistry; the treatments below are SI sensitivity checks.

## SI microsolvation sensitivity

`atcgu reproduce` applies four thermochemical treatments to the 12 primary clusters and the 12 continuum references (six bases, six continuum nucleotides), from the released frequencies at 298.0 K (`src/atcgu/analysis/thermochemistry_sensitivity.py`):

- `harmonic`: Gaussian's harmonic G.
- `floor_50`, `floor_100`: the Cramer/Truhlar treatment. Each real mode below the cutoff is raised to 50 or 100 cm^-1 for its vibrational entropy only; G changes by T × [S_harmonic − S_raised].
- `grimme_qrrho_100`: Grimme's entropy-only quasi-RRHO. Each mode's entropy interpolates between the harmonic oscillator and a free rotor with damping 1/[1 + (100 cm^-1 / ν)^4] and average moment of inertia 1 × 10^-44 kg m^2.

E, ZPE and H are unchanged and frequencies are unscaled. Every record must have all-real modes, a complete frequency file and T within 0.25 K of 298.0 K.

The displayed SI placement score keeps Figure 1C's harmonic reference fixed for every family and entropy treatment:

```text
X(i,p,t) = 627.509474 × [G_t(cluster_i,p) − G_t(base_i)
             − (G_harmonic(continuum_nucleotide_adenine) − G_harmonic(base_adenine))
             − 2 G_harmonic(water)] + C
```

Here `C = −8.5396331666 kcal mol⁻¹` is the same 298.0 K standard-state/water-reservoir correction used in Figure 1C. Thus the cluster and corresponding base receive treatment `t`, while the continuum-adenine zero, isolated-water reservoir and `C` remain harmonic and fixed. The fixed water term is deliberate: the SI changes only the cluster/base vibrational entropy and does not introduce treatment-dependent water energies or zeros. For each family and treatment the filled endpoint is the lower of splitA and splitB (ties by candidate ID). Connecting lines compare two sampled optimized placements; they are not uncertainty intervals or a converged hydration ensemble.

The entropy-shift panel reports `min_p X(i,p,t) − min_p X(i,p,harmonic)`. Each minimum is selected independently, so the quantity includes a placement switch where one occurs; it is not necessarily a same-structure correction. Treatment-specific adenine-relative scores are retained in the tables as ranking diagnostics, but they are not the displayed scale.

Separately, the 64-combination design fixes one placement per family in the order adenine, cytosine, imidazole, melamine, TAP, guanine (bit 0 = splitA, 1 = splitB). For each combination the tables retain the candidate identities, resulting order, cytosine − adenine gap, whether A and C are the lowest two, and how many of the 15 family pairs keep the **treatment-matched continuum** sign. Adding the fixed common zero cannot change these rankings. The rank-occupancy figure counts how often each family occupies ranks 1–6 among those 64 enumerated choices; these are counts, not probabilities, populations, confidence levels or equally likely physical states. Results are in `data/published/si_microsolvation/` and `results/supplement/`.

## Identity and validation limits

The molecular registry checks formula and distance-inferred connectivity for CHNOP structures. Bond orders, tautomer assignments and absolute stereochemistry require inspection of the intended structures. Starting SDFs retain supplied bond and stereochemical annotations.
