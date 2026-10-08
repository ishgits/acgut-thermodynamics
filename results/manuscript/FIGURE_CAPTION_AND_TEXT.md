# Figure 1 caption and supporting text

## Figure caption

**Figure 1. Relative nucleotide-formation Gibbs energies and sensitivity to computational treatment.** Reaction families are named for their free bases; the products are the corresponding nucleotides. All values derive from harmonic thermochemistry at 298.0 K and are in kcal mol⁻¹. (A) Relative Gibbs energies for the balanced net reaction base + ribose + H₃PO₄ → nucleotide + 2 H₂O across 21 families, calculated at B3LYP/6-311++G(2df,2p) with aqueous IEFPCM using the lowest valid sampled Gibbs-energy structure for each species. Each reaction is referenced to the adenine-family reaction; lower relative values indicate more favorable formation within the modeled states. Green symbols and labels identify canonical families, and shaded rows identify the ten families examined in B. (B) Sensitivity to the four functional/continuum-solvent combinations shown in the legend. At the three alternative methods, the selected B3LYP/IEFPCM nucleotide structures were reoptimized and evaluated with method-matched reactants and adenine references; each method has its own adenine-family zero. (C) Product-side local-water sensitivity for the six lowest-scoring families in A. Two starting arrangements per nucleotide each place one explicit water at a phosphate OH and one at a base acceptor; both were optimized without constraints at B3LYP/IEFPCM. Open diamonds show continuum-only scores, filled circles the lower-G two-water placement, and open circles the other placement. Connecting lines show the difference between sampled placements, not statistical uncertainty. All C scores retain the harmonic continuum-only adenine nucleotide-minus-base score as the common zero. Two-water scores use the harmonic isolated-water reservoir and the standard-state/water-reservoir correction (−8.540 kcal mol⁻¹ at 298.0 K); reactants retain their continuum treatment. AICA, 5-aminoimidazole-4-carboxamide; TAP (C5-linked), C5-linked 2,4,6-triaminopyrimidine; Barbituric-C, C-linked barbituric acid. Other C- and N-linked labels specify the modeled linkage.

## Primary harmonic result

Selecting the lower of the two primary hydrated placements gives cytosine, adenine, imidazole, melamine, TAP and guanine. Under the common-zero water convention at 298 K, the selected adenine displacement is +0.669 kcal/mol and cytosine moves by −0.129 kcal/mol from its continuum score. Placement differences span 0.430–3.817 kcal/mol.

## Placement details

Each start places one water at a phosphate OH and one at a base acceptor. splitA/splitB identify starting arrangements; their targets are listed below. Construction and the unconstrained optimization protocol are summarized in the [calculation setup](../../docs/CALCULATIONS.md). A manuscript link will be added when it is public.

| Family | splitA base target | splitB base target | Phosphate side |
|---|---|---|---|
| Adenine | Ring N, input atom 10 | Different ring N, input atom 12 | Different OH groups |
| Cytosine | Ring N, input atom 11 | Carbonyl O, input atom 6 | Different OH groups |
| Guanine | Carbonyl O, input atom 8 | Ring N, input atom 12 | Different OH groups |
| TAP | Ring N, input atom 1 | Different ring N, input atom 2 | Different OH groups |
| Melamine | Ring N, input atom 1 | Different ring N, input atom 2 | Different OH groups |
| Imidazole | Ring N, input atom 2 | Same ring N; different water position/orientation | Different OH groups |

Atom numbers are one-based input indices. Placement manifests use zero-based indices.
