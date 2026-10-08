# SI microsolvation figures and supporting text

## Fixed-reference placement and entropy figure

**SI Figure. Two-water placement sensitivity across entropy treatments on the Figure 1C fixed-reference scale.** Scores are B3LYP/6-311++G(2df,2p), aqueous IEFPCM quantities in kcal mol⁻¹. In every treatment the displayed score subtracts the same harmonic continuum-adenine nucleotide-minus-base quantity and two harmonic isolated-water Gibbs energies, then applies the same −8.539633 kcal mol⁻¹ standard-state/water-reservoir correction used in Figure 1C. Only each cluster and its corresponding base receive the named entropy treatment. Each of the four panels joins the two sampled placements per family: circles denote splitA and squares splitB. Connecting lines compare two optimized sampled placements, not statistical uncertainty or a converged hydration ensemble. The per-family minima keep the selected order C < A < imidazole < melamine < TAP < G under every treatment. Adenine is not forced to zero in any panel.

## Placement-choice sensitivity (text point)

Enumerating all 64 splitA/splitB placement choices per treatment (`09_placement_combinations.csv`, counts in `10_rank_occupancy.csv`) shows the selected order is stable while the middle ranks are placement-sensitive: only 32–36 of the 64 vectors put adenine and cytosine in the lowest two ranks. This separates robustness of the selected ordering from uncertainty associated with the limited placement sample; no figure is drawn for it.

## Summary

| Treatment | Selected order | C − A (kcal mol⁻¹) | Pairs retained from treatment-matched continuum (of 15) | Distinct orders over 64 combinations | Combinations with A and C first two |
|---|---|---|---|---|---|
| Harmonic | C < A < imidazole < melamine < TAP < G | −0.442 | 13 | 24 | 32 |
| 50 cm⁻¹ floor | C < A < imidazole < melamine < TAP < G | −0.071 | 14 | 18 | 32 |
| 100 cm⁻¹ floor | C < A < imidazole < melamine < TAP < G | −0.011 | 15 | 15 | 36 |
| Grimme qRRHO (100 cm⁻¹) | C < A < imidazole < melamine < TAP < G | −0.239 | 15 | 18 | 32 |

The selected order is stable under all four entropy treatments. The occupancy panels nevertheless show substantial sensitivity to alternative sampled placements: only 32–36 of the 64 enumerated vectors put adenine and cytosine in the lowest two ranks. This contrast separates robustness of the selected ordering from uncertainty associated with the limited placement sample.

## Sources

- Fixed convention: `data/published/figure1/03_fixed_reference_convention.csv`
- All 48 displayed placement scores: `data/published/si_microsolvation/05_fixed_reference_placement_scores.csv`
- Selected-score shifts: `08_selected_fixed_reference_shifts.csv`
- All 256 placement/treatment combinations: `09_placement_combinations.csv` (`choice_bits` lists adenine, cytosine, imidazole, melamine, TAP, guanine; 0 = splitA, 1 = splitB)
- Rank counts: `10_rank_occupancy.csv`
- Plotted data: `results/supplement/si_microsolvation_plotted_data.csv`

The unequal-sampling historical-basin analysis is not part of this release.
