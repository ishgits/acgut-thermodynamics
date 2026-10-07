# SI microsolvation figures and supporting text

## Fixed-reference placement and entropy figure

**SI Figure. Two-water placement and entropy-treatment sensitivity on the Figure 1C fixed-reference scale.** Scores are B3LYP/6-311++G(2df,2p), aqueous IEFPCM quantities in kcal mol⁻¹. In every treatment the displayed score subtracts the same harmonic continuum-adenine nucleotide-minus-base quantity and two harmonic isolated-water Gibbs energies, then applies the same −8.539633 kcal mol⁻¹ standard-state/water-reservoir correction used in Figure 1C. Only each cluster and its corresponding base receive the named entropy treatment. (A) Four connected splitA/splitB pairs per family. Filled markers identify the lower-G sampled placement selected independently for that treatment; open markers identify the other sampled placement. Circles denote splitA and squares splitB. Connecting lines compare two optimized sampled placements, not statistical uncertainty or a converged hydration ensemble. (B) `min_p X(i,p,t) − min_p X(i,p,harmonic)` for each non-harmonic treatment. The minima are reselected independently; outlined symbols identify a change from the harmonic placement, so these shifts do not always compare the same optimized structure. Adenine is not forced to zero in either panel.

## Rank-occupancy figure

**SI Figure. Rank occupancy over enumerated two-water placement choices.** Each treatment panel counts, for every family and rank, how many of the 64 splitA/splitB choice vectors place that family at ranks 1–6. Cell color and labels use one shared count scale across all four panels. Black outlines mark the rank produced when the lower-G sampled placement is selected independently for every family. These are counts over enumerated choices—not probabilities, populations, confidence levels or equally likely physical states. The ranking comparison retains the treatment-matched continuum comparator used for `pairwise_retained`; the fixed display zero does not change any combination's order.

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
- Plotted data: `results/supplement/si_microsolvation_plotted_data.csv` and `si_microsolvation_ranking_plotted_data.csv`

The unequal-sampling historical-basin analysis is not part of this release.
