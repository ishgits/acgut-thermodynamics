# Microsolvation reference update

## Before and after

The previous SI panel A formed an adenine-relative value separately for every one of the 64 placement vectors. That made adenine identically zero, mixed variation in the reference placement with variation in the plotted family, and allowed the independently selected family point to fall away from a displayed range endpoint.

The revised panel evaluates every splitA/splitB endpoint against one fixed harmonic Figure 1C reference: the continuum-adenine nucleotide-minus-base score, two isolated harmonic waters and the −8.5396331666 kcal mol⁻¹ correction. Only the cluster and its corresponding base change entropy treatment. The 64-vector analysis remains separate and its rankings are unchanged.

## Harmonic cross-check against Figure 1C

| Family | splitA | splitB | Lower-G placement |
|---|---:|---:|---|
| Cytosine | 2.405896 | 0.226373 | splitB |
| Adenine | 4.485969 | 0.668710 | splitB |
| Imidazole | 4.386177 | 3.485375 | splitB |
| Melamine | 5.667398 | 3.950972 | splitB |
| TAP | 4.271274 | 4.701451 | splitA |
| Guanine | 7.568689 | 5.265874 | splitB |

Values are kcal mol⁻¹ on the fixed-reference scale. The 12 endpoints match the Figure 1C candidate coordinates within 10⁻⁸ kcal mol⁻¹.

## Selection and interpretation checks

- Harmonic and qRRHO select splitA only for TAP; every other family selects splitB.
- Both entropy-floor treatments select splitB for all six families, so TAP switches relative to harmonic.
- All four selected rankings are Cytosine < Adenine < Imidazole < Melamine < TAP < Guanine.
- Each of the 24 filled markers is the lower endpoint of its pair.
- The 256 treatment/choice records and their candidate IDs, bit strings, orders and diagnostics are retained. The common-zero change does not alter any rank order.
- Rank-occupancy rows sum to 64 for every family and treatment. Adenine and cytosine occupy the lowest two ranks in 32, 32, 36 and 32 combinations for harmonic, 50 cm⁻¹ floor, 100 cm⁻¹ floor and qRRHO, respectively.

The stable selected ordering should therefore be read separately from the broader sensitivity across the two sampled placements per family. Neither the connecting pairs nor the 64 enumerated choices constitute a hydration population or convergence estimate.
