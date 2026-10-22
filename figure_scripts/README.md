# Figure scripts

Each script here re-renders one of the study's figures **exactly** from a
reproduction's tables — the same tested renderer the `acgut` commands use,
behind a small runnable file you can copy and adapt for your own figures.

## Usage

First produce the tables once:

```
acgut reproduce --out outputs/my_review
```

Then render any figure from those tables:

```
python figure_scripts/figure1.py --tables outputs/my_review/tables --out my_figures/
python figure_scripts/si_entropy.py --tables outputs/my_review/tables --out my_figures/
python figure_scripts/decomposition.py --tables outputs/my_review/tables --out my_figures/
```

`--dpi` overrides the PNG resolution (default: `config/figures.yaml`).

## What's what

| Script | Figure | Tested renderer |
|---|---|---|
| `figure1.py` | Figure 1 (panels A–C) | `src/atcgu/plotting/figure1.py` |
| `si_entropy.py` | SI entropy-treatment figure (4-panel) | `src/atcgu/plotting/si_microsolvation.py` |
| `decomposition.py` | Figure S1: three reaction steps in common Figure 1A family order | `src/atcgu/plotting/decomposition.py` |

Each script writes PNG, PDF and SVG with undated, reproducible metadata.
Tweak the script, keep the renderer import — or vendor the renderer into your
own project once you're happy with it.

Figure S1 retains `decomposition.png/pdf/svg` and writes `figure_s1_plotted_data.csv` (63 points). Figure 2 is a conceptual/manual asset; its final editable source has not been verified. These wrappers may overwrite same-named artwork in an existing directory, unlike the main CLI's fresh-directory rule.

The S1 benchmark was created once from the recovered final standalone Figure S1 script's guarded loader and the released ranking table, independently of the adapted renderer. It retains three reactions, the net-reaction family order, canonical flags and 63 values; it is not regenerated during tests.
