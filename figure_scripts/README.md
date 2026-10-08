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
| `decomposition.py` | Nucleoside-formation vs phosphorylation rankings (SI candidate) | `src/atcgu/plotting/decomposition.py` |

Each script writes PNG, PDF and SVG with undated, reproducible metadata.
Tweak the script, keep the renderer import — or vendor the renderer into your
own project once you're happy with it.
