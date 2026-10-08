"""Re-render the decomposition figure exactly from a reproduction's tables.

The net reaction split into nucleoside formation and phosphorylation, 21
families each, ranked relative to adenine.

Usage:
    python figure_scripts/decomposition.py --tables outputs/my_review/tables --out my_figures/
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import argparse
import yaml

from atcgu.plotting.decomposition import render_decomposition


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tables", type=Path, required=True, help="tables/ from an acgut reproduce run")
    p.add_argument("--out", type=Path, required=True, help="output directory for the figure")
    p.add_argument("--dpi", type=int, help="PNG resolution; default figures.yaml png_dpi")
    args = p.parse_args()
    config = yaml.safe_load((ROOT / "config/figures.yaml").read_text())
    dpi = args.dpi or int(config["figure1"]["png_dpi"])
    args.out.mkdir(parents=True, exist_ok=True)
    print("wrote", render_decomposition(args.tables.resolve(), args.out.resolve(), dpi, config))


if __name__ == "__main__":
    main()
