"""Public commands; paths are resolved from a checkout, never the author's workspace."""
from pathlib import Path
import argparse

from .paths import check_output, fresh


def repository_root(value: str | None) -> Path:
    candidates = [Path(value).resolve()] if value else [Path.cwd(), *Path.cwd().parents, Path(__file__).resolve().parents[2]]
    for root in candidates:
        if (root / "config/molecules.csv").is_file() and (root / "data/calculations/manifest.csv").is_file():
            return root
    raise ValueError("Run from the repository checkout or specify --root /path/to/checkout")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Reproduce and verify the neutral nucleotide thermochemistry study")
    p.add_argument("--root", help="Checkout containing config/ and data/")
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("reproduce", help="Recompute the continuum analysis, Figure 1 and the SI sensitivity, and compare with the published results")
    r.add_argument("--out", required=True, type=Path); r.add_argument("--no-figures", action="store_true")
    r.add_argument("--dpi", type=int, help="PNG resolution; default figures.yaml png_dpi")
    r = sub.add_parser("plot", help="Render Figure 1, both SI figures and the decomposition figure from a reproduction's tables")
    r.add_argument("--tables", type=Path, required=True); r.add_argument("--out", type=Path, required=True)
    r.add_argument("--dpi", type=int, help="PNG resolution; default figures.yaml png_dpi")
    v = sub.add_parser("validate", help="Check calculation records, source structures and sampling checksums")
    v.add_argument("--out", required=True, type=Path)
    a = sub.add_parser("audit", help="Independent arithmetic checks using standard-library Python")
    a.add_argument("--out", required=True, type=Path)
    args = p.parse_args(argv)
    try:
        root = repository_root(args.root)
        out = check_output(root, args.out)
        if args.command == "reproduce":
            from .runner import run
            run(root, out, figures=not args.no_figures, dpi=args.dpi)
        elif args.command == "plot":
            import yaml
            from .plotting.figure1 import render_figure1
            from .plotting.si_microsolvation import render_si_microsolvation, render_si_ranking_occupancy
            from .plotting.decomposition import render_decomposition
            fresh(out)
            config = yaml.safe_load((root / "config/figures.yaml").read_text())
            dpi = args.dpi or int(config["figure1"]["png_dpi"])
            render_figure1(args.tables.resolve(), out, dpi, config)
            if (args.tables / "si").is_dir():
                render_si_microsolvation(args.tables.resolve(), out, dpi, config)
                render_si_ranking_occupancy(args.tables.resolve(), out, dpi, config)
            render_decomposition(args.tables.resolve(), out, dpi, config)
        elif args.command == "validate":
            from .validation.inventory import validate
            validate(root, out)
        elif args.command == "audit":
            from .validation.arithmetic import audit
            audit(root, out)
    except (ValueError, RuntimeError, FileExistsError, FileNotFoundError, KeyError) as exc:
        p.exit(1, f"Study workflow failed: {exc}\n")
    print(f"Completed {args.command}: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
