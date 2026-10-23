"""Public commands; paths are resolved from a checkout, never the author's workspace."""
from pathlib import Path
import argparse
import csv
import json

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
    r = sub.add_parser("plot", help="Render Figure 1, the SI entropy figure and Figure S1 from a reproduction's tables")
    r.add_argument("--tables", type=Path, required=True); r.add_argument("--out", type=Path, required=True)
    r.add_argument("--dpi", type=int, help="PNG resolution; default figures.yaml png_dpi")
    v = sub.add_parser("validate", help="Check calculation records, source structures, sampling checksums and candidate selection")
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
            with (out / "reproduction_checks.csv").open(newline="") as h:
                rows = list(csv.DictReader(h))
            passed = sum(str(row["passed"]).lower() == "true" for row in rows)
            completion = (f"{passed}/{len(rows)} published-result comparisons passed.\n"
                          f"Report: {out / 'VALIDATION_REPORT.md'}")
        elif args.command == "plot":
            import yaml
            from .plotting.figure1 import render_figure1
            from .plotting.si_microsolvation import render_si_microsolvation
            from .plotting.decomposition import render_decomposition
            fresh(out)
            config = yaml.safe_load((root / "config/figures.yaml").read_text())
            dpi = args.dpi or int(config["figure1"]["png_dpi"])
            tables = args.tables.resolve()
            render_figure1(tables, out, dpi, config)
            # Only a missing SI input is a supported partial run. If the SI
            # tables are present, any validation or rendering error from the
            # SI renderer propagates as a failure instead of a skip.
            si_input = tables / "si" / "05_fixed_reference_placement_scores.csv"
            si_missing = not si_input.is_file()
            if not si_missing:
                render_si_microsolvation(tables, out, dpi, config)
            render_decomposition(tables, out, dpi, config)
            if si_missing:
                missing = "tables/si/" if not (tables / "si").is_dir() else "tables/si/05_fixed_reference_placement_scores.csv"
                completion = ("Partial completion: generated Figure 1 and Figure S1.\n"
                              f"Skipped SI entropy figure: missing {missing}.\n"
                              f"Output: {out}")
            else:
                completion = (f"Generated Figure 1, SI entropy figure and Figure S1.\n"
                              f"Output: {out}")
        elif args.command == "validate":
            from .validation.inventory import validate
            validate(root, out)
            summary = json.loads((out / "summary.json").read_text())
            total = int(summary["checks"]); passed = total - int(summary["failed_checks"])
            completion = (f"{passed}/{total} checks passed.\n"
                          f"Summary: {out / 'summary.json'}\n"
                          f"Detailed checks: {out / 'record_checks.csv'}")
        elif args.command == "audit":
            from .validation.arithmetic import audit
            audit(root, out)
            summary = json.loads((out / "summary.json").read_text())
            total = int(summary["checks"]); passed = total - int(summary["failed_checks"])
            completion = (f"{passed}/{total} checks passed.\n"
                          f"Summary: {out / 'summary.json'}\n"
                          f"Detailed checks: {out / 'independent_arithmetic_checks.csv'}")
    except (ValueError, RuntimeError, FileExistsError, FileNotFoundError, KeyError) as exc:
        p.exit(1, f"Study workflow failed: {exc}\n")
    print(completion)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
