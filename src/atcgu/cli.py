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
    p = argparse.ArgumentParser(description="Reproduce and extend the neutral nucleotide thermochemistry study")
    p.add_argument("--root", help="Checkout containing config/ and data/")
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("reproduce", help="Recompute the continuum analysis, Figure 1 and the SI sensitivity, and compare with the published results")
    r.add_argument("--out", required=True, type=Path); r.add_argument("--no-figures", action="store_true")
    r.add_argument("--dpi", type=int, help="PNG resolution; default figures.yaml png_dpi")
    r = sub.add_parser("plot", help="Render Figure 1 and, when tables/si exists, both SI figures from a reproduction's tables")
    r.add_argument("--tables", type=Path, required=True); r.add_argument("--out", type=Path, required=True)
    r.add_argument("--dpi", type=int, help="PNG resolution; default figures.yaml png_dpi")
    v = sub.add_parser("validate", help="Check released bundles, source structures and sampling checksums")
    v.add_argument("--out", required=True, type=Path)
    a = sub.add_parser("analyze", help="Compute continuum tables for a modified registry; no comparison with published results")
    a.add_argument("--out", required=True, type=Path)
    a = sub.add_parser("audit", help="Independent arithmetic checks using standard-library Python")
    a.add_argument("--out", required=True, type=Path)
    e = sub.add_parser("extract", help="Extract newly supplied Gaussian logs into a fresh dataset")
    e.add_argument("logs", nargs="+", type=Path); e.add_argument("--out", required=True, type=Path)
    b = sub.add_parser("build", help="Construct the configured neutral linkages; requires the structures extra")
    b.add_argument("--registry", type=Path); b.add_argument("--out", required=True, type=Path)
    b = sub.add_parser("acquire", help="Fetch the configured PubChem CIDs; requires network")
    b.add_argument("--registry", type=Path); b.add_argument("--out", required=True, type=Path)
    g = sub.add_parser("prepare-gaussian", help="Prepare default-Opt rerun inputs from bundled optimized geometries")
    g.add_argument("--out", required=True, type=Path); g.add_argument("--cpus", type=int, default=8); g.add_argument("--memory", default="6GB")
    g = sub.add_parser("prepare-input", help="Prepare a new molecule or a transferred optimized geometry at a configured method")
    g.add_argument("geometry", type=Path); g.add_argument("--method", required=True); g.add_argument("--name", required=True)
    g.add_argument("--out", required=True, type=Path); g.add_argument("--charge", type=int, default=0); g.add_argument("--multiplicity", type=int, default=1)
    g.add_argument("--cpus", type=int, default=8); g.add_argument("--memory", default="6GB")
    c = sub.add_parser("prepare-crest", help="Prepare a CREST command; does not execute CREST")
    c.add_argument("seed", type=Path); c.add_argument("--out", required=True, type=Path)
    c.add_argument("--charge", type=int, default=0); c.add_argument("--multiplicity", type=int, default=1); c.add_argument("--threads", type=int, default=8)
    c = sub.add_parser("prepare-refinement", help="Select from an existing CREST ensemble and prepare DFT starts")
    c.add_argument("seed", type=Path); c.add_argument("ensemble", type=Path)
    c.add_argument("--moiety", required=True, help="A moiety with a selection rule in config/study.json sampling")
    c.add_argument("--out", required=True, type=Path); c.add_argument("--charge", type=int, default=0); c.add_argument("--multiplicity", type=int, default=1)
    c = sub.add_parser("prepare-waters", help="Construct the primary two-water split-site starts")
    c.add_argument("geometry", type=Path); c.add_argument("--family", required=True)
    c.add_argument("--folded", action="store_true", help="Use selected folded A/C construction rules")
    c.add_argument("--out", required=True, type=Path)
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
            fresh(out)
            config = yaml.safe_load((root / "config/figures.yaml").read_text())
            dpi = args.dpi or int(config["figure1"]["png_dpi"])
            render_figure1(args.tables.resolve(), out, dpi, config)
            if (args.tables / "si").is_dir():
                render_si_microsolvation(args.tables.resolve(), out, dpi, config)
                render_si_ranking_occupancy(args.tables.resolve(), out, dpi, config)
        elif args.command == "analyze":
            from .analysis.workflow import run_continuum
            fresh(out); run_continuum(root, out)
        elif args.command == "validate":
            from .validation.inventory import validate
            validate(root, out)
        elif args.command == "audit":
            from .validation.arithmetic import audit
            audit(root, out)
        elif args.command == "extract":
            from .extraction.workflow import extract
            extract(args.logs, out)
        elif args.command == "build":
            from .building.molecules import build_from_registry
            build_from_registry(root, args.registry or root / "config/construction.csv", out)
        elif args.command == "acquire":
            from .acquisition.pubchem import download_by_cid
            download_by_cid(args.registry or root / "config/pubchem.csv", out)
        elif args.command == "prepare-gaussian":
            from .building.gaussian import prepare_study_inputs
            prepare_study_inputs(root, out, args.cpus, args.memory)
        elif args.command == "prepare-input":
            from .building.gaussian import prepare_input
            prepare_input(root, args.geometry, args.method, args.name, out, args.charge, args.multiplicity, args.cpus, args.memory)
        elif args.command == "prepare-crest":
            from .sampling.prepare import prepare_crest
            prepare_crest(args.seed, out, args.charge, args.multiplicity, args.threads)
        elif args.command == "prepare-refinement":
            from .sampling.prepare import prepare_refinement, sampling_policy
            prepare_refinement(args.seed, args.ensemble, args.moiety, out, sampling_policy(root), args.charge, args.multiplicity)
        elif args.command == "prepare-waters":
            from .sampling.prepare import prepare_waters
            prepare_waters(args.geometry, args.family, out, args.folded)
    except (ValueError, RuntimeError, FileExistsError, FileNotFoundError, KeyError) as exc:
        p.exit(1, f"Study workflow failed: {exc}\n")
    print(f"Completed {args.command}: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
