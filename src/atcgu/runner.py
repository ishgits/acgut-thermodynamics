"""One fresh reproduction run, compared with the published results, with code/configuration/data provenance."""
from __future__ import annotations
from datetime import datetime, timezone
import importlib.metadata
import json
import platform
from pathlib import Path
import yaml

from .constants import check_study_conventions
from .paths import check_output, fresh
from .records import digest
from .analysis.workflow import run_continuum
from .analysis.microsolvation import run_primary
from .analysis.microsolvation_si import run_si
from .validation.regression import compare_run

KEY_TABLES = ["tables/07_reaction_energies.csv", "tables/09_relative_rankings.csv",
              "tables/micro_06_primary_selected_candidates.csv", "tables/si/06_selected_candidates.csv",
              "tables/si/05_fixed_reference_placement_scores.csv", "tables/si/08_selected_fixed_reference_shifts.csv",
              "tables/si/09_placement_combinations.csv", "tables/si/10_rank_occupancy.csv",
              "tables/si/thermochemical_summary.csv"]


def run(root: Path, output: Path, *, figures: bool = True, dpi: int | None = None) -> Path:
    started_utc = datetime.now(timezone.utc).isoformat()
    root = root.resolve()
    output = check_output(root, output)
    figure_config = yaml.safe_load((root / "config/figures.yaml").read_text())
    study = json.loads((root / "config/study.json").read_text())
    check_study_conventions(study, figure_config)
    dpi = int(figure_config["figure1"]["png_dpi"]) if dpi is None else dpi
    fresh(output)
    tables = output / "tables"
    tables.mkdir()
    outputs = run_continuum(root, tables)
    primary = run_primary(root, tables, outputs["05_selected_species"])
    outputs.update(run_si(root, tables, primary))
    outputs.update({name: frame for name, frame in primary.items() if name != "records"})
    panel_a = figure_config["figure1"]
    rankings = outputs["09_relative_rankings"]
    rankings.loc[rankings.analysis_id.eq(panel_a["panel_a_analysis_id"]) & rankings.reaction_id.eq(panel_a["reaction_id"])].to_csv(tables / "figure1_panel_a.csv", index=False)
    plotted = si_plotted = decomposition_plotted = None
    if figures:
        from .plotting.figure1 import render_figure1
        from .plotting.si_microsolvation import render_si_microsolvation
        from .plotting.decomposition import render_decomposition
        plotted = render_figure1(tables, output / "figures", dpi, figure_config)
        si_plotted = render_si_microsolvation(tables, output / "figures", dpi, figure_config)
        decomposition_plotted = render_decomposition(tables, output / "figures", dpi, figure_config)
    checks = compare_run(root, outputs, plotted, si_plotted)
    checks.to_csv(output / "reproduction_checks.csv", index=False)
    passed = bool(checks.passed.all())
    hashes = {str(p.relative_to(root)): digest(p) for folder in [root / "src/atcgu", root / "config"]
              for p in sorted(folder.rglob("*")) if p.is_file() and "__pycache__" not in p.parts}
    hashes["data/calculations/manifest.csv"] = digest(root / "data/calculations/manifest.csv")
    for p in sorted((root / "data/published").rglob("*.csv")):
        hashes[str(p.relative_to(root))] = digest(p)
    packages = {}
    for name in ["numpy", "pandas", "matplotlib", "PyYAML"]:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    metadata = {"started_utc": started_utc, "finished_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
                "packages": packages, "input_hashes": hashes, "figure_dpi": dpi if figures else None,
                "source_logs_required": False, "reference_family": study["reference_family"],
                "scope": study["scope"], "comparisons": len(checks), "all_regressions_passed": passed}
    (output / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    valid = outputs["01_input_inventory"].is_usable.sum()
    lines = ["# Reproduction report", "",
             f"{int(checks.passed.sum())} of {len(checks)} published-result comparisons passed." if passed else "Some comparisons FAILED; see the table below.", "",
             f"- Continuum records: {len(outputs['01_input_inventory'])}; usable: {valid}.",
             f"- Reaction quantities: {len(outputs['07_reaction_energies'])}.",
             f"- Primary water candidates: {len(outputs['micro_05_primary_candidate_energies'])}.",
             f"- SI fixed-reference placement scores: {len(outputs['si_05_fixed_reference_placement_scores'])}; placement combinations: {len(outputs['si_09_placement_combinations'])}.",
             "- Coverage, balanced elements/charges, path closure, state/method, complete modes and structural criteria checked.", "",
             "## Comparisons", "", "| Published table | Rows | Max numeric difference | Tolerance | Passed |", "|---|---|---|---|---|"]
    lines += [f"| `{r.check_id}` | {r.rows} | {r.maximum_numeric_difference:.3g} | {r.tolerance:g} | {r.passed} |" for r in checks.itertuples()]
    lines += ["", "## Where to look", ""]
    if figures:
        lines += ["- Figure 1: `figures/figure1.png` (plotted points: `figures/figure1_plotted_data.csv`)",
                  "- SI fixed-reference figure: `figures/si_microsolvation_sensitivity.png` (plotted points: `figures/si_microsolvation_plotted_data.csv`)",
                  "- Decomposition figure: `figures/decomposition.png`"]
    else:
        lines += ["- Figures were skipped (`--no-figures`); every table comparison above still ran."]
    lines += [f"- `{name}`" for name in KEY_TABLES]
    (output / "VALIDATION_REPORT.md").write_text("\n".join(lines) + "\n")
    if not passed:
        raise ValueError("Published-result comparisons failed; inspect reproduction_checks.csv")
    return output
