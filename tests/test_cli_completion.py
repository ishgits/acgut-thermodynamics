"""CLI completion messages: actual counts, and plot's supported partial run."""
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from atcgu.cli import main
from atcgu.runner import run

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def tables(tmp_path_factory):
    out = run(ROOT, tmp_path_factory.mktemp("reproduction"), figures=False)
    return out / "tables"


def _copy_tables(tables, tmp_path):
    target = tmp_path / "tables"
    shutil.copytree(tables, target)
    return target


def test_reproduce_no_figures_reports_actual_count(root, tmp_path, capsys):
    out = tmp_path / "reproduction"
    assert main(["--root", str(root), "reproduce", "--no-figures", "--out", str(out)]) == 0
    text = capsys.readouterr().out
    assert "25/25 published-result comparisons passed." in text
    assert f"Report: {out.resolve() / 'VALIDATION_REPORT.md'}" in text


def test_validate_reports_counts_and_paths(root, tmp_path, capsys):
    out = tmp_path / "validation"
    assert main(["--root", str(root), "validate", "--out", str(out)]) == 0
    total = json.loads((out / "summary.json").read_text())["checks"]
    text = capsys.readouterr().out
    assert f"{total}/{total} checks passed." in text
    assert f"Summary: {out.resolve() / 'summary.json'}" in text
    assert f"Detailed checks: {out.resolve() / 'record_checks.csv'}" in text


def test_audit_reports_counts_and_paths(root, tmp_path, capsys):
    out = tmp_path / "audit"
    assert main(["--root", str(root), "audit", "--out", str(out)]) == 0
    total = json.loads((out / "summary.json").read_text())["checks"]
    text = capsys.readouterr().out
    assert f"{total}/{total} checks passed." in text
    assert f"Summary: {out.resolve() / 'summary.json'}" in text
    assert f"Detailed checks: {out.resolve() / 'independent_arithmetic_checks.csv'}" in text


def test_plot_full_reports_all_figures(tables, tmp_path, capsys):
    out = tmp_path / "figures"
    args = ["--root", str(ROOT), "plot", "--tables", str(_copy_tables(tables, tmp_path)),
            "--out", str(out), "--dpi", "72"]
    assert main(args) == 0
    text = capsys.readouterr().out
    assert "Generated Figure 1, SI entropy figure and Figure S1." in text
    assert f"Output: {out.resolve()}" in text
    for name in ["figure1.png", "si_microsolvation_sensitivity.png", "decomposition.png"]:
        assert (out / name).is_file()


def test_plot_missing_si_is_partial_and_preserves_outputs(tables, tmp_path, capsys):
    copied = _copy_tables(tables, tmp_path)
    shutil.rmtree(copied / "si")
    out = tmp_path / "figures"
    args = ["--root", str(ROOT), "plot", "--tables", str(copied), "--out", str(out), "--dpi", "72"]
    assert main(args) == 0
    text = capsys.readouterr().out
    assert "Partial completion: generated Figure 1 and Figure S1." in text
    assert "Skipped SI entropy figure: missing tables/si/." in text
    assert f"Output: {out.resolve()}" in text
    assert (out / "figure1.png").is_file() and (out / "decomposition.png").is_file()
    assert not (out / "si_microsolvation_sensitivity.png").exists()


def test_plot_invalid_si_is_not_reported_as_skipped(tables, tmp_path, capsys):
    copied = _copy_tables(tables, tmp_path)
    path = copied / "si" / "05_fixed_reference_placement_scores.csv"
    frame = pd.read_csv(path)
    frame.iloc[:-1].to_csv(path, index=False)
    out = tmp_path / "figures"
    args = ["--root", str(ROOT), "plot", "--tables", str(copied), "--out", str(out), "--dpi", "72"]
    with pytest.raises(SystemExit) as error:
        main(args)
    assert error.value.code == 1
    captured = capsys.readouterr()
    assert "Skipped SI entropy figure" not in captured.out + captured.err
    assert "Study workflow failed" in captured.err
    assert (out / "figure1.png").is_file()
