"""Prepare all included study stages in a fresh directory. No external jobs run."""
from pathlib import Path
import argparse
import csv

from atcgu.cli import repository_root
from atcgu.paths import check_output, fresh
from atcgu.records import within
from atcgu.building.gaussian import prepare_study_inputs
from atcgu.sampling.prepare import prepare_crest, prepare_refinement, prepare_waters, sampling_policy


def prepare(root: Path, output: Path):
    fresh(check_output(root, output))
    prepare_study_inputs(root, output / "gaussian")
    policy = sampling_policy(root)
    with (root / "config/sampling.csv").open() as h:
        for row in csv.DictReader(h):
            name = row["sampling_run_id"].lower()
            seed = within(root, row["seed_path"])
            prepare_crest(seed, output / "crest" / name, int(row["charge"]), int(row["multiplicity"]), threads=int(row["study_threads"]))
            prepare_refinement(seed, within(root, row["ensemble_path"]), row["moiety"], output / "refinement" / name,
                               policy, int(row["charge"]), int(row["multiplicity"]))
    with (root / "config/microsolvation_references.csv").open() as h:
        for row in csv.DictReader(h):
            if row["role"] == "continuum_nucleotide":
                geometry = within(root, row["record_path"]).parent / "optimized.xyz"
                prepare_waters(geometry, row["family"], output / "waters" / row["family"], row["folded_start"] == "True")
    print(f"Prepared all study inputs in {output}; no calculations submitted.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--root"); p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(); prepare(repository_root(args.root), args.out.resolve())
