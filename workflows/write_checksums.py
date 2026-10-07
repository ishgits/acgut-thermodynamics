"""Record the SHA-256 of every released file in provenance/file_checksums.csv."""
from pathlib import Path
import csv
import hashlib

ROOT = Path(__file__).resolve().parents[1]
IGNORED = {".git", ".venv", "outputs", "archive", "__pycache__", ".pytest_cache", ".ipynb_checkpoints", "build", "dist"}
TARGET = "provenance/file_checksums.csv"


def released_files(root=ROOT):
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root)
        if not p.is_file() or any(x in IGNORED or x.endswith(".egg-info") for x in rel.parts):
            continue
        if p.name == ".DS_Store" or p.suffix in {".pyc", ".pyo"} or str(rel) == TARGET:
            continue
        if p.suffix in {".log", ".out", ".err", ".chk", ".fchk", ".rwf"}:
            raise ValueError(f"Raw calculation output must not be released: {rel}")
        yield p


if __name__ == "__main__":
    target = ROOT / TARGET
    with target.open("w", newline="") as h:
        w = csv.writer(h); w.writerow(["path", "sha256", "size_bytes"])
        for p in released_files():
            w.writerow([str(p.relative_to(ROOT)), hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_size])
    print(f"Wrote {target.name}. Verify that every changed file was intentional before publication.")
