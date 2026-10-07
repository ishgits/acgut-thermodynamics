"""Verify every released file against provenance/file_checksums.csv, including missing/extra files. Standard library only."""
import csv
import hashlib
from write_checksums import ROOT, TARGET, released_files

if __name__ == "__main__":
    rows = list(csv.DictReader((ROOT / TARGET).open()))
    expected = {r["path"] for r in rows}
    if len(rows) != len(expected):
        raise ValueError("Duplicate release checksum paths")
    actual = {str(p.relative_to(ROOT)) for p in released_files()}
    errors = [f"Missing: {p}" for p in expected - actual] + [f"Unexpected: {p}" for p in actual - expected]
    for row in rows:
        p = (ROOT / row["path"]).resolve()
        if not p.is_relative_to(ROOT) or not p.is_file():
            errors.append(f"Invalid path: {row['path']}"); continue
        if hashlib.sha256(p.read_bytes()).hexdigest() != row["sha256"] or p.stat().st_size != int(row["size_bytes"]):
            errors.append(f"Changed: {row['path']}")
    if errors:
        raise SystemExit("Release integrity failed:\n" + "\n".join(errors))
    print(f"All {len(rows)} released files match their recorded hashes.")
