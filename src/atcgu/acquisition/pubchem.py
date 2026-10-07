"""Study acquisition adapter. Adapted from doi:10.5281/zenodo.21365151.

Use explicit, configured CIDs. No live requests occur during ordinary reproduction.
"""
from __future__ import annotations
import csv
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

from ..records import digest
from ..paths import fresh


def download_by_cid(registry: Path, output: Path, retries: int = 4) -> Path:
    fresh(output)
    manifest = []
    for row in csv.DictReader(registry.open()):
        cid = int(row["cid"])
        if cid < 1:
            raise ValueError("PubChem CID must be a positive integer")
        name = row.get("species_id") or f"cid_{cid}"
        if any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in name):
            raise ValueError("Invalid species ID")
        error = None
        for record_type in ["3d", "2d"]:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type={record_type}"
            for attempt in range(retries):
                try:
                    time.sleep(0.25)
                    with urlopen(Request(url, headers={"User-Agent": "atcgu-study-reproduction/0.1"}), timeout=60) as response:
                        body = response.read()
                    if b"M  END" not in body:
                        raise ValueError("Response is not an SDF")
                    break
                except Exception as exc:
                    error = exc
                    if attempt + 1 < retries:
                        time.sleep(2 ** attempt)
            else:
                continue
            path = output / f"{name}.sdf"
            if path.exists():
                raise ValueError("Duplicate acquisition species ID")
            path.write_bytes(body)
            manifest.append({"species_id": name, "cid": cid, "record_type": record_type,
                             "source_url": url, "sha256": digest(path), "retrieved_unix": time.time()})
            break
        else:
            raise RuntimeError(f"PubChem download failed for CID {cid}") from error
    (output / "acquisition_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return output
