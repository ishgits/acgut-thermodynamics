"""Water-fragment and phosphate/base-contact checks for the analysis tables."""
from __future__ import annotations
from pathlib import Path

from ..records import read_xyz
from ..structure import analyze_structure


def analyze_microsolvation_structure(path: str | Path) -> dict[str, object]:
    """Classify a cluster geometry; contact counts are numbers of distinct waters."""
    try:
        z, xyz = read_xyz(Path(path))
    except ValueError as exc:
        return {"structural_status": "failed", "structural_failure": str(exc)}
    metrics = analyze_structure(z, xyz)
    if metrics["structural_parse"] == "FAIL":
        return {"structural_status": "failed", "structural_failure": metrics["reason"]}
    return {
        "structural_status": "valid", "structural_failure": "",
        "atom_count_from_geometry": metrics["natoms"],
        "fragment_sizes": metrics["fragment_sizes"],
        "intact_water_count": metrics["intact_water_count"], "direct_fold": metrics["direct_fold"],
        "near_fold": metrics["near_fold"], "fold_state": metrics["fold_state"],
        "phosphate_to_water_hb_count": len({c["water_id"] for c in metrics["phosphate_to_water_contacts"]}),
        "water_to_base_hb_count": len({c["water_id"] for c in metrics["water_to_base_contacts"]}),
        "bridging_water_count": metrics["bridging_water_count"],
    }
