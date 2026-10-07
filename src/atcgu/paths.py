"""Output-directory guards shared by every command."""
from __future__ import annotations
from pathlib import Path


def check_output(root: Path, output: Path) -> Path:
    """Inside the checkout, output may only go under ``outputs/``; elsewhere is fine."""
    root, output = root.resolve(), output.resolve()
    if output.is_relative_to(root) and not output.is_relative_to(root / "outputs"):
        raise ValueError(f"Inside the checkout, write outputs under {root / 'outputs'}")
    return output


def fresh(path: Path) -> Path:
    """Create an empty output directory; refuse to reuse a non-empty one."""
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Choose a fresh output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path
