"""Convenience wrapper over the same tested public command."""
from pathlib import Path
import sys
from atcgu.cli import main

raise SystemExit(main(["--root", str(Path(__file__).resolve().parents[1]), "reproduce", *sys.argv[1:]]))
