#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
crest seed.xyz --gfn2 --alpb water --ewin 6.0 --T 16 --chrg 0 --uhf 0
