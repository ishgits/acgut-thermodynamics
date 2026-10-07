#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 || ! -f "$1" || "$1" != *.com ]]; then
  echo 'Usage: bash run_gaussian.sh /path/to/reviewed_input.com' >&2
  exit 2
fi
input_dir=$(cd -- "$(dirname -- "$1")" && pwd)
input_name=$(basename -- "$1")
cd -- "$input_dir"
log_name="${input_name%.com}.log"
if [[ -e "$log_name" ]]; then
  echo "Refusing to overwrite $log_name" >&2
  exit 1
fi
command -v g16 >/dev/null || { echo "g16 not found; load Gaussian 16" >&2; exit 127; }
g16 < "$input_name" > "$log_name"
