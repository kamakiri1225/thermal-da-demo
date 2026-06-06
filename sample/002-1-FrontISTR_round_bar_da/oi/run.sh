#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export FISTR1="${FISTR1:-/home/kamakiri/local/frontistr/bin/fistr1}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-frontistr-da}"

python3 da_main.py
