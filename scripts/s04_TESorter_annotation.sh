#!/usr/bin/env bash
set -euo pipefail

# Run TEsorter on a single TE FASTA (plants database)
# Usage:
#   bash scripts/run_tesorter_plants.sh path/to/TE.fasta
#
# Output:
#   results/tesorter/<basename>/

FASTA="${1:-}"

if [[ -z "$FASTA" ]]; then
  echo "Usage: $0 <TE_fasta>" >&2
  exit 1
fi

if [[ ! -f "$FASTA" ]]; then
  echo "ERROR: FASTA not found: $FASTA" >&2
  exit 1
fi

if ! command -v TEsorter >/dev/null 2>&1; then
  echo "ERROR: TEsorter not found in PATH. Install it (e.g., conda/mamba install -c bioconda tesorter)." >&2
  exit 1
fi

mkdir -p results/tesorter

base="$(basename "$FASTA")"
base="${base%.gz}"
base="${base%.*}"

outdir="results/tesorter/${base}"
mkdir -p "$outdir"

echo "[INFO] Input : $FASTA"
echo "[INFO] Outdir: $outdir"
echo "[INFO] Cmd   : TEsorter \"$FASTA\" -db rexdb-plant"

( cd "$outdir" && TEsorter "$FASTA" -db rexdb-plant )

echo "[OK] Done."