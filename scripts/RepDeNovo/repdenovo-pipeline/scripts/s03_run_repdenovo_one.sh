#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --repdenovo-dir /path/to/REPdenovo --config configs/generated/S1.conf --reads reads/S1.txt [--python /usr/bin/python]"
  exit 1
}

PYTHON_BIN="python"
REPDENOVO_DIR=""
CONFIG=""
READS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repdenovo-dir) REPDENOVO_DIR="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --reads) READS="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    *) usage ;;
  esac
done

[[ -n "$REPDENOVO_DIR" && -n "$CONFIG" && -n "$READS" ]] || usage

MAIN="${REPDENOVO_DIR%/}/main.py"
[[ -f "$MAIN" ]] || { echo "ERROR: main.py not found at $MAIN"; exit 2; }

mkdir -p logs
ts="$(date -u +%Y%m%dT%H%M%SZ)"

echo "[${ts}] Running REPdenovo Assembly..."
"$PYTHON_BIN" "$MAIN" -c Assembly -g "$CONFIG" -r "$READS" 2>&1 | tee "logs/run_assembly_${ts}.log"

# If your REPdenovo version uses a different stage name, adjust here.
# echo "[${ts}] Running REPdenovo Scaffolding..."
# "$PYTHON_BIN" "$MAIN" -c Scaffolding -g "$CONFIG" -r "$READS" 2>&1 | tee "logs/run_scaffolding_${ts}.log"

echo "Done."