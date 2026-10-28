#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Run PASTEC / PASTEClassifier on a single consensus FASTA
#
# Requirements:
#   - python2.7 available (PASTEC is python2-based in many setups)
#   - PASTEClassifier.py available (local install OR container mount)
#   - a config file (see config/PASTEClassifier.cfg.template)
#
# Usage:
#   bash scripts/run_pastec.sh -i path/to/consensi.fasta -c config/PASTEClassifier.cfg -p /path/to/PASTEClassifier.py
#
# Example:
#   bash scripts/run_pastec.sh -i data/example/consensi.fasta -c config/PASTEClassifier.cfg -p /opt/PASTEC_linux-x64-2.0/bin/PASTEClassifier.py
# ------------------------------------------------------------

FASTA=""
CFG=""
PASTECLASSIFIER_PY=""
OUTDIR="results/pastec"
STEP="0"
VERBOSITY="3"
DO_CLEAN="false"

usage() {
  cat <<EOF
Usage:
  run_pastec.sh -i <consensi.fasta> -c <PASTEClassifier.cfg> -p <PASTEClassifier.py> [--outdir <dir>] [--step 0|1|2] [--verbosity 1-4] [--clean]

Required:
  -i    Input consensus FASTA
  -c    PASTEClassifier configuration file
  -p    Path to PASTEClassifier.py

Optional:
  --outdir     Output directory (default: results/pastec/<basename>/)
  --step       0=all, 1=detect features, 2=classification only (default: 0)  (see tool docs) 
  --verbosity  Verbosity level 1-4 (default: 3)
  --clean      Remove temporary files (passes --clean)
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) FASTA="$2"; shift 2;;
    -c) CFG="$2"; shift 2;;
    -p) PASTECLASSIFIER_PY="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --step) STEP="$2"; shift 2;;
    --verbosity) VERBOSITY="$2"; shift 2;;
    --clean) DO_CLEAN="true"; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1;;
  esac
done

# Checks
[[ -z "$FASTA" || -z "$CFG" || -z "$PASTECLASSIFIER_PY" ]] && { echo "ERROR: missing required args" >&2; usage; exit 1; }
[[ ! -f "$FASTA" ]] && { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
[[ ! -f "$CFG" ]] && { echo "ERROR: config not found: $CFG" >&2; exit 1; }
[[ ! -f "$PASTECLASSIFIER_PY" ]] && { echo "ERROR: PASTEClassifier.py not found: $PASTECLASSIFIER_PY" >&2; exit 1; }

if ! command -v python2.7 >/dev/null 2>&1; then
  echo "ERROR: python2.7 not found in PATH (PASTEC often requires python2.7)." >&2
  exit 1
fi

mkdir -p "$OUTDIR"

base="$(basename "$FASTA")"
base="${base%.gz}"
base="${base%.*}"

run_dir="${OUTDIR}/${base}"
mkdir -p "$run_dir"

log="${run_dir}/pastec.log"

echo "[INFO] Input FASTA : $FASTA"
echo "[INFO] Config      : $CFG"
echo "[INFO] Runner dir  : $run_dir"
echo "[INFO] Step        : $STEP"
echo "[INFO] Verbosity   : $VERBOSITY"
echo "[INFO] Log         : $log"

# Run inside run_dir so outputs stay grouped
cd "$run_dir"

cmd=(python2.7 "$PASTECLASSIFIER_PY" -i "$FASTA" -C "$CFG" -S "$STEP" -v "$VERBOSITY")

if [[ "$DO_CLEAN" == "true" ]]; then
  cmd+=(--clean)
fi

echo "[INFO] Command: ${cmd[*]}" > "$log"
"${cmd[@]}" >> "$log" 2>&1

echo "[OK] Done."