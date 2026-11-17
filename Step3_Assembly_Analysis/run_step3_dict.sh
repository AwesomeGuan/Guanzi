#!/usr/bin/env bash
# ============================================================
# run_step3_dict.sh - Build DefenseFinder dictionaries
#
# Usage:
# bash Step3_Assembly_Analysis/run_step3_dict.sh \
# --input output/step2_raw/defensefinder_results \
# --output output/step3_raw
# ============================================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT="$PROJECT_ROOT/output/step2_raw/defensefinder_results"
OUTPUT="$PROJECT_ROOT/output/step3_raw"

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
 --input) INPUT="$2"; shift 2 ;;
 --output) OUTPUT="$2"; shift 2 ;;
 -h|--help) echo "Usage: bash Step3_Assembly_Analysis/run_step3_dict.sh [--input PATH] [--output PATH]"; exit 0 ;;
 *) echo "Unknown option: $1" >&2; exit 2 ;;
 esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


log() { echo "[Step3-Dict] $*"; }
log "=== Step3: Building DefenseFinder Dictionaries ==="
log "Input: $INPUT"
log "Output: $OUTPUT"
mkdir -p "$OUTPUT"

log "S01: genes dict (hmmer) ..."
run_stage "Step3/S01" --device CPU -- python -u "$SCRIPT_DIR/S01_build_genes_dict.py" \
 --defensefinder "$INPUT" \
 --output "$OUTPUT/defensefinder_genes_dict.pkl"

log "S02: systems dict ..."
run_stage "Step3/S02" --device CPU -- python -u "$SCRIPT_DIR/S02_build_systems_dict.py" \
 --defensefinder "$INPUT" \
 --output "$OUTPUT/defensefinder_systems_dict.pkl"

log "S03: incomplete genes dict ..."
run_stage "Step3/S03" --device CPU -- python -u "$SCRIPT_DIR/S03_build_incomplete_genes_dict.py" \
 --genes_dict "$OUTPUT/defensefinder_genes_dict.pkl" \
 --systems_dict "$OUTPUT/defensefinder_systems_dict.pkl" \
 --output "$OUTPUT/defensefinder_incomplete_dict.pkl"

log "=== Done: $OUTPUT/ ==="
