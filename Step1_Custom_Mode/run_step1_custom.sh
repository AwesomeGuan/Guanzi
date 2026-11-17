#!/usr/bin/env bash
# ============================================================
# Custom Mode downstream processing (S02-S04)
#
# Reads FoldSeek raw output from --input, runs filtering + mapping.
# The FoldSeek search (S01) must already be complete.
#
# Usage:
# bash run_step1_custom.sh --input output/step1_raw/foldseek_search_results --output output/step1_mapping
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT_DIR=""
OUTPUT_DIR="$PROJECT_ROOT/output/step1_mapping"
RAW_DIR="$PROJECT_ROOT/output/step1_raw"
MIN_QCOV=0.70
MIN_TCOV=0.70
FIDENT_BLACKLIST_THRESHOLD=0.30

usage() {
 cat <<EOF
Usage: bash run_step1_custom.sh [OPTIONS]

Options:
 --input PATH FoldSeek raw output from S01 (required)
 --output PATH Final mapping output dir (default: output/step1_mapping)
 --min_qcov FLOAT Minimum FoldSeek query coverage (default: 0.70)
 --min_tcov FLOAT Minimum FoldSeek target coverage (default: 0.70)
 --fident_blacklist_threshold FLOAT Target blacklist boundary (default: 0.30)
EOF
 exit 1
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
 --input) INPUT_DIR="$2"; shift 2 ;;
 --output) OUTPUT_DIR="$2"; shift 2 ;;
 --min_qcov) MIN_QCOV="$2"; shift 2 ;;
 --min_tcov) MIN_TCOV="$2"; shift 2 ;;
 --fident_blacklist_threshold) FIDENT_BLACKLIST_THRESHOLD="$2"; shift 2 ;;
 -h|--help) usage ;;
 *) echo "Unknown: $1"; usage ;;
 esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=custom
source "$PROJECT_ROOT/scripts/stage_logging.sh"


if [[ -z "$INPUT_DIR" ]]; then
 echo "ERROR: --input is required"
 usage
fi
if [[ ! -d "$INPUT_DIR" ]]; then
 echo "ERROR: $INPUT_DIR not found. Run run_foldseek_search.sh first."
 exit 1
fi

FIDENT_OUT="$RAW_DIR/fident_filter"
COVERAGE_OUT="$RAW_DIR/coverage_filter"

log() { echo "[Custom] $*"; }

log "=== Custom Mode: S02-S04 ==="
log " Input: $INPUT_DIR"
log " Output: $OUTPUT_DIR"

# --- S02: coverage filter ---
log "[1/3] coverage filter ..."
run_stage "Step1-Custom/S02" --device CPU --info "qcov>=${MIN_QCOV:-0.70} tcov>=${MIN_TCOV:-0.70}" -- python -u "$SCRIPT_DIR/S02_coverage_filter.py" \
    --input "$INPUT_DIR" \
    --output "$COVERAGE_OUT" \
    --min_qcov "$MIN_QCOV" \
    --min_tcov "$MIN_TCOV"

# --- S03: target-level fident blacklist and deduplication ---
log "[2/3] target-level fident blacklist and target deduplication ..."
run_stage "Step1-Custom/S03" --device CPU --info "target_blacklist: any fident>${FIDENT_BLACKLIST_THRESHOLD:-0.30}; then max-TM-score deduplication" -- python -u "$SCRIPT_DIR/S03_fident_filter.py" \
    --input "$COVERAGE_OUT" \
    --output "$FIDENT_OUT" \
    --fident_blacklist_threshold "$FIDENT_BLACKLIST_THRESHOLD"

# --- S04: Generate final mapping ---
log "[3/3] Generate canonical Query/Mining protein mapping ..."
run_stage "Step1-Custom/S04" --device CPU -- python -u "$SCRIPT_DIR/S04_custom_mapping.py" \
 --isozyme_dir "$FIDENT_OUT" \
 --output_dir "$OUTPUT_DIR"

log "=== Custom Mode downstream processing complete ==="
echo "Mapping: $OUTPUT_DIR/"
