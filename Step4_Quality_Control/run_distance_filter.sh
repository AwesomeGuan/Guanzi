#!/usr/bin/env bash
# Stage 1: Genomic distance filtering (S01)
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT="$PROJECT_ROOT/output/step3_output"
RULES="$PROJECT_ROOT/data/DefenseFinder_rules"
OUTPUT="$PROJECT_ROOT/output/step4_raw/distance_filter"
MAX_GENE_GAP="10"

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
        --mode) MODE="$2"; shift 2 ;;
 --input) INPUT="$2"; shift 2 ;;
 --rules) RULES="$2"; shift 2 ;;
 --output) OUTPUT="$2"; shift 2 ;;
 --max_gene_gap) MAX_GENE_GAP="$2"; shift 2 ;;
 -h|--help)
  echo "Usage: bash Step4_Quality_Control/run_distance_filter.sh [options] [--mode custom|replication]"
  echo "  --input PATH          Step3 output directory"
  echo "  --rules PATH          DefenseFinder rules directory"
  echo "  --output PATH         Distance-filter output directory"
  echo "  --max_gene_gap N      Maximum adjacent gene-index difference (default: 10)"
  exit 0 ;;
 *) echo "Unknown argument: $1" >&2; exit 2 ;;
 esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


log() { echo "[Distance] $*"; }
log "=== Genomic Distance Filtering ==="
[[ -d "$INPUT" ]] || { echo "ERROR: Step3 output directory not found: $INPUT" >&2; exit 1; }
[[ -d "$RULES" ]] || { echo "ERROR: rules directory not found: $RULES" >&2; exit 1; }
run_stage "Step4/S01" --device CPU --info "max_gene_gap=$MAX_GENE_GAP" -- python -u "$SCRIPT_DIR/S01_distance_filter.py" \
  --input "$INPUT" \
  --rules "$RULES" \
  --output "$OUTPUT" \
  --max_gene_gap "$MAX_GENE_GAP"
log "Done: $OUTPUT/"
