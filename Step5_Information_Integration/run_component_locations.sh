#!/usr/bin/env bash
# Stage 1: Resolve Mining component coordinates from the Prodigal GFF tree (S01).
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

STEP4="$PROJECT_ROOT/output/step4_output"
PRODIGAL_DIR=""
OUTPUT="$PROJECT_ROOT/output/step5_raw/gene_location/component_locations.csv"

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
        --mode) MODE="$2"; shift 2 ;;
 --step4) STEP4="$2"; shift 2 ;;
 --prodigal_dir) PRODIGAL_DIR="$2"; shift 2 ;;
 --output) OUTPUT="$2"; shift 2 ;;
 -h|--help)
  echo "Usage: bash Step5_Information_Integration/run_component_locations.sh [options] [--mode custom|replication]"
  echo "  --step4 PATH          Step 4 formal output directory (default: output/step4_output)"
  echo "  --prodigal_dir PATH   Prodigal GFF tree (default: mode-specific)"
  echo "  --output FILE         Component-location table (default: output/step5_raw/gene_location/component_locations.csv)"
  exit 0 ;;
 *) echo "Unknown argument: $1" >&2; exit 2 ;;
 esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }
if [[ -z "$PRODIGAL_DIR" ]]; then
    if [[ "$MODE" == "replication" ]]; then
        PRODIGAL_DIR="$PROJECT_ROOT/output/rep_prodigal_raw"
    else
        PRODIGAL_DIR="$PROJECT_ROOT/output/prodigal_raw"
    fi
fi
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


log() { echo "[Step5-S01] $*"; }
log "=== Component Locations ==="
[[ -d "$STEP4" ]] || { echo "ERROR: Step 4 output directory not found: $STEP4" >&2; exit 1; }
[[ -d "$PRODIGAL_DIR" ]] || { echo "ERROR: Prodigal directory not found: $PRODIGAL_DIR" >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT")"
run_stage "Step5/S01" --device CPU -- python -u "$SCRIPT_DIR/S01_build_component_locations.py" \
  --step4 "$STEP4" \
  --prodigal_dir "$PRODIGAL_DIR" \
  --output "$OUTPUT"
log "Done: $OUTPUT"
