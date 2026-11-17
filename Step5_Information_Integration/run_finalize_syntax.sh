#!/usr/bin/env bash
# Stage 2: Append transcription order and syntax to the Comparison tables (S02).
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

STEP4="$PROJECT_ROOT/output/step4_output"
LOCATIONS="$PROJECT_ROOT/output/step5_raw/gene_location/component_locations.csv"
OUTPUT="$PROJECT_ROOT/output/step5_output"

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
        --mode) MODE="$2"; shift 2 ;;
 --step4) STEP4="$2"; shift 2 ;;
 --locations) LOCATIONS="$2"; shift 2 ;;
 --output) OUTPUT="$2"; shift 2 ;;
 -h|--help)
  echo "Usage: bash Step5_Information_Integration/run_finalize_syntax.sh [options] [--mode custom|replication]"
  echo "  --step4 PATH       Step 4 formal output directory (default: output/step4_output)"
  echo "  --locations FILE   Component-location table written by S01 (default: output/step5_raw/gene_location/component_locations.csv)"
  echo "  --output PATH      Final output directory (default: output/step5_output)"
  exit 0 ;;
 *) echo "Unknown argument: $1" >&2; exit 2 ;;
 esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


log() { echo "[Step5-S02] $*"; }
log "=== Final Tables and Syntax ==="
[[ -d "$STEP4" ]] || { echo "ERROR: Step 4 output directory not found: $STEP4" >&2; exit 1; }
[[ -s "$LOCATIONS" ]] || { echo "ERROR: component-location table not found: $LOCATIONS" >&2; exit 1; }
mkdir -p "$OUTPUT"
run_stage "Step5/S02" --device CPU -- python -u "$SCRIPT_DIR/S02_finalize_syntax.py" \
  --step4 "$STEP4" \
  --locations "$LOCATIONS" \
  --output "$OUTPUT"
log "Done: $OUTPUT/"
