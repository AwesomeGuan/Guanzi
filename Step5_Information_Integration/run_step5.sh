#!/usr/bin/env bash
# Run the complete Step 5 information-integration chapter.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
PRODIGAL_DIR=""
TMALIGN_DIR="$PROJECT_ROOT/output/step4_output"
STEP5_OUTPUT="$PROJECT_ROOT/output/step5_output"

usage() {
    cat <<EOF
Usage: bash Step5_Information_Integration/run_step5.sh --mode custom|replication [OPTIONS]
  --prodigal_dir PATH  Override the mode-specific Prodigal directory
  --tmalign_dir PATH   Step 4 formal output (default: output/step4_output)
  --output PATH        Final output directory (default: output/step5_output)
EOF
}

ORIGINAL_ARGS=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --prodigal_dir) PRODIGAL_DIR="$2"; shift 2 ;;
        --tmalign_dir) TMALIGN_DIR="$2"; shift 2 ;;
        --output) STEP5_OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }
if [[ "$MODE" == "custom" ]]; then
    PRODIGAL_DIR="${PRODIGAL_DIR:-$PROJECT_ROOT/output/prodigal_raw}"
else
    PRODIGAL_DIR="${PRODIGAL_DIR:-$PROJECT_ROOT/output/rep_prodigal_raw}"
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step5" "CPU" "${GPUS:-}"

source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step5"
source "$PROJECT_ROOT/scripts/stage_state.sh"
mkdir -p "$MARKER_DIR" "$PROJECT_ROOT/output/step5_raw/gene_location" "$STEP5_OUTPUT"
LOCATIONS="$PROJECT_ROOT/output/step5_raw/gene_location/component_locations.csv"

[[ -s "$TMALIGN_DIR/comparison_qualified.csv" ]] \
    && [[ -s "$TMALIGN_DIR/comparison_disqualified.csv" ]] \
    || { echo "ERROR: Step 4 structural evidence is incomplete under $TMALIGN_DIR" >&2; exit 1; }
stop_if_no_candidates "Step5" "input" "formal_instance_rows" csv --path "$TMALIGN_DIR/assembly_qualified.csv" --path "$TMALIGN_DIR/assembly_disqualified.csv" --required-column mining_contig_system_index
[[ -d "$PRODIGAL_DIR" ]] || { echo "ERROR: Prodigal directory not found: $PRODIGAL_DIR" >&2; exit 1; }

# The location table is reused only while it still covers every Step 4 component
# key, so a Step 4 rerun that introduces new instances rebuilds it instead of
# silently reusing an incomplete table.
if stage_done 01_locations && [[ -s "$LOCATIONS" ]] \
        && python "$SCRIPT_DIR/S01_build_component_locations.py" \
        --step4 "$TMALIGN_DIR" --prodigal_dir "$PRODIGAL_DIR" \
        --output "$LOCATIONS" --check >/dev/null 2>&1; then
    log_skip "Step5/S01" --device CPU
    echo "Step 5 locations: already complete"
else
    invalidate_stages 01_locations
    run_stage "Step5/S01" --device CPU -- python -u "$SCRIPT_DIR/S01_build_component_locations.py" \
        --step4 "$TMALIGN_DIR" --prodigal_dir "$PRODIGAL_DIR" --output "$LOCATIONS"
    [[ -s "$LOCATIONS" ]] || exit 1
    mark_stage_done 01_locations
fi

[[ -s "$LOCATIONS" ]] || { echo "ERROR: location table not found: $LOCATIONS" >&2; exit 1; }
COMPLETE=1
for kind in mapping assembly comparison; do
    for result in qualified disqualified; do
        [[ -s "$STEP5_OUTPUT/${kind}_${result}.csv" ]] || COMPLETE=0
    done
done
if stage_done 02_finalize && [[ "$COMPLETE" -eq 1 ]]; then
    log_skip "Step5/S02" --device CPU
    echo "Step 5 finalization: already complete"
else
    invalidate_stages 02_finalize
    run_stage "Step5/S02" --device CPU -- python -u "$SCRIPT_DIR/S02_finalize_syntax.py" \
        --step4 "$TMALIGN_DIR" --locations "$LOCATIONS" --output "$STEP5_OUTPUT"
    for kind in mapping assembly comparison; do
        for result in qualified disqualified; do
            [[ -s "$STEP5_OUTPUT/${kind}_${result}.csv" ]] || exit 1
        done
    done
    mark_stage_done 02_finalize
fi

echo "Step 5 complete for $MODE mode: $STEP5_OUTPUT"
