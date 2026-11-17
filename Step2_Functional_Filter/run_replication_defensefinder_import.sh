#!/usr/bin/env bash
# Import the precomputed GTDB DefenseFinder results used by Replication Mode.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
SOURCE="$PROJECT_ROOT/data/gtdb/defensefinder_results"
MAPPING_DIR="$PROJECT_ROOT/output/rep_step1_mapping"
PRODIGAL="$PROJECT_ROOT/output/rep_prodigal_raw"
OUTPUT="$PROJECT_ROOT/output/step2_raw/defensefinder_results"
ORIGINAL_ARGS=("$@")

usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_replication_defensefinder_import.sh --mode replication [OPTIONS]

Required:
  --mode replication    Import the Replication-Mode reference results

Optional:
  --source DIR          Precomputed GTDB DefenseFinder tree
  --mapping_dir DIR     Replication Step 1 mapping directory
  --prodigal DIR        Materialized GTDB Prodigal directory
  --output DIR          DefenseFinder output directory
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --source) SOURCE="$2"; shift 2 ;;
        --mapping_dir) MAPPING_DIR="$2"; shift 2 ;;
        --prodigal) PRODIGAL="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "replication" ]] || { echo "ERROR: this entry point is only for Replication Mode" >&2; usage; exit 2; }

if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi

source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
source "$PROJECT_ROOT/scripts/stage_state.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
mkdir -p "$MARKER_DIR"

log_script_start "Step2/S04-Replication" "CPU"
stop_if_no_candidates "Step2/S04-Replication" "input" "mapping_rows" csv --path "$MAPPING_DIR" --required-column mining_protein_id

if stage_done 04_defensefinder && \
    python "$SCRIPT_DIR/step2_candidate_utils.py" check --marker "$(stage_marker 04_defensefinder)" \
        --input "$SOURCE" --input "$MAPPING_DIR" --input "$PRODIGAL" --output "$OUTPUT" && \
    python -c 'import importlib,sys; sys.path.insert(0,sys.argv[1]); m=importlib.import_module("S04_replication_import_defensefinder"); sys.exit(0 if m.import_complete(*sys.argv[2:]) else 1)' \
        "$SCRIPT_DIR" "$SOURCE" "$MAPPING_DIR" "$OUTPUT"; then
    log_skip "Step2/S04-Replication" --device CPU
else
    invalidate_stages 04_defensefinder
    run_stage "Step2/S04-Replication" --device CPU -- python -u "$SCRIPT_DIR/S04_replication_import_defensefinder.py" \
        --source "$SOURCE" --mapping_dir "$MAPPING_DIR" \
        --prodigal "$PRODIGAL" --output "$OUTPUT"
    [[ -d "$OUTPUT" ]] || { echo "ERROR: imported DefenseFinder output was not created: $OUTPUT" >&2; exit 1; }
    mark_stage_done 04_defensefinder
fi

