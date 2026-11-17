#!/usr/bin/env bash
# Remove candidates already annotated by DefenseFinder.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
INPUT="$PROJECT_ROOT/output/step2_raw/eggnog_filtered"
DEFENSEFINDER="$PROJECT_ROOT/output/step2_raw/defensefinder_results"
OUTPUT="$PROJECT_ROOT/output/step2_raw/defensefinder_filtered"
FINAL_OUTPUT="$PROJECT_ROOT/output"
ORIGINAL_ARGS=("$@")

usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_defensefinder_filter.sh --mode custom|replication [OPTIONS]

Required:
  --mode MODE             Select Custom or Replication state

Optional:
  --input DIR             eggNOG-filtered mapping directory
  --defensefinder DIR     DefenseFinder result directory
  --output DIR            Filtered working directory
  --final_output DIR      Final summary directory (default: output)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --defensefinder) DEFENSEFINDER="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --final_output) FINAL_OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }

if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi

source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
source "$PROJECT_ROOT/scripts/stage_state.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
mkdir -p "$MARKER_DIR"

log_script_start "Step2/S05" "CPU"
stop_if_no_candidates "Step2/S05" "input" "eggnog_retained_rows" csv --path "$INPUT" --required-column mining_protein_id
[[ -d "$DEFENSEFINDER" ]] || { echo "ERROR: DefenseFinder result directory not found: $DEFENSEFINDER" >&2; exit 1; }

if stage_done 05_novelty_filter && \
    python "$SCRIPT_DIR/step2_candidate_utils.py" check --marker "$(stage_marker 05_novelty_filter)" \
        --input "$INPUT" --input "$DEFENSEFINDER" --output "$OUTPUT" \
        --output "$FINAL_OUTPUT/isozyme_candidates_novel.csv" \
        --output "$FINAL_OUTPUT/isozyme_candidates_defensefinder_annotated.csv" --mirror "$INPUT" "$OUTPUT"; then
    log_skip "Step2/S05" --device CPU
else
    invalidate_stages 05_novelty_filter 06_classify
    python "$SCRIPT_DIR/step2_candidate_utils.py" clear-classification --root "$FINAL_OUTPUT"
    run_stage "Step2/S05" --device CPU -- python -u "$SCRIPT_DIR/S05_defensefinder_filter.py" \
        --input "$INPUT" --defensefinder "$DEFENSEFINDER" \
        --output "$OUTPUT" --final_output "$FINAL_OUTPUT"
fi

NOVEL="$FINAL_OUTPUT/isozyme_candidates_novel.csv"
ROWS="$(candidate_count csv --path "$NOVEL" --required-column mining_protein_id)"
echo "[Step2/S05] SUMMARY novel_mapping_rows=$ROWS"
mark_stage_done 05_novelty_filter
if [[ "$ROWS" -eq 0 ]]; then
    echo "[Step2/S05] STOP mode=$MODE reason=no_candidates after=defensefinder_filter"
fi
