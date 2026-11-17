#!/usr/bin/env bash
# Run Custom-Mode DefenseFinder with mode-selected default paths.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
INPUT="$PROJECT_ROOT/output/prodigal_raw"
REFERENCE="$PROJECT_ROOT/output/step2_raw/eggnog_filtered"
OUTPUT="$PROJECT_ROOT/output/step2_raw/defensefinder_results"
MODELS_DIR=""
ORIGINAL_ARGS=("$@")

usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_defensefinder.sh --mode custom [OPTIONS]

Required:
  --mode custom        Run DefenseFinder for Custom Mode

Optional:
  --input DIR          Prodigal protein directory (default: output/prodigal_raw)
  --reference DIR      eggNOG-filtered mapping directory
  --output DIR         DefenseFinder output directory
  --models_dir DIR     Override the DefenseFinder models directory
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --models_dir) MODELS_DIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" ]] || { echo "ERROR: run_defensefinder.sh is only for Custom Mode" >&2; usage; exit 2; }

if [[ "${CONDA_DEFAULT_ENV:-}" != "defensefinder" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defensefinder bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi

source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
source "$PROJECT_ROOT/scripts/stage_state.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
mkdir -p "$MARKER_DIR"

log_script_start "Step2/S04-Custom" "CPU"
stop_if_no_candidates "Step2/S04-Custom" "input" "eggnog_retained_rows" csv --path "$REFERENCE" --required-column mining_protein_id

# Each current contig validates its own input and result before reuse.
invalidate_stages 04_defensefinder
COMMAND=(python -u "$SCRIPT_DIR/S04_defensefinder_run.py" --input "$INPUT" --reference "$REFERENCE" --output "$OUTPUT")
if [[ -n "$MODELS_DIR" ]]; then
    COMMAND+=(--models_dir "$MODELS_DIR")
fi
run_stage "Step2/S04-Custom" --device CPU -- "${COMMAND[@]}"
[[ -d "$OUTPUT" ]] || { echo "ERROR: DefenseFinder output was not created: $OUTPUT" >&2; exit 1; }
mark_stage_done 04_defensefinder

