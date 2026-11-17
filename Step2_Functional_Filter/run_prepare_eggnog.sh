#!/usr/bin/env bash
# Prepare the Step 2 FASTA for local eggNOG-mapper annotation.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
INPUT=""
PROTEIN_FAA=""
OUTPUT="$PROJECT_ROOT/output/step2_raw"
ORIGINAL_ARGS=("$@")

usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_prepare_eggnog.sh --mode custom|replication [OPTIONS]

Required:
  --mode MODE          Select the Step 1 input layout

Optional:
  --input DIR          Override the Step 1 mapping directory
  --protein_faa DIR    Override the Prodigal protein directory
  --output DIR         Step 2 working directory (default: output/step2_raw)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --protein_faa) PROTEIN_FAA="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }

if [[ "$MODE" == "custom" ]]; then
    INPUT="${INPUT:-$PROJECT_ROOT/output/step1_mapping}"
    PROTEIN_FAA="${PROTEIN_FAA:-$PROJECT_ROOT/output/prodigal_raw}"
else
    INPUT="${INPUT:-$PROJECT_ROOT/output/rep_step1_mapping}"
    PROTEIN_FAA="${PROTEIN_FAA:-$PROJECT_ROOT/output/rep_prodigal_raw}"
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi

source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
source "$PROJECT_ROOT/scripts/stage_state.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
mkdir -p "$MARKER_DIR" "$OUTPUT"
UPLOAD="$OUTPUT/eggnog_upload.faa"

log_script_start "Step2/S01" "CPU"
stop_if_no_candidates "Step2/S01" "input" "mapping_rows" csv --path "$INPUT" --required-column mining_protein_id

if stage_done 01_prepare_eggnog && [[ -s "$UPLOAD" ]] && \
    python "$SCRIPT_DIR/step2_candidate_utils.py" check --marker "$(stage_marker 01_prepare_eggnog)" \
        --input "$INPUT" --input "$PROTEIN_FAA" --output "$UPLOAD"; then
    log_skip "Step2/S01" --device CPU
else
    invalidate_stages 01_prepare_eggnog 02_eggnog_annotation 03_eggnog_filter 05_novelty_filter 06_classify
    [[ "$MODE" != custom ]] || invalidate_stages 04_defensefinder
    run_stage "Step2/S01" --device CPU -- python -u "$SCRIPT_DIR/S01_prepare_eggNOG.py" \
        --input "$INPUT" --protein_faa "$PROTEIN_FAA" --output "$OUTPUT"
fi

stop_if_no_candidates "Step2/S01" "prepare" "upload_protein_records" fasta --path "$UPLOAD"
mark_stage_done 01_prepare_eggnog
echo "Step 2 eggNOG upload FASTA is ready: $UPLOAD"
