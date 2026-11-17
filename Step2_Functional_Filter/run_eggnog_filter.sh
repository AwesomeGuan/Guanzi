#!/usr/bin/env bash
# Filter candidates with the validated local or explicitly supplied annotation table.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
INPUT=""
PRODIGAL=""
EGGNOG_TSV="$PROJECT_ROOT/output/step2_raw/eggnog_results.tsv"
DESCRIPTION_KEYWORDS="$PROJECT_ROOT/data/eggNOG_Description_out_function.csv"
PFAM_KEYWORDS="$PROJECT_ROOT/data/eggNOG_PFAMs_out_function.csv"
OUTPUT="$PROJECT_ROOT/output/step2_raw/eggnog_filtered"
FINAL_OUTPUT="$PROJECT_ROOT/output"
ORIGINAL_ARGS=("$@")

usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_eggnog_filter.sh --mode custom|replication [OPTIONS]

Required:
  --mode MODE                  Select the Step 1 input layout

Optional:
  --input DIR                  Override the Step 1 mapping directory
  --prodigal DIR               Override the Prodigal protein directory
  --eggnog_results FILE        eggNOG annotation table (default: output/step2_raw/eggnog_results.tsv)
  --description_keywords FILE  Description keyword table
  --pfam_keywords FILE         PFAM keyword table
  --output DIR                 Filtered working directory
  --final_output DIR           Final summary directory (default: output)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --prodigal) PRODIGAL="$2"; shift 2 ;;
        --eggnog_results) EGGNOG_TSV="$2"; shift 2 ;;
        --description_keywords) DESCRIPTION_KEYWORDS="$2"; shift 2 ;;
        --pfam_keywords) PFAM_KEYWORDS="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --final_output) FINAL_OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }

if [[ "$MODE" == "custom" ]]; then
    INPUT="${INPUT:-$PROJECT_ROOT/output/step1_mapping}"
    PRODIGAL="${PRODIGAL:-$PROJECT_ROOT/output/prodigal_raw}"
else
    INPUT="${INPUT:-$PROJECT_ROOT/output/rep_step1_mapping}"
    PRODIGAL="${PRODIGAL:-$PROJECT_ROOT/output/rep_prodigal_raw}"
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi

source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
source "$PROJECT_ROOT/scripts/stage_state.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
mkdir -p "$MARKER_DIR"

log_script_start "Step2/S03" "CPU"
[[ -s "$EGGNOG_TSV" ]] || { echo "ERROR: eggNOG result not found or empty: $EGGNOG_TSV" >&2; exit 1; }
stop_if_no_candidates "Step2/S03" "input" "mapping_rows" csv --path "$INPUT" --required-column mining_protein_id

# Validate annotation semantics even when a previous marker exists.
python -c 'import sys; sys.path.insert(0, sys.argv[1]); from eggnog_annotation_utils import read_annotation_rows; read_annotation_rows(sys.argv[2])' "$SCRIPT_DIR" "$EGGNOG_TSV"
if stage_done 03_eggnog_filter && \
    python "$SCRIPT_DIR/step2_candidate_utils.py" check --marker "$(stage_marker 03_eggnog_filter)" \
        --input "$INPUT" --input "$EGGNOG_TSV" --input "$DESCRIPTION_KEYWORDS" --input "$PFAM_KEYWORDS" \
        --output "$OUTPUT" --output "$FINAL_OUTPUT/isozyme_candidates_eggnog_annotated.csv" --mirror "$INPUT" "$OUTPUT"; then
    log_skip "Step2/S03" --device CPU
else
    invalidate_stages 03_eggnog_filter 05_novelty_filter 06_classify
    [[ "$MODE" != custom ]] || invalidate_stages 04_defensefinder
    python "$SCRIPT_DIR/step2_candidate_utils.py" clear-novelty --root "$FINAL_OUTPUT"
    run_stage "Step2/S03" --device CPU -- python -u "$SCRIPT_DIR/S03_eggNOG_filter.py" \
        --input "$INPUT" --eggNOG_tsv "$EGGNOG_TSV" \
        --description_keywords "$DESCRIPTION_KEYWORDS" --pfam_keywords "$PFAM_KEYWORDS" \
        --output "$OUTPUT" --prodigal "$PRODIGAL" --final_output "$FINAL_OUTPUT"
fi

RETAINED="$(candidate_count csv --path "$OUTPUT" --required-column mining_protein_id)"
echo "[Step2/S03] SUMMARY retained_mapping_rows=$RETAINED"
mark_stage_done 03_eggnog_filter
if [[ "$RETAINED" -eq 0 ]]; then
    echo "[Step2/S03] STOP mode=$MODE reason=no_candidates after=eggnog_filter"
fi

