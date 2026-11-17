#!/usr/bin/env bash
# Classify novel mappings and supplement Mining and Query location context.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
INPUT="$PROJECT_ROOT/output/isozyme_candidates_novel.csv"
STYLE_FILE="$PROJECT_ROOT/data/defense_system_style.csv"
OUTPUT="$PROJECT_ROOT/output/step2_output"
PRODIGAL_RAW=""
RULES="$PROJECT_ROOT/data/DefenseFinder_rules"
QUERY_LOCATION="$PROJECT_ROOT/data/query_data/query_location"
ORIGINAL_ARGS=("$@")

usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_classify_mining_systems.sh --mode custom|replication [OPTIONS]

Required:
  --mode MODE             Select the Prodigal input layout

Optional:
  --input FILE            Novel candidate table
  --style_file FILE       Defense-system style table
  --output DIR            Classified Step 2 output directory
  --prodigal_raw DIR      Override the mode-selected Prodigal directory
  --rules DIR             DefenseFinder rule directory
  --query_location DIR    Normalized Query-location directory
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --style_file) STYLE_FILE="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --prodigal_raw) PRODIGAL_RAW="$2"; shift 2 ;;
        --rules) RULES="$2"; shift 2 ;;
        --query_location) QUERY_LOCATION="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }
if [[ "$MODE" == "custom" ]]; then
    PRODIGAL_RAW="${PRODIGAL_RAW:-$PROJECT_ROOT/output/prodigal_raw}"
else
    PRODIGAL_RAW="${PRODIGAL_RAW:-$PROJECT_ROOT/output/rep_prodigal_raw}"
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

classification_complete() {
    local group
    for group in single multi complex cas; do
        [[ -d "$OUTPUT/step2_mapping_$group" ]] || return 1
    done
}

log_script_start "Step2/S06-S08" "CPU"
stop_if_no_candidates "Step2/S06-S08" "input" "novel_mapping_rows" csv --path "$INPUT" --required-column mining_protein_id

if stage_done 06_classify && classification_complete && \
    python "$SCRIPT_DIR/step2_candidate_utils.py" check --marker "$(stage_marker 06_classify)" \
        --input "$INPUT" --input "$STYLE_FILE" --input "$RULES" --input "$QUERY_LOCATION" \
        --output "$OUTPUT/step2_mapping_single" --output "$OUTPUT/step2_mapping_multi" \
        --output "$OUTPUT/step2_mapping_complex" --output "$OUTPUT/step2_mapping_cas"; then
    log_skip "Step2/S06-S08" --device CPU
else
    invalidate_stages 06_classify
    run_stage "Step2/S06" --device CPU -- python -u "$SCRIPT_DIR/S06_classify_defensefinder.py" \
        --input "$INPUT" --style_file "$STYLE_FILE" --output "$OUTPUT"

    for GROUP in single multi complex cas; do
        GROUP_DIR="$OUTPUT/step2_mapping_$GROUP"
        if [[ ! -d "$GROUP_DIR" ]] || ! compgen -G "$GROUP_DIR/*.csv" >/dev/null; then
            continue
        fi
        run_stage "Step2/S07-$GROUP" --device CPU -- python -u "$SCRIPT_DIR/S07_add_genome_id.py" \
            --input "$GROUP_DIR" --prodigal "$PRODIGAL_RAW"
        run_stage "Step2/S08-$GROUP" --device CPU -- python -u "$SCRIPT_DIR/S08_add_query_location.py" \
            --input "$GROUP_DIR" --query_location "$QUERY_LOCATION" --rules_dir "$RULES"
    done
    classification_complete || { echo "ERROR: classification output is incomplete: $OUTPUT" >&2; exit 1; }
    mark_stage_done 06_classify
fi

for GROUP in single multi complex cas; do
    ROWS="$(candidate_count csv --path "$OUTPUT/step2_mapping_$GROUP" --required-column mining_protein_id)"
    echo "[Step2/S06-S08] SUMMARY ${GROUP}_mapping_rows=$ROWS"
done
