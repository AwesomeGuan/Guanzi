#!/usr/bin/env bash
# Complete CPU-only Step2 chapter, including local eggNOG annotation.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
THREADS=12
ORIGINAL_ARGS=("$@")
usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_step2.sh --mode custom|replication [--threads N]
Required:
  --mode MODE          Select Custom or Replication inputs and DefenseFinder source
Optional:
  --threads N          eggNOG/DIAMOND CPUs (default: 12)
EOF
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == custom || "$MODE" == replication ]] || { usage; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be positive" >&2; exit 2; }
if [[ "${CONDA_DEFAULT_ENV:-}" != defense_system ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
source "$PROJECT_ROOT/scripts/stage_state.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
log_script_start Step2 CPU
echo "[Step2] INFO eggnog_threads=$THREADS"
INPUT="$PROJECT_ROOT/output/step1_mapping"
[[ "$MODE" != replication ]] || INPUT="$PROJECT_ROOT/output/rep_step1_mapping"
ROWS="$(candidate_count csv --path "$INPUT" --required-column mining_protein_id)"
echo "[Step2] SUMMARY mapping_rows=$ROWS"
if [[ "$ROWS" -eq 0 ]]; then
    invalidate_stages 01_prepare_eggnog 02_eggnog_annotation 03_eggnog_filter 04_defensefinder 05_novelty_filter 06_classify
    python "$SCRIPT_DIR/step2_candidate_utils.py" clear-eggnog --root "$PROJECT_ROOT/output"
    echo "[Step2] STOP mode=$MODE reason=no_candidates after=input"
    exit 0
fi
bash "$SCRIPT_DIR/run_prepare_eggnog.sh" --mode "$MODE"
ROWS="$(candidate_count fasta --path "$PROJECT_ROOT/output/step2_raw/eggnog_upload.faa")"
echo "[Step2] SUMMARY upload_protein_records=$ROWS"
if [[ "$ROWS" -eq 0 ]]; then
    invalidate_stages 02_eggnog_annotation 03_eggnog_filter 05_novelty_filter 06_classify
    python "$SCRIPT_DIR/step2_candidate_utils.py" clear-eggnog --root "$PROJECT_ROOT/output"
    echo "[Step2] STOP mode=$MODE reason=no_candidates after=prepare"
    exit 0
fi
bash "$SCRIPT_DIR/run_eggnog_annotation.sh" --threads "$THREADS"
bash "$SCRIPT_DIR/run_eggnog_filter.sh" --mode "$MODE"
stop_if_no_candidates Step2 eggnog_filter retained_mapping_rows csv --path "$PROJECT_ROOT/output/step2_raw/eggnog_filtered" --required-column mining_protein_id
if [[ "$MODE" == custom ]]; then
    bash "$SCRIPT_DIR/run_defensefinder.sh" --mode custom
else
    bash "$SCRIPT_DIR/run_replication_defensefinder_import.sh" --mode replication
fi
bash "$SCRIPT_DIR/run_defensefinder_filter.sh" --mode "$MODE"
stop_if_no_candidates Step2 novelty_filter novel_mapping_rows csv --path "$PROJECT_ROOT/output/isozyme_candidates_novel.csv" --required-column mining_protein_id
bash "$SCRIPT_DIR/run_classify_mining_systems.sh" --mode "$MODE"
echo "[Step2] SUMMARY output=$PROJECT_ROOT/output/step2_output"
