#!/usr/bin/env bash
# Run or validate the local eggNOG stage; native files and canonical TSV are separate.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
THREADS=12
ANNOTATIONS=""
ORIGINAL_ARGS=("$@")
usage() {
    cat <<'EOF'
Usage: bash Step2_Functional_Filter/run_eggnog_annotation.sh [OPTIONS]
Optional:
  --threads N          eggNOG/DIAMOND CPUs (default: 12)
  --annotations FILE   Import a completed native v2 result for the current FASTA
                       Caller must confirm identical input and documented settings
EOF
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --annotations) ANNOTATIONS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be positive" >&2; exit 2; }
if [[ "${CONDA_DEFAULT_ENV:-}" != eggnog ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" eggnog bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_state.sh"
source "$PROJECT_ROOT/scripts/stage_logging.sh"
[[ -f "$PROJECT_ROOT/output/.active_mode" ]] || { echo "ERROR: run Step2 preparation first to select the active mode" >&2; exit 1; }
MODE="$(tr -d '[:space:]' < "$PROJECT_ROOT/output/.active_mode")"
[[ "$MODE" == custom || "$MODE" == replication ]] || { echo "ERROR: invalid output/.active_mode" >&2; exit 1; }
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step2"
mkdir -p "$MARKER_DIR"
log_script_start "Step2/S02" CPU
echo "[Step2/S02] INFO threads=$THREADS"
stop_if_no_candidates "Step2/S02" input upload_protein_records fasta --path "$PROJECT_ROOT/output/step2_raw/eggnog_upload.faa"
TSV="$PROJECT_ROOT/output/step2_raw/eggnog_results.tsv"
NATIVE="$PROJECT_ROOT/output/step2_raw/eggnog_results/query.emapper.annotations"
if [[ -z "$ANNOTATIONS" ]] && stage_done 02_eggnog_annotation && [[ -s "$NATIVE" && -s "$TSV" ]] && \
    python -c 'import sys; sys.path.insert(0,sys.argv.pop(1)); from S02_annotate_eggNOG import outputs_match; sys.exit(0 if outputs_match(*sys.argv[1:]) else 1)' \
        "$SCRIPT_DIR" "$NATIVE" "$TSV" "$PROJECT_ROOT/output/step2_raw/eggnog_upload.faa"; then
    log_skip "Step2/S02" --device CPU
    exit 0
fi
digest() {
    python -c 'import hashlib,sys; from pathlib import Path; p=Path(sys.argv[1]); print(hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "missing")' "$TSV"
}
BEFORE="$(digest)"
invalidate_stages 02_eggnog_annotation
COMMAND=(python -u "$SCRIPT_DIR/S02_annotate_eggNOG.py" --threads "$THREADS")
[[ -z "$ANNOTATIONS" ]] || COMMAND+=(--annotations "$ANNOTATIONS")
if ! run_stage "Step2/S02" --device CPU -- "${COMMAND[@]}"; then
    invalidate_stages 03_eggnog_filter 05_novelty_filter 06_classify
    [[ "$MODE" != custom ]] || invalidate_stages 04_defensefinder
    exit 1
fi
if [[ "$BEFORE" != "$(digest)" ]]; then
    invalidate_stages 03_eggnog_filter 05_novelty_filter 06_classify
    [[ "$MODE" != custom ]] || invalidate_stages 04_defensefinder
fi
mark_stage_done 02_eggnog_annotation
