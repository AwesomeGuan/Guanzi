#!/usr/bin/env bash
# ============================================================
# run_esmfold_prediction.sh — predict one PDB per Prodigal protein
#
# Runs esmfold_predict.py in the dedicated esmfold conda environment.
#
# Usage:
# bash Step1_Custom_Mode/run_esmfold_prediction.sh \
# --protein_faa output/prodigal_raw --db_pdb output/custom_pdb \
# --model_path esmfold_v1 --gpus 0,1 [options]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---------- defaults ----------
PROTEIN_FAA=""
DB_PDB=""
GPUS=""
MODEL_PATH=""
MAX_SEQ_LEN=3500
THREAD_PER_GPU=1
STATUS="$PROJECT_ROOT/output/step1_raw/esmfold_status.csv"
INCLUDE_PARENT=""

usage() {
 cat <<EOF
Usage: bash Step1_Custom_Mode/run_esmfold_prediction.sh [OPTIONS]

Required:
 --protein_faa PATH Directory containing .faa files (e.g. prodigal_raw/)
 --db_pdb PATH Output directory for .pdb files
 --model_path PATH Path to ESMFold model (e.g. esmfold_v1, or local path)
 --gpus IDS Comma-separated GPU IDs (use 0 for GPU 0)

Optional:
 --max_seq_len N Max AA length (default: 3500)
 --thread_per_gpu N Workers per GPU (default: 1)
 --mode MODE Log mode: custom or replication (default: unspecified)
 --log_context NAME Optional stage label for shared ESMFold calls (default: ESMFold)
 --status FILE Status CSV (default: output/step1_raw/esmfold_status.csv)
 --include_parent NAME Process only FAA files directly under directories named NAME

Example:
 bash Step1_Custom_Mode/run_esmfold_prediction.sh \\
 --protein_faa output/prodigal_raw --db_pdb output/custom_pdb \\
 --model_path esmfold_v1 --gpus 0,5,6
EOF
 exit 1
}

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
        --log_context) LOG_CONTEXT="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
 --protein_faa) PROTEIN_FAA="$2"; shift 2 ;;
 --db_pdb) DB_PDB="$2"; shift 2 ;;
 --gpus) GPUS="$2"; shift 2 ;;
 --model_path) MODEL_PATH="$2"; shift 2 ;;
 --max_seq_len) MAX_SEQ_LEN="$2"; shift 2 ;;
 --thread_per_gpu) THREAD_PER_GPU="$2"; shift 2 ;;
 --status) STATUS="$2"; shift 2 ;;
 --include_parent) INCLUDE_PARENT="$2"; shift 2 ;;
 -h|--help) usage ;;
 *) echo "Unknown: $1"; usage ;;
 esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }
if [[ "${CONDA_DEFAULT_ENV:-}" != "esmfold" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" esmfold bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


# ---------- validate ----------
if [[ -z "$PROTEIN_FAA" ]] || [[ -z "$DB_PDB" ]] || [[ -z "$MODEL_PATH" ]] || [[ -z "$GPUS" ]]; then
 echo "ERROR: --protein_faa, --db_pdb, --model_path, and --gpus are required"
 usage
fi

if [[ ! -d "$PROTEIN_FAA" ]]; then
 echo "ERROR: --protein_faa directory not found: $PROTEIN_FAA"
 exit 1
fi

python "$PROJECT_ROOT/scripts/preflight.py" gpu --gpus "$GPUS" --torch

mkdir -p "$DB_PDB" "$PROJECT_ROOT/output/step1_raw"

# ---------- activate conda ----------
log() { echo "[${LOG_CONTEXT:-ESMFold}] $*"; }

log "Protein FAA: $PROTEIN_FAA"
log "PDB output: $DB_PDB"
log "Model: $MODEL_PATH"
log "GPUs: $GPUS"


# ---------- run ----------
cd "$PROJECT_ROOT"

export PYTHONUNBUFFERED=1

ESMFOLD_ARGS=( \
 --protein_faa "$PROTEIN_FAA" \
 --db_pdb "$DB_PDB" \
 --model_path "$MODEL_PATH" \
 --gpus "$GPUS" \
 --max_seq_len "$MAX_SEQ_LEN" \
 --thread_per_gpu "$THREAD_PER_GPU" \
 --status "$STATUS"
)
if [[ -n "$INCLUDE_PARENT" ]]; then
 ESMFOLD_ARGS+=(--include_parent "$INCLUDE_PARENT")
fi
run_stage "${LOG_CONTEXT:-ESMFold}" --device GPU --gpus "$GPUS" --processes-per-gpu "$THREAD_PER_GPU" --info "max_seq_len=$MAX_SEQ_LEN" -- python -u "$SCRIPT_DIR/esmfold_predict.py" "${ESMFOLD_ARGS[@]}"

# ---------- post-run ----------
log "prediction complete"
