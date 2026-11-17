#!/usr/bin/env bash
# ============================================================
# FoldSeek structural search - standalone S01 runner
#
# Searches ALL query protein PDBs against user's FoldSeek DB.
# Output goes to --output (per-system subdirectories).
# This is the GPU-intensive step; run once per DB.
#
# Usage:
# # Normal mode:
# bash run_foldseek_search.sh --query_pdb data/query_data/query_pdb --db_prefix mydb --gpus 0 --num_threads 12 --max_rounds 5 --output output/step1_raw/foldseek_search_results
#
# # Replication Mode uses Step1_Replication_Mode/run_foldseek_search.sh.
# bash run_foldseek_search.sh --query_pdb data/query_data/query_pdb --db_prefix afdb50 --gpus 0 --num_threads 12 --max_rounds 5 --output output/step1_raw/foldseek_search_afdb50_gtdb_raw
#
# Script auto-appends _pad: foldseek/{db_prefix}_pad
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DB_PREFIX=""
QUERY_PDB="$PROJECT_ROOT/data/query_data/query_pdb"
GPUS=""
NUM_THREADS="12"
MAX_ROUNDS="5"
OUTPUT="$PROJECT_ROOT/output/step1_raw/foldseek_search_results"

usage() {
 cat <<EOF
Usage: bash run_foldseek_search.sh [OPTIONS]

Options:
 --db_prefix NAME FoldSeek DB prefix (required; uses foldseek/{db_prefix}_pad)
 --query_pdb PATH Query PDB dir (default: data/query_data/query_pdb)
 --gpus IDS Comma-separated GPU IDs for CUDA_VISIBLE_DEVICES (required)
 --num_threads N Number of parallel foldseek processes (default: 12)
 --max_rounds N Max search rounds; 0 means unlimited (default: 5)
 --output PATH Raw search output dir (default: output/step1_raw/foldseek_search_results)
EOF
 exit 1
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
 --db_prefix) DB_PREFIX="$2"; shift 2 ;;
 --query_pdb) QUERY_PDB="$2"; shift 2 ;;
 --gpus) GPUS="$2"; shift 2 ;;
 --num_threads) NUM_THREADS="$2"; shift 2 ;;
 --max_rounds) MAX_ROUNDS="$2"; shift 2 ;;
 --output) OUTPUT="$2"; shift 2 ;;
 -h|--help) usage ;;
 *) echo "Unknown: $1"; usage ;;
 esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "foldseek" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" foldseek bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=custom
source "$PROJECT_ROOT/scripts/stage_logging.sh"


if [[ -z "$DB_PREFIX" ]] || [[ -z "$GPUS" ]]; then
 echo "ERROR: --db_prefix and --gpus are required"
 usage
fi

python "$PROJECT_ROOT/scripts/preflight.py" gpu --gpus "$GPUS" \
 --database "$PROJECT_ROOT/foldseek/${DB_PREFIX}_pad"

log() { echo "[FoldSeek] $*"; }

log "=== FoldSeek Search (S01) ==="
log " Query PDB: $QUERY_PDB"
log " DB prefix: $DB_PREFIX"
log " GPU: $GPUS"
log " Num threads: $NUM_THREADS"
log " Output: $OUTPUT"
if [[ "$MAX_ROUNDS" != "0" ]]; then
 log " Max rounds: $MAX_ROUNDS"
fi

S01_ARGS=(
 --query_pdb "$QUERY_PDB"
 --db_prefix "$DB_PREFIX"
 --output "$OUTPUT"
 --num_threads "$NUM_THREADS"
)
if [[ "$MAX_ROUNDS" != "0" ]]; then
 S01_ARGS+=(--max_rounds "$MAX_ROUNDS")
fi

# ---------- activate conda ----------

CUDA_VISIBLE_DEVICES="$GPUS" run_stage "Step1-Custom/S01" --device GPU --gpus "$GPUS" --counted \
 --info "TM-score>=0.7 lDDT>=0.7 parallel_searches=$NUM_THREADS" \
 -- python -u "$SCRIPT_DIR/S01_foldseek_search.py" "${S01_ARGS[@]}"

log "=== FoldSeek search done ==="
log "Raw output: $OUTPUT/"
