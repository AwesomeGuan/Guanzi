#!/usr/bin/env bash
# Stage 2 Replication branch: predict structures only for Mining components
# retained after CD-HIT. Runs after run_cdhit_cluster.sh and before TM-align.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT="$PROJECT_ROOT/output/step4_raw/cdhit_cluster/qma_faa_clustered"
OUTPUT="$PROJECT_ROOT/output/step4_raw/mining_pdb"
STATUS="$PROJECT_ROOT/output/step4_raw/mining_pdb_status.csv"
MODEL_PATH=""
GPUS=""
MAX_SEQ_LEN=3500
THREAD_PER_GPU=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        --model_path) MODEL_PATH="$2"; shift 2 ;;
        --gpus) GPUS="$2"; shift 2 ;;
        --max_seq_len) MAX_SEQ_LEN="$2"; shift 2 ;;
        --thread_per_gpu) THREAD_PER_GPU="$2"; shift 2 ;;
        -h|--help)
            cat <<EOF
Usage: bash Step4_Quality_Control/run_mining_esmfold.sh --model_path PATH --gpus IDS [OPTIONS]
  --max_seq_len N      Maximum sequence length (default: 3500)
  --thread_per_gpu N   ESMFold model processes per GPU (default: 1)
EOF
            exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 2 ;;
    esac
done
if [[ -z "$MODEL_PATH" || -z "$GPUS" ]]; then
    echo "ERROR: --model_path and --gpus are required" >&2
    exit 2
fi

MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"
[[ -d "$INPUT" ]] || { echo "ERROR: clustered QMA FAA tree not found: $INPUT" >&2; exit 1; }
bash "$PROJECT_ROOT/Step1_Custom_Mode/run_esmfold_prediction.sh" --mode replication --log_context Step4/Mining-ESMFold \
    --protein_faa "$INPUT" --db_pdb "$OUTPUT" \
    --model_path "$MODEL_PATH" --gpus "$GPUS" \
    --max_seq_len "$MAX_SEQ_LEN" --thread_per_gpu "$THREAD_PER_GPU" \
    --status "$STATUS" --include_parent mining_assembly
