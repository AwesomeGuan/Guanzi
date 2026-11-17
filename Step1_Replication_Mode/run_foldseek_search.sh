#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"

GPUS=""
NUM_THREADS=12
MAX_ROUNDS=5
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"

usage() {
    cat <<EOF
Usage: bash Step1_Replication_Mode/run_foldseek_search.sh --gpus IDS [OPTIONS]

Required:
  --gpus IDS       Comma-separated GPU IDs (use 0 for GPU 0)

Optional:
  --query_pdb PATH Query PDB directory (default: data/query_data/query_pdb)
  --num_threads N  Parallel FoldSeek processes (default: 12)
  --max_rounds N   Maximum search rounds; 0 means unlimited (default: 5)
EOF
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gpus) GPUS="$2"; shift 2 ;;
        --query_pdb) QUERY_PDB="$2"; shift 2 ;;
        --num_threads) NUM_THREADS="$2"; shift 2 ;;
        --max_rounds) MAX_ROUNDS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "foldseek" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" foldseek bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"


if [[ -z "$GPUS" ]]; then
    echo "ERROR: --gpus is required" >&2
    usage
    exit 1
fi

if stage_done 01_foldseek_search && python "$VALIDATOR" foldseek-search \
        --query "$QUERY_PDB" --output "$SEARCH_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S01" --device GPU --gpus "$GPUS"
    echo "FoldSeek already complete"
    exit 0
fi

python "$PROJECT_ROOT/scripts/preflight.py" gpu --gpus "$GPUS" \
    --database "$PROJECT_ROOT/foldseek/${AFDB50_PREFIX}_pad"

mkdir -p "$REP_RAW"
mkdir -p "$MARKER_DIR"
invalidate_stages 01_foldseek_search 02_coverage_filter 03_fident_filter \
    04_extract_afdb50_ids 05_download_isozyme_pdb 06_extract_isozyme_faa \
    07_query_homology_exclusion 08_gtdb_mapping 09_extract_candidate_contigs \
    10_materialize_gtdb_prodigal 11_build_replication_mapping

CUDA_VISIBLE_DEVICES="$GPUS" run_stage "Step1-Replication/S01" --device GPU --gpus "$GPUS" --counted \
    --info "TM-score>=0.7 lDDT>=0.7 parallel_searches=$NUM_THREADS" \
    -- python -u "$SCRIPT_DIR/S01_foldseek_search.py" \
    --query_pdb "$QUERY_PDB" \
    --db_prefix "$AFDB50_PREFIX" \
    --output "$SEARCH_OUT" \
    --num_threads "$NUM_THREADS" \
    --max_rounds "$MAX_ROUNDS"

python "$VALIDATOR" foldseek-search --query "$QUERY_PDB" --output "$SEARCH_OUT"
mark_stage_done 01_foldseek_search
echo "FoldSeek complete"
