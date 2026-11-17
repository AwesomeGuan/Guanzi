#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

GPUS=""
QUERY_PDB=""
NUM_THREADS=12
DOWNLOAD_THREADS=8
MAX_ROUNDS=5
FOLDSEEK_MIN_QCOV=0.70
FOLDSEEK_MIN_TCOV=0.70
FOLDSEEK_FIDENT_BLACKLIST_THRESHOLD=0.30
GTDB_MIN_IDENTITY=90
GTDB_MIN_QCOV=90
GTDB_MIN_TCOV=90
GTDB_EVALUE=1e-5
GTDB_MAX_TARGET_SEQS=0

usage() {
    cat <<EOF
Usage: bash Step1_Replication_Mode/run_replication_mode.sh --gpus IDS [OPTIONS]

Required:
  --gpus IDS       Comma-separated GPU IDs (use 0 for GPU 0)

Optional:
  --query_pdb PATH Query PDB directory (default: data/query_data/query_pdb)
  --num_threads N  Parallel FoldSeek processes (default: 12)
  --download_threads N Concurrent AlphaFold downloads (default: 8)
  --max_rounds N   Maximum FoldSeek search rounds; 0 means unlimited (default: 5)
  --foldseek_min_qcov FLOAT  Minimum FoldSeek query coverage (default: 0.70)
  --foldseek_min_tcov FLOAT  Minimum FoldSeek target coverage (default: 0.70)
  --foldseek_fident_blacklist_threshold FLOAT
                            Target blacklist boundary; values above it are removed (default: 0.30)
  --gtdb_min_identity FLOAT  Minimum GTDB DIAMOND identity percent (default: 90)
  --gtdb_min_qcov FLOAT      Minimum GTDB query coverage percent (default: 90)
  --gtdb_min_tcov FLOAT      Minimum GTDB target coverage percent (default: 90)
  --gtdb_evalue FLOAT        Maximum GTDB DIAMOND E-value (default: 1e-5)
  --gtdb_max_target_seqs N   Maximum GTDB targets; 0 means unlimited (default: 0)
EOF
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gpus) GPUS="$2"; shift 2 ;;
        --query_pdb) QUERY_PDB="$2"; shift 2 ;;
        --num_threads) NUM_THREADS="$2"; shift 2 ;;
        --download_threads) DOWNLOAD_THREADS="$2"; shift 2 ;;
        --max_rounds) MAX_ROUNDS="$2"; shift 2 ;;
        --foldseek_min_qcov) FOLDSEEK_MIN_QCOV="$2"; shift 2 ;;
        --foldseek_min_tcov) FOLDSEEK_MIN_TCOV="$2"; shift 2 ;;
        --foldseek_fident_blacklist_threshold) FOLDSEEK_FIDENT_BLACKLIST_THRESHOLD="$2"; shift 2 ;;
        --gtdb_min_identity) GTDB_MIN_IDENTITY="$2"; shift 2 ;;
        --gtdb_min_qcov) GTDB_MIN_QCOV="$2"; shift 2 ;;
        --gtdb_min_tcov) GTDB_MIN_TCOV="$2"; shift 2 ;;
        --gtdb_evalue) GTDB_EVALUE="$2"; shift 2 ;;
        --gtdb_max_target_seqs) GTDB_MAX_TARGET_SEQS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step1-Replication" "mixed" "${GPUS:-}"


if [[ -z "$GPUS" ]]; then
    echo "ERROR: --gpus is required" >&2
    usage
    exit 1
fi

log() {
    echo "[Replication] === $*"
}

log "Step1 Replication Mode orchestration"

REP_RAW="$PROJECT_ROOT/output/rep_step1_raw"


FOLDSEEK_ARGS=(--gpus "$GPUS" --num_threads "$NUM_THREADS" --max_rounds "$MAX_ROUNDS")
if [[ -n "$QUERY_PDB" ]]; then
    FOLDSEEK_ARGS+=(--query_pdb "$QUERY_PDB")
fi

bash "$SCRIPT_DIR/run_foldseek_search.sh" "${FOLDSEEK_ARGS[@]}"

stop_if_no_candidates "Step1-Replication" "S01" "valid_foldseek_hits" foldseek --path "$REP_RAW/foldseek_search_results"

bash "$SCRIPT_DIR/run_afdb50_filter.sh" \
    --min_qcov "$FOLDSEEK_MIN_QCOV" \
    --min_tcov "$FOLDSEEK_MIN_TCOV" \
    --fident_blacklist_threshold "$FOLDSEEK_FIDENT_BLACKLIST_THRESHOLD"

stop_if_no_candidates "Step1-Replication" "S02-S04" "candidate_ids" ids --path "$REP_RAW/afdb50_isozyme_ids"

bash "$SCRIPT_DIR/run_afdb50_pdb_download.sh" --threads "$DOWNLOAD_THREADS"

stop_if_no_candidates "Step1-Replication" "S05-S06" "successful_pdbs" pdb-status --path "$REP_RAW/afdb50_pdb_download_status.csv" --pdb-root "$REP_RAW/afdb50_isozyme_pdb"

bash "$SCRIPT_DIR/run_query_exclusion.sh"

stop_if_no_candidates "Step1-Replication" "S07" "retained_protein_records" fasta --path "$REP_RAW/diamond_filtered_faa"

bash "$SCRIPT_DIR/run_gtdb_mapping.sh" \
    --threads "$NUM_THREADS" \
    --min_identity "$GTDB_MIN_IDENTITY" \
    --min_qcov "$GTDB_MIN_QCOV" \
    --min_tcov "$GTDB_MIN_TCOV" \
    --evalue "$GTDB_EVALUE" \
    --max_target_seqs "$GTDB_MAX_TARGET_SEQS"

stop_if_no_candidates "Step1-Replication" "S08" "gtdb_mapping_rows" csv --path "$REP_RAW/gtdb_mapping" --pattern "*/*.csv" --required-column gtdb_protein_id

bash "$SCRIPT_DIR/run_gtdb_prodigal_materialization.sh"

stop_if_no_candidates "Step1-Replication" "S09-S10" "materialized_protein_records" fasta --path "$PROJECT_ROOT/output/rep_prodigal_raw" --pattern "**/*.faa"

bash "$SCRIPT_DIR/run_mapping_extraction.sh"

log "Replication Mode complete"
