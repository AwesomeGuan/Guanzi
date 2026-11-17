#!/usr/bin/env bash
# Run the complete Custom Mode workflow with stage-level resume markers.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/custom_paths.sh"

GENOME_DIR=""
MODEL_PATH=""
GPUS=""
NUM_THREADS=12
MAX_SEQ_LEN=3500
THREAD_PER_GPU=1
MAX_ROUNDS=5
MIN_QCOV=0.70
MIN_TCOV=0.70
FIDENT_BLACKLIST_THRESHOLD=0.30
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"

usage() {
    cat <<EOF
Usage: bash Step1_Custom_Mode/run_custom_mode.sh [OPTIONS]

Required for a fresh run:
  --genome_dir PATH       Directory containing MAG/genome FASTA files
  --model_path PATH       ESMFold model path
  --gpus IDS              Comma-separated GPU IDs (use 0 for GPU 0)

Optional:
  --query_pdb PATH        Query PDB directory (default: data/query_data/query_pdb)
  --db_prefix NAME        FoldSeek database prefix (default: custom)
  --num_threads N         Parallel FoldSeek searches (default: 12)
  --max_seq_len N         Maximum ESMFold sequence length (default: 3500)
  --thread_per_gpu N      ESMFold workers per GPU (default: 1)
  --max_rounds N          FoldSeek search rounds; 0 means unlimited (default: 5)
  --min_qcov FLOAT        Minimum FoldSeek query coverage (default: 0.70)
  --min_tcov FLOAT        Minimum FoldSeek target coverage (default: 0.70)
  --fident_blacklist_threshold FLOAT
                          Remove a target if any coverage-passing hit is above
                          this FoldSeek identity (default: 0.30)
EOF
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome_dir) GENOME_DIR="$2"; shift 2 ;;
        --model_path) MODEL_PATH="$2"; shift 2 ;;
        --query_pdb) QUERY_PDB="$2"; shift 2 ;;
        --db_prefix) DB_PREFIX="$2"; shift 2 ;;
        --gpus) GPUS="$2"; shift 2 ;;
        --num_threads) NUM_THREADS="$2"; shift 2 ;;
        --max_seq_len) MAX_SEQ_LEN="$2"; shift 2 ;;
        --thread_per_gpu) THREAD_PER_GPU="$2"; shift 2 ;;
        --max_rounds) MAX_ROUNDS="$2"; shift 2 ;;
        --min_qcov) MIN_QCOV="$2"; shift 2 ;;
        --min_tcov) MIN_TCOV="$2"; shift 2 ;;
        --fident_blacklist_threshold) FIDENT_BLACKLIST_THRESHOLD="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=custom
source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step1-Custom" "mixed" "${GPUS:-}"



PRODIGAL_REUSABLE=0
if stage_done 01_prodigal && python "$VALIDATOR" prodigal --output "$PRODIGAL_OUT" >/dev/null 2>&1; then
    PRODIGAL_REUSABLE=1
fi
ESMFOLD_REUSABLE=0
if stage_done 02_esmfold && \
        python "$VALIDATOR" esmfold --output "$DB_PDB_OUT" \
        --status "$CUSTOM_RAW/esmfold_status.csv" \
        --input "$PRODIGAL_OUT" >/dev/null 2>&1; then
    ESMFOLD_REUSABLE=1
fi

if [[ "$PRODIGAL_REUSABLE" -eq 0 && -z "$GENOME_DIR" ]]; then
    echo "ERROR: --genome_dir is required until Prodigal is complete" >&2
    usage
    exit 1
fi
if [[ "$ESMFOLD_REUSABLE" -eq 0 && -z "$MODEL_PATH" ]]; then
    echo "ERROR: --model_path is required until ESMFold is complete" >&2
    usage
    exit 1
fi
if [[ -z "$GPUS" ]]; then
    echo "ERROR: --gpus is required" >&2
    usage
    exit 1
fi

log() { echo "[Custom] $*"; }
mkdir -p "$CUSTOM_RAW" "$MARKER_DIR"

if [[ "$PRODIGAL_REUSABLE" -eq 1 ]]; then
    stop_if_no_candidates "Step1-Custom" "Prodigal" "protein_records" fasta --path "$PRODIGAL_OUT" --pattern "**/*.faa"
    log_skip "Step1-Custom/Prodigal" --device CPU
    log "[1/7] Prodigal already complete"
else
    log "[1/7] Predicting CDS and proteins with Prodigal"
    invalidate_stages 01_prodigal 02_esmfold 03_foldseek_db 04_foldseek_search \
        05_coverage_filter 06_fident_filter 07_mapping
    bash "$SCRIPT_DIR/run_prodigal.sh" \
        --genome_dir "$GENOME_DIR" --outdir "$PRODIGAL_OUT"
    stop_if_no_candidates "Step1-Custom" "Prodigal" "protein_records" fasta --path "$PRODIGAL_OUT" --pattern "**/*.faa"
    python "$VALIDATOR" prodigal --output "$PRODIGAL_OUT"
    mark_stage_done 01_prodigal
fi


if [[ "$ESMFOLD_REUSABLE" -eq 1 ]]; then
    log_skip "Step1-Custom/ESMFold" --device GPU --gpus "$GPUS"
    log "[2/7] ESMFold already complete"
else
    log "[2/7] Predicting one PDB per Prodigal protein"
    bash "$PROJECT_ROOT/scripts/run_in_env.sh" esmfold python \
        "$PROJECT_ROOT/scripts/preflight.py" gpu --gpus "$GPUS" --torch
    invalidate_stages 02_esmfold 03_foldseek_db 04_foldseek_search \
        05_coverage_filter 06_fident_filter 07_mapping
    bash "$SCRIPT_DIR/run_esmfold_prediction.sh" --mode "$MODE" \
        --log_context Step1-Custom/ESMFold \
        --protein_faa "$PRODIGAL_OUT" --db_pdb "$DB_PDB_OUT" \
        --model_path "$MODEL_PATH" \
        --gpus "$GPUS" --max_seq_len "$MAX_SEQ_LEN" \
        --thread_per_gpu "$THREAD_PER_GPU"
    python "$VALIDATOR" esmfold --output "$DB_PDB_OUT" \
        --status "$CUSTOM_RAW/esmfold_status.csv" --input "$PRODIGAL_OUT"
    mark_stage_done 02_esmfold
fi

stop_if_no_candidates "Step1-Custom" "ESMFold" "available_pdbs" pdb --path "$DB_PDB_OUT"


if stage_done 03_foldseek_db && python "$VALIDATOR" foldseek-db \
        --prefix "$FOLDSEEK_DIR/$DB_PREFIX" >/dev/null 2>&1; then
    log_skip "Step1-Custom/FoldSeek-DB" --device CPU
    log "[3/7] FoldSeek database already complete"
else
    log "[3/7] Building FoldSeek database"
    invalidate_stages 03_foldseek_db 04_foldseek_search 05_coverage_filter \
        06_fident_filter 07_mapping
    bash "$SCRIPT_DIR/run_foldseek_db.sh" \
        --db_pdb "$DB_PDB_OUT" --db_prefix "$DB_PREFIX" \
        --threads "$NUM_THREADS"
    python "$VALIDATOR" foldseek-db --prefix "$FOLDSEEK_DIR/$DB_PREFIX"
    mark_stage_done 03_foldseek_db
fi

if stage_done 04_foldseek_search && python "$VALIDATOR" foldseek-search \
        --query "$QUERY_PDB" --output "$SEARCH_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Custom/S01" --device GPU --gpus "$GPUS"
    log "[4/7] FoldSeek search already complete"
else
    log "[4/7] Searching preset Query PDBs"
    python "$PROJECT_ROOT/scripts/preflight.py" gpu --gpus "$GPUS" \
        --database "$PROJECT_ROOT/foldseek/${DB_PREFIX}_pad"
    invalidate_stages 04_foldseek_search 05_coverage_filter 06_fident_filter 07_mapping
    bash "$SCRIPT_DIR/run_foldseek_search.sh" \
        --query_pdb "$QUERY_PDB" --db_prefix "$DB_PREFIX" \
        --gpus "$GPUS" --num_threads "$NUM_THREADS" \
        --max_rounds "$MAX_ROUNDS" --output "$SEARCH_OUT"
    python "$VALIDATOR" foldseek-search --query "$QUERY_PDB" --output "$SEARCH_OUT"
    mark_stage_done 04_foldseek_search
fi

stop_if_no_candidates "Step1-Custom" "S01" "valid_foldseek_hits" foldseek --path "$SEARCH_OUT"


if stage_done 05_coverage_filter && \
        python "$VALIDATOR" coverage --input "$SEARCH_OUT" --output "$COVERAGE_OUT" \
        --min-qcov "$MIN_QCOV" --min-tcov "$MIN_TCOV" >/dev/null 2>&1; then
    log_skip "Step1-Custom/S02" --device CPU
    log "[5/7] coverage filter already complete"
else
    log "[5/7] Applying coverage filter"
    invalidate_stages 05_coverage_filter 06_fident_filter 07_mapping
    run_stage "Step1-Custom/S02" --device CPU --info "qcov>=${MIN_QCOV:-0.70} tcov>=${MIN_TCOV:-0.70}" -- python -u "$SCRIPT_DIR/S02_coverage_filter.py" \
        --input "$SEARCH_OUT" --output "$COVERAGE_OUT" \
        --min_qcov "$MIN_QCOV" --min_tcov "$MIN_TCOV"
    python "$VALIDATOR" coverage --input "$SEARCH_OUT" --output "$COVERAGE_OUT" \
        --min-qcov "$MIN_QCOV" --min-tcov "$MIN_TCOV"
    mark_stage_done 05_coverage_filter
fi

stop_if_no_candidates "Step1-Custom" "S02" "coverage_rows" csv --path "$COVERAGE_OUT" --required-column subject_id


if stage_done 06_fident_filter && \
        python "$VALIDATOR" fident --input "$COVERAGE_OUT" --output "$FIDENT_OUT" \
        --threshold "$FIDENT_BLACKLIST_THRESHOLD" >/dev/null 2>&1; then
    log_skip "Step1-Custom/S03" --device CPU
    log "[6/7] fident filter already complete"
else
    log "[6/7] Applying target-level fident blacklist and target deduplication"
    invalidate_stages 06_fident_filter 07_mapping
    run_stage "Step1-Custom/S03" --device CPU --info "target_blacklist: any fident>${FIDENT_BLACKLIST_THRESHOLD:-0.30}; then max-TM-score deduplication" -- python -u "$SCRIPT_DIR/S03_fident_filter.py" \
        --input "$COVERAGE_OUT" --output "$FIDENT_OUT" \
        --fident_blacklist_threshold "$FIDENT_BLACKLIST_THRESHOLD"
    python "$VALIDATOR" fident --input "$COVERAGE_OUT" --output "$FIDENT_OUT" \
        --threshold "$FIDENT_BLACKLIST_THRESHOLD"
    mark_stage_done 06_fident_filter
fi

stop_if_no_candidates "Step1-Custom" "S03" "candidate_targets" csv --path "$FIDENT_OUT" --required-column subject_id


if stage_done 07_mapping && python "$VALIDATOR" custom-mapping \
        --input "$FIDENT_OUT" --output "$MAPPING_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Custom/S04" --device CPU
    log "[7/7] mapping already complete"
else
    log "[7/7] Writing Query-Mining mappings"
    invalidate_stages 07_mapping
    run_stage "Step1-Custom/S04" --device CPU -- python -u "$SCRIPT_DIR/S04_custom_mapping.py" \
        --isozyme_dir "$FIDENT_OUT" --output_dir "$MAPPING_OUT"
    python "$VALIDATOR" custom-mapping --input "$FIDENT_OUT" --output "$MAPPING_OUT"
    mark_stage_done 07_mapping
fi

log "Custom Mode complete: $MAPPING_OUT"
