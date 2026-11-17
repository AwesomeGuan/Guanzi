#!/usr/bin/env bash
# Run the complete Step 4 quality-control chapter.
#
# Stage mapping. S numbers are contiguous and follow the execution order; each
# stage runs its S numbers in ascending order.
#   Stage 1  distance filter      S01
#   Stage 2  clustering           S02-S04 (+ Replication Mining ESMFold)
#   Stage 3  TM-align structural  S05-S08
#   Stage 4  homologous frequency S09 (Replication only)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
PRODIGAL_DIR=""
MODEL_PATH=""
GPUS=""
MAX_SEQ_LEN=3500
THREAD_PER_GPU=1
THREADS=12
WORKERS=4
SHARD_SIZE=128
MAX_GENE_GAP=10

usage() {
    cat <<EOF
Usage: bash Step4_Quality_Control/run_step4.sh --mode custom|replication [OPTIONS]

Required:
  --mode MODE          Select Custom or Replication inputs

Optional:
  --prodigal_dir PATH  Override the mode-specific Prodigal directory
  --model_path PATH    ESMFold model; required for fresh Replication structures
  --gpus IDS           GPU IDs; required for fresh Replication structures
  --max_seq_len N      Maximum ESMFold sequence length (default: 3500)
  --thread_per_gpu N   ESMFold model processes per GPU (default: 1)
  --threads N          CPU threads for CD-HIT, TM-align workers and each DIAMOND worker (default: 12)
  --workers N          Concurrent DIAMOND workers (default: 4)
  --shard_size N       Unique queries per DIAMOND shard (default: 128)
  --max_gene_gap N     Maximum adjacent gene-index difference (default: 10)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --prodigal_dir) PRODIGAL_DIR="$2"; shift 2 ;;
        --model_path) MODEL_PATH="$2"; shift 2 ;;
        --gpus) GPUS="$2"; shift 2 ;;
        --max_seq_len) MAX_SEQ_LEN="$2"; shift 2 ;;
        --thread_per_gpu) THREAD_PER_GPU="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --workers) WORKERS="$2"; shift 2 ;;
        --shard_size) SHARD_SIZE="$2"; shift 2 ;;
        --max_gene_gap) MAX_GENE_GAP="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }
if [[ "$MODE" == "custom" ]]; then
    PRODIGAL_DIR="${PRODIGAL_DIR:-$PROJECT_ROOT/output/prodigal_raw}"
else
    PRODIGAL_DIR="${PRODIGAL_DIR:-$PROJECT_ROOT/output/rep_prodigal_raw}"
fi


source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step4" "mixed" "${GPUS:-}"
echo "[Step4] INFO threads=$THREADS workers=$WORKERS"
source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step4"
source "$PROJECT_ROOT/scripts/stage_state.sh"
mkdir -p "$MARKER_DIR"
RAW="$PROJECT_ROOT/output/step4_raw"
FORMAL="$PROJECT_ROOT/output/step4_output"
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"
REP_MINING_PDB="$RAW/mining_pdb"
REP_MINING_STATUS="$RAW/mining_pdb_status.csv"

stop_if_no_candidates "Step4" "input" "assembled_instance_rows" csv --path "$PROJECT_ROOT/output/step3_output" --pattern "*/assembly.csv" --required-column mining_contig_system_index


if stage_done 01_distance && [[ -s "$RAW/distance_filter/step4_distance_filtered_mapping.csv" ]] \
        && [[ -s "$RAW/distance_filter/step4_distance_filtered_assembly.csv" ]]; then
    log_skip "Step4/S01" --device CPU
    echo "Step 4 distance filter: already complete"
else
    invalidate_stages 01_distance
    bash "$SCRIPT_DIR/run_distance_filter.sh" --mode "$MODE" --max_gene_gap "$MAX_GENE_GAP"
    [[ -s "$RAW/distance_filter/step4_distance_filtered_mapping.csv" ]] || exit 1
    [[ -s "$RAW/distance_filter/step4_distance_filtered_assembly.csv" ]] || exit 1
    mark_stage_done 01_distance
fi

stop_if_no_candidates "Step4" "S01" "distance_filtered_instance_rows" csv --path "$RAW/distance_filter/step4_distance_filtered_assembly.csv" --required-column mining_contig_system_index


if stage_done 02_cluster && [[ -s "$RAW/cdhit_cluster/clustered_mapping.csv" ]] \
        && [[ -s "$RAW/cdhit_cluster/clustered_assembly.csv" ]]; then
    log_skip "Step4/S02-S04" --device CPU
    echo "Step 4 clustering: already complete"
else
    invalidate_stages 02_cluster 02_replication_esmfold
    bash "$SCRIPT_DIR/run_cdhit_cluster.sh" \
        --mode "$MODE" --prodigal_dir "$PRODIGAL_DIR" --threads "$THREADS"
    [[ -s "$RAW/cdhit_cluster/clustered_mapping.csv" ]] || exit 1
    [[ -s "$RAW/cdhit_cluster/clustered_assembly.csv" ]] || exit 1
    mark_stage_done 02_cluster
fi

stop_if_no_candidates "Step4" "S02-S04" "clustered_instance_rows" csv --path "$RAW/cdhit_cluster/clustered_assembly.csv" --required-column mining_contig_system_index


if [[ "$MODE" == "replication" ]]; then
    # Stage 2 Replication branch. The legacy marker name 03_replication_esmfold
    # is still honored so existing runs are not forced to recompute structures.
    if { stage_done 02_replication_esmfold || stage_done 03_replication_esmfold; } \
            && python "$VALIDATOR" esmfold \
            --input "$RAW/cdhit_cluster/qma_faa_clustered" \
            --include-parent mining_assembly --output "$REP_MINING_PDB" \
            --status "$REP_MINING_STATUS" >/dev/null 2>&1; then
        log_skip "Step4/Mining-ESMFold" --device GPU --gpus "$GPUS"
        echo "Step 4 Replication Mining ESMFold: already complete"
    else
        if [[ -z "$MODEL_PATH" || -z "$GPUS" ]]; then
            echo "ERROR: --model_path and --gpus are required for Replication Mining ESMFold" >&2
            exit 2
        fi
        bash "$PROJECT_ROOT/scripts/run_in_env.sh" esmfold \
            python "$PROJECT_ROOT/scripts/preflight.py" gpu --gpus "$GPUS" --torch
        invalidate_stages 02_replication_esmfold 03_replication_esmfold 03_tmalign 04_frequency
        bash "$SCRIPT_DIR/run_mining_esmfold.sh" \
            --model_path "$MODEL_PATH" --gpus "$GPUS" \
            --max_seq_len "$MAX_SEQ_LEN" --thread_per_gpu "$THREAD_PER_GPU"
        python "$VALIDATOR" esmfold \
            --input "$RAW/cdhit_cluster/qma_faa_clustered" \
            --include-parent mining_assembly --output "$REP_MINING_PDB" \
            --status "$REP_MINING_STATUS"
        mark_stage_done 02_replication_esmfold
    fi
fi

TMALIGN_FORMAL_ROWS=0
if [[ -s "$FORMAL/assembly_qualified.csv" && -s "$FORMAL/assembly_disqualified.csv" ]]; then
    TMALIGN_FORMAL_ROWS="$(candidate_count csv \
        --path "$FORMAL/assembly_qualified.csv" \
        --path "$FORMAL/assembly_disqualified.csv" \
        --required-column mining_contig_system_index)"
fi
if stage_done 03_tmalign && [[ -s "$FORMAL/comparison_qualified.csv" ]] \
        && [[ -s "$FORMAL/comparison_disqualified.csv" ]] \
        && [[ "$TMALIGN_FORMAL_ROWS" -gt 0 ]]; then
    log_skip "Step4/S05-S08" --device CPU
    echo "Step 4 TM-align: already complete"
else
    invalidate_stages 03_tmalign 04_frequency
    if [[ "$MODE" == "replication" ]]; then
        python "$VALIDATOR" esmfold \
            --input "$RAW/cdhit_cluster/qma_faa_clustered" \
            --include-parent mining_assembly --output "$REP_MINING_PDB" \
            --status "$REP_MINING_STATUS" >/dev/null 2>&1 || {
                echo "ERROR: Replication Mining structures are incomplete; run run_mining_esmfold.sh" >&2
                exit 1
            }
        MINING_PDB="$REP_MINING_PDB"
    else
        MINING_PDB="$PROJECT_ROOT/output/custom_pdb"
    fi
    bash "$SCRIPT_DIR/run_tmalign.sh" \
        --mode "$MODE" --mining_pdb "$MINING_PDB" --threads "$THREADS"
    [[ -s "$FORMAL/comparison_qualified.csv" ]] || exit 1
    [[ -s "$FORMAL/comparison_disqualified.csv" ]] || exit 1
    TMALIGN_FORMAL_ROWS="$(candidate_count csv \
        --path "$FORMAL/assembly_qualified.csv" \
        --path "$FORMAL/assembly_disqualified.csv" \
        --required-column mining_contig_system_index)"
    if [[ "$TMALIGN_FORMAL_ROWS" -eq 0 ]]; then
        echo "ERROR: TM-align produced zero evaluated instances from non-empty clustered input" >&2
        exit 1
    fi
    mark_stage_done 03_tmalign
fi

for RESULT in qualified disqualified; do
    INSTANCE_ROWS="$(candidate_count csv --path "$FORMAL/assembly_$RESULT.csv" --required-column mining_contig_system_index)"
    echo "[Step4] SUMMARY ${RESULT}_instance_rows=$INSTANCE_ROWS"
done
stop_if_no_candidates "Step4" "S05-S08" "formal_instance_rows" csv --path "$FORMAL/assembly_qualified.csv" --path "$FORMAL/assembly_disqualified.csv" --required-column mining_contig_system_index


if [[ "$MODE" == "replication" ]]; then
    if stage_done 04_frequency && [[ -s "$FORMAL/homologous_frequency/homologous_frequency_all.csv" ]]; then
        log_skip "Step4/S09" --device CPU
        echo "Step 4 homologous frequency: already complete"
    else
        invalidate_stages 04_frequency
        bash "$SCRIPT_DIR/run_replication_homologous_frequency.sh" \
            --threads "$THREADS" --workers "$WORKERS" --shard_size "$SHARD_SIZE"
        [[ -s "$FORMAL/homologous_frequency/homologous_frequency_all.csv" ]] || exit 1
        mark_stage_done 04_frequency
    fi
fi

echo "Step 4 complete for $MODE mode"
