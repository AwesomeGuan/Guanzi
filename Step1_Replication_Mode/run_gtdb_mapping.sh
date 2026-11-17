#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
THREADS=12
MIN_IDENTITY=90
MIN_QCOV=90
MIN_TCOV=90
EVALUE=1e-5
MAX_TARGET_SEQS=0
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --min_identity) MIN_IDENTITY="$2"; shift 2 ;;
        --min_qcov) MIN_QCOV="$2"; shift 2 ;;
        --min_tcov) MIN_TCOV="$2"; shift 2 ;;
        --evalue) EVALUE="$2"; shift 2 ;;
        --max_target_seqs) MAX_TARGET_SEQS="$2"; shift 2 ;;
        -h|--help)
            cat <<EOF
Usage: bash Step1_Replication_Mode/run_gtdb_mapping.sh [OPTIONS]
  --threads N             DIAMOND threads (default: 12)
  --min_identity FLOAT    Minimum identity percent (default: 90)
  --min_qcov FLOAT        Minimum query coverage percent (default: 90)
  --min_tcov FLOAT        Minimum target coverage percent (default: 90)
  --evalue FLOAT          Maximum E-value (default: 1e-5)
  --max_target_seqs N     Maximum targets; 0 means unlimited (default: 0)
EOF
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"


bash "$SCRIPT_DIR/run_prepare_gtdb_diamond_db.sh" --threads "$THREADS"

HAD_MARKER=0
if stage_done 08_gtdb_mapping; then
    HAD_MARKER=1
fi
if [[ "$HAD_MARKER" -eq 1 ]] && \
        python "$VALIDATOR" gtdb --input "$DIAMOND_OUT" --output "$GTDB_MAPPING_OUT" \
        --min-identity "$MIN_IDENTITY" --min-qcov "$MIN_QCOV" \
        --min-tcov "$MIN_TCOV" --evalue "$EVALUE" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S08" --device CPU
    echo "GTDB mapping already complete"
else
    invalidate_stages 08_gtdb_mapping
    GTDB_ARGS=(
        --query_dir "$DIAMOND_OUT" --diamond_db_dir "$GTDB_DIAMOND_ROOT"
        --output "$GTDB_MAPPING_OUT"
        --threads "$THREADS"
        --min_identity "$MIN_IDENTITY" --min_qcov "$MIN_QCOV"
        --min_tcov "$MIN_TCOV" --evalue "$EVALUE"
        --max_target_seqs "$MAX_TARGET_SEQS"
    )
    if [[ "$HAD_MARKER" -eq 0 ]]; then
        GTDB_ARGS+=(--reset)
    fi
    run_stage "Step1-Replication/S08" --device CPU --threads "$THREADS" --counted --info "identity>=$MIN_IDENTITY% qcov>=$MIN_QCOV% tcov>=$MIN_TCOV% evalue<=$EVALUE max_target_seqs=$MAX_TARGET_SEQS" -- python -u "$SCRIPT_DIR/S08_gtdb_mapping.py" "${GTDB_ARGS[@]}"
    python "$VALIDATOR" gtdb --input "$DIAMOND_OUT" --output "$GTDB_MAPPING_OUT" \
        --min-identity "$MIN_IDENTITY" --min-qcov "$MIN_QCOV" \
        --min-tcov "$MIN_TCOV" --evalue "$EVALUE"
    mark_stage_done 08_gtdb_mapping
fi
