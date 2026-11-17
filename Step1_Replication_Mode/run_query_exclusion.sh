#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --query_location) QUERY_LOCATION="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash Step1_Replication_Mode/run_query_exclusion.sh [--query_location PATH]"
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"


if stage_done 07_query_homology_exclusion && python "$VALIDATOR" query-exclusion \
        --input "$AFDB50_FAA_OUT" --output "$DIAMOND_OUT" \
        --markers "$REP_RAW/.task_markers/diamond" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S07" --device CPU
    echo "DIAMOND query exclusion already complete"
    exit 0
fi

mkdir -p "$REP_RAW"
mkdir -p "$MARKER_DIR"
invalidate_stages 07_query_homology_exclusion


run_stage "Step1-Replication/S07" --device CPU --counted --info "identity>=70% qcov>=70% tcov>=70%" -- python -u "$SCRIPT_DIR/S07_query_homology_exclusion.py" \
    --query_location "$QUERY_LOCATION" \
    --subject_dir "$AFDB50_FAA_OUT" \
    --database_dir "$QUERY_DIAMOND_DB" \
    --output "$DIAMOND_OUT" \
    --marker_dir "$REP_RAW/.task_markers/diamond"

python "$VALIDATOR" query-exclusion --input "$AFDB50_FAA_OUT" \
    --output "$DIAMOND_OUT" --markers "$REP_RAW/.task_markers/diamond"
mark_stage_done 07_query_homology_exclusion
echo "DIAMOND query exclusion complete"
