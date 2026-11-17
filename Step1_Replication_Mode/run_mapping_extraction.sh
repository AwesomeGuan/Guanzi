#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            echo "Usage: bash Step1_Replication_Mode/run_mapping_extraction.sh"
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"


mkdir -p "$REP_RAW"
mkdir -p "$MARKER_DIR"


if stage_done 11_build_replication_mapping && python "$VALIDATOR" final \
        --gtdb "$GTDB_MAPPING_OUT" --foldseek "$FIDENT_OUT" \
        --prodigal "$PRODIGAL_OUT" --output "$MAPPING_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S11" --device CPU
    echo "Replication mapping already complete"
else
    invalidate_stages 11_build_replication_mapping
    run_stage "Step1-Replication/S11" --device CPU -- python -u "$SCRIPT_DIR/S11_build_replication_mapping.py" \
        --gtdb_mapping "$GTDB_MAPPING_OUT" --foldseek_results "$FIDENT_OUT" \
        --prodigal "$PRODIGAL_OUT" --output "$MAPPING_OUT"
    python "$VALIDATOR" final --gtdb "$GTDB_MAPPING_OUT" \
        --foldseek "$FIDENT_OUT" --prodigal "$PRODIGAL_OUT" --output "$MAPPING_OUT"
    mark_stage_done 11_build_replication_mapping
fi

echo "Mapping extraction complete"
