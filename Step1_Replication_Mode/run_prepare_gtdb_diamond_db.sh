#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
THREADS=12
FORCE=0

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help)
            echo "Usage: bash Step1_Replication_Mode/run_prepare_gtdb_diamond_db.sh [--threads N] [--force]"
            echo "  --threads N  DIAMOND threads (default: 12)"
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"



args=(--protein_faa_reps "$PROTEIN_FAA_REPS" --output "$GTDB_DIAMOND_ROOT" --threads "$THREADS")
[[ "$FORCE" -eq 1 ]] && args+=(--force)
run_stage "Step1-Replication/GTDB-DB" --device CPU --threads "$THREADS" -- python -u "$SCRIPT_DIR/prepare_gtdb_diamond_db.py" "${args[@]}"
