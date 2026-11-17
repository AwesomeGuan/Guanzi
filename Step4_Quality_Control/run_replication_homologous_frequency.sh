#!/usr/bin/env bash
# Stage 4: GTDB homologous-frequency filter (S09, Replication only)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
STEP4_OUTPUT="$PROJECT_ROOT/output/step4_output"
WORK_DIR="$PROJECT_ROOT/output/step4_raw/homologous_frequency"
QMA="$PROJECT_ROOT/output/step4_raw/cdhit_cluster/qma_faa_clustered"
DB_ROOT="$PROJECT_ROOT/data/gtdb/diamond_db"
THREADS=12
WORKERS=4
SHARD_SIZE=128
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
  --step4_output) STEP4_OUTPUT="$2"; shift 2 ;;
  --diamond_db_dir) DB_ROOT="$2"; shift 2 ;;
  --threads) THREADS="$2"; shift 2 ;;
  --workers) WORKERS="$2"; shift 2 ;;
  --shard_size) SHARD_SIZE="$2"; shift 2 ;;
  --work_dir) WORK_DIR="$2"; shift 2 ;;
  -h|--help)
   echo "Usage: bash Step4_Quality_Control/run_replication_homologous_frequency.sh [options]"
   echo "  --step4_output PATH    Step4 formal output directory"
   echo "  --diamond_db_dir PATH  GTDB DIAMOND database directory"
   echo "  --threads N            Threads per DIAMOND worker (default: 12)"
   echo "  --workers N            Concurrent DIAMOND workers (default: 4)"
   echo "  --shard_size N         Unique queries per DIAMOND shard (default: 128)"
   echo "  --work_dir PATH        Internal shard/cache directory"
   exit 0 ;;
  *) echo "Unknown argument: $1" >&2; exit 2 ;;
 esac
done

if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"


run_stage "Step4/S09" --device CPU --threads "$THREADS" --counted --info "sensitive identity>=50% qcov>=70% tcov>=70% min_count=5" -- python -u "$SCRIPT_DIR/S09_replication_frequency.py" \
 --step4_output "$STEP4_OUTPUT" --qma "$QMA" \
 --bacteria_db "$DB_ROOT/gtdb_r220_bacteria" \
 --archaea_db "$DB_ROOT/gtdb_r220_archaea" \
 --output "$STEP4_OUTPUT/homologous_frequency" \
 --work_dir "$WORK_DIR" --threads "$THREADS" \
 --workers "$WORKERS" --shard_size "$SHARD_SIZE"
