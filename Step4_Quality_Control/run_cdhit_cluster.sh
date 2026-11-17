#!/usr/bin/env bash
# Stage 2: CD-HIT clustering (S02-S04)
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RAW="$PROJECT_ROOT/output/step4_raw"
INPUT="$RAW/distance_filter"
OUTPUT="$RAW/cdhit_cluster"
PRODIGAL=""   # derived from --mode unless --prodigal_dir overrides it
QUERY_LOC="$PROJECT_ROOT/data/query_data/query_location"
THREADS=12

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --prodigal_dir) PRODIGAL="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash Step4_Quality_Control/run_cdhit_cluster.sh --mode custom|replication [OPTIONS]"
            echo "  --mode MODE          Input mode; selects the default Prodigal tree"
            echo "  --prodigal_dir PATH  Prodigal protein tree override; default is derived from --mode"
            echo "  --threads N          CD-HIT threads per process (default: 12)"
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }

# --prodigal_dir overrides; otherwise the standard protein tree is derived from
# --mode so the mode alone selects the input path.
if [[ -z "$PRODIGAL" ]]; then
    case "$MODE" in
        custom) PRODIGAL="$PROJECT_ROOT/output/prodigal_raw" ;;
        replication) PRODIGAL="$PROJECT_ROOT/output/rep_prodigal_raw" ;;
        *) echo "ERROR: --mode custom|replication is required unless --prodigal_dir is given" >&2; exit 2 ;;
    esac
fi
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


log() { echo "[CD-HIT] $*"; }
log "Prodigal: $PRODIGAL"
mkdir -p "$OUTPUT"

# Stage 2 owns this working tree: clear the previous run so stale instances or
# components cannot re-enter the CD-HIT pool or the paired PDB trees.
for STALE in qma_faa qma_faa_clustered clustering_temp \
             clustered_mapping.csv clustered_assembly.csv; do
    rm -rf -- "$OUTPUT/$STALE"
done
log "Cleared previous Stage 2 working tree: $OUTPUT"


log "S02: Building canonical Query/Mining assembly FAA tree ..."
run_stage "Step4/S02" --device CPU -- python -u "$SCRIPT_DIR/S02_extract_mining_faa.py" \
 --prodigal_dir "$PRODIGAL" \
 --assembly_csv "$INPUT/step4_distance_filtered_assembly.csv" \
 --qma_faa "$OUTPUT/qma_faa" \
 --mapping_csv "$INPUT/step4_distance_filtered_mapping.csv" \
 --query_location "$QUERY_LOC"

log "S03: CD-HIT clustering ..."
run_stage "Step4/S03" --device CPU --counted --threads "$THREADS" --info "identity>=0.5 cd_hit_threads=$THREADS" -- python -u "$SCRIPT_DIR/S03_cd_hit_clustering.py" \
 --cdhit_output "$OUTPUT/clustering_temp" \
 --qma_faa "$OUTPUT/qma_faa" \
 --threads "$THREADS"

log "S04: Selecting representative instances and formal tables ..."
run_stage "Step4/S04" --device CPU -- python -u "$SCRIPT_DIR/S04_build_clustered_tables.py" \
 --cdhit_output "$OUTPUT/clustering_temp" \
 --qma_faa "$OUTPUT/qma_faa" \
 --mapping_csv "$INPUT/step4_distance_filtered_mapping.csv" \
 --assembly_csv "$INPUT/step4_distance_filtered_assembly.csv" \
 --output "$OUTPUT"

log "Done: $OUTPUT/"
