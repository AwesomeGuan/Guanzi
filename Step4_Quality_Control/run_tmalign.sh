#!/usr/bin/env bash
# Stage 3: TM-align structural QC (S05-S08)
set -eu
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RAW="$PROJECT_ROOT/output/step4_raw"
CDHIT="$RAW/cdhit_cluster"
OUTPUT="$PROJECT_ROOT/output/step4_output/tmalign"
FORMAL_OUTPUT="$PROJECT_ROOT/output/step4_output"
QUERY_PDB="$PROJECT_ROOT/data/query_data/query_pdb"
MINING_PDB=""   # derived from --mode unless --mining_pdb overrides it
THREADS=12

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --mining_pdb) MINING_PDB="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash Step4_Quality_Control/run_tmalign.sh --mode custom|replication [OPTIONS]"
            echo "  --mode MODE         Input mode; selects the default Mining PDB directory"
            echo "  --mining_pdb PATH   Mining PDB directory override; default is derived from --mode"
            echo "  --output PATH       TM-align working directory"
            echo "  --threads N         Concurrent TM-align workers (default: 12)"
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }

# --mining_pdb overrides; otherwise the standard Mining structure library is
# derived from --mode so the mode alone selects the input path.
if [[ -z "$MINING_PDB" ]]; then
    case "$MODE" in
        custom) MINING_PDB="$PROJECT_ROOT/output/custom_pdb" ;;
        replication) MINING_PDB="$PROJECT_ROOT/output/step4_raw/mining_pdb" ;;
        *) echo "ERROR: --mode custom|replication is required unless --mining_pdb is given" >&2; exit 2 ;;
    esac
fi
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"


log() { echo "[TMalign] $*"; }
log "Mining PDB: $MINING_PDB"
mkdir -p "$OUTPUT"


log "S05: Preparing paired component PDBs ..."
run_stage "Step4/S05" --device CPU -- python -u "$SCRIPT_DIR/S05_prepare_pdb_pairs.py" \
 --qma "$CDHIT/qma_faa_clustered" \
 --mining_pdb "$MINING_PDB" --query_pdb "$QUERY_PDB" \
 --mining_output "$OUTPUT/mining_assembly_pdb" \
 --query_output "$OUTPUT/query_assembly_pdb"

log "S06: Run TM-align ..."
run_stage "Step4/S06" --device CPU --counted --threads "$THREADS" -- python -u "$SCRIPT_DIR/S06_tmalign_run.py" \
 --mining_assembly_pdb "$OUTPUT/mining_assembly_pdb" \
 --query_assembly_pdb "$OUTPUT/query_assembly_pdb" \
 --tmalign_results "$OUTPUT/tmalign_raw" \
 --threads "$THREADS"

log "S07: Building canonical component Comparison table ..."
run_stage "Step4/S07" --device CPU -- python -u "$SCRIPT_DIR/S07_build_comparison.py" \
 --qma "$CDHIT/qma_faa_clustered" \
 --assembly "$CDHIT/clustered_assembly.csv" \
 --tmalign_results "$OUTPUT/tmalign_raw" \
 --output "$OUTPUT/comparison.csv"

log "S08: Writing qualified/disqualified Mapping, Assembly and Comparison ..."
run_stage "Step4/S08" --device CPU --info "candidate_TM-score>=0.7 sequence_identity<=0.3" -- python -u "$SCRIPT_DIR/S08_finalize_step4.py" \
 --mapping "$CDHIT/clustered_mapping.csv" \
 --assembly "$CDHIT/clustered_assembly.csv" \
 --comparison "$OUTPUT/comparison.csv" \
 --output "$FORMAL_OUTPUT"

log "Done: $OUTPUT/"
