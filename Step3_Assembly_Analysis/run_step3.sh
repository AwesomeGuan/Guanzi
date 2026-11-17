#!/usr/bin/env bash
# Run Step 3 dictionaries and assembly with mode-specific Prodigal inputs.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODE=""
THREADS=12

usage() {
    cat <<EOF
Usage: bash Step3_Assembly_Analysis/run_step3.sh --mode custom|replication [OPTIONS]
  --threads N      DIAMOND threads for query-guided assembly (default: 12)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done
[[ "$MODE" == "custom" || "$MODE" == "replication" ]] || { usage; exit 2; }


source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step3" "CPU" "${GPUS:-}"
source "$PROJECT_ROOT/scripts/pipeline_mode.sh"
select_pipeline_mode "$MODE" "$PROJECT_ROOT/output"
MARKER_DIR="$PROJECT_ROOT/output/.done/$MODE/step3"
source "$PROJECT_ROOT/scripts/stage_state.sh"
mkdir -p "$MARKER_DIR"

if [[ "$MODE" == "custom" ]]; then
    PRODIGAL="$PROJECT_ROOT/output/prodigal_raw"
else
    PRODIGAL="$PROJECT_ROOT/output/rep_prodigal_raw"
fi
DICT_DIR="$PROJECT_ROOT/output/step3_raw"
ASSEMBLY="$PROJECT_ROOT/output/step3_output"

stop_if_no_candidates "Step3" "input" "multi_mapping_rows" csv --path "$PROJECT_ROOT/output/step2_output/step2_mapping_multi" --required-column mining_protein_id


if stage_done 01_dictionaries \
        && [[ -s "$DICT_DIR/defensefinder_genes_dict.pkl" ]] \
        && [[ -s "$DICT_DIR/defensefinder_systems_dict.pkl" ]] \
        && [[ -s "$DICT_DIR/defensefinder_incomplete_dict.pkl" ]]; then
    log_skip "Step3/S01-S03" --device CPU
    echo "Step 3 dictionaries: already complete"
else
    invalidate_stages 01_dictionaries
    bash "$SCRIPT_DIR/run_step3_dict.sh"
    [[ -s "$DICT_DIR/defensefinder_genes_dict.pkl" ]] || exit 1
    [[ -s "$DICT_DIR/defensefinder_systems_dict.pkl" ]] || exit 1
    [[ -s "$DICT_DIR/defensefinder_incomplete_dict.pkl" ]] || exit 1
    mark_stage_done 01_dictionaries
fi

if stage_done 02_assembly \
        && [[ -d "$ASSEMBLY/independent_novel_hit_assembly" ]] \
        && [[ -d "$ASSEMBLY/partial_complementation" ]] \
        && [[ -d "$ASSEMBLY/complete_expansion" ]] \
        && [[ -d "$ASSEMBLY/query_guided_assembly" ]]; then
    log_skip "Step3/S04-S07" --device CPU
    echo "Step 3 assembly: already complete"
else
    invalidate_stages 02_assembly
    bash "$SCRIPT_DIR/run_step3_assembly.sh" --mode "$MODE" \
        --prodigal "$PRODIGAL" --threads "$THREADS"
    [[ -d "$ASSEMBLY/query_guided_assembly" ]] || exit 1
    mark_stage_done 02_assembly
fi

ASSEMBLY_ROWS="$(candidate_count csv --path "$ASSEMBLY" --pattern "*/assembly.csv" --required-column mining_contig_system_index)"
echo "[Step3] SUMMARY assembled_instance_rows=$ASSEMBLY_ROWS"


echo "Step 3 complete for $MODE mode"
