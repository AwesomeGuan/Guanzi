#!/usr/bin/env bash
# ============================================================
# run_step3_assembly.sh - System assembly
#
# Usage:
# bash Step3_Assembly_Analysis/run_step3_assembly.sh \
# --mode replication \
# --input output/step2_output/step2_mapping_multi \
# --dict_dir output/step3_raw \
# --rules data/DefenseFinder_rules \
# --query_location data/query_data/query_location \
# --threads 12 \
# --output output/step3_output
#
# --mode selects the default protein tree (custom -> output/prodigal_raw,
# replication -> output/rep_prodigal_raw); --prodigal overrides it and is only
# needed when the tree lives outside the standard location.
# --diamond_bin is optional and resolved from defense_system when omitted.
# ============================================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT="$PROJECT_ROOT/output/step2_output/step2_mapping_multi"
DICT_DIR="$PROJECT_ROOT/output/step3_raw"
RULES="$PROJECT_ROOT/data/DefenseFinder_rules"
OUTPUT="$PROJECT_ROOT/output/step3_output"
QUERY_LOCATION="$PROJECT_ROOT/data/query_data/query_location"
PRODIGAL=""   # derived from --mode unless --prodigal overrides it
THREADS=12
DIAMOND_BIN=""

MODE=unspecified
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
        --mode) MODE="$2"; shift 2 ;;
 --input) INPUT="$2"; shift 2 ;;
 --dict_dir) DICT_DIR="$2"; shift 2 ;;
 --rules) RULES="$2"; shift 2 ;;
 --query_location) QUERY_LOCATION="$2"; shift 2 ;;
 --prodigal) PRODIGAL="$2"; shift 2 ;;
 --diamond_bin) DIAMOND_BIN="$2"; shift 2 ;;
 --threads) THREADS="$2"; shift 2 ;;
 --output) OUTPUT="$2"; shift 2 ;;
 -h|--help)
   echo "Usage: bash Step3_Assembly_Analysis/run_step3_assembly.sh [OPTIONS]"
   echo "  --mode MODE       Input mode: custom or replication (selects the default protein tree)"
   echo "  --threads N       DIAMOND threads (default: 12)"
   echo "  --prodigal PATH   Protein tree override; default is derived from --mode"
   echo "  --input PATH --dict_dir PATH --rules PATH --query_location PATH --output PATH"
   echo "  --diamond_bin PATH  Optional DIAMOND executable override"
   exit 0 ;;

 *) echo "Unknown option: $1" >&2; exit 2 ;;
 esac
done


[[ "$MODE" == "custom" || "$MODE" == "replication" || "$MODE" == "unspecified" ]] || { echo "ERROR: invalid --mode" >&2; exit 2; }

# --prodigal overrides; otherwise the standard protein tree is derived from
# --mode so the mode alone selects the input path.
if [[ -z "$PRODIGAL" ]]; then
    case "$MODE" in
        custom) PRODIGAL="$PROJECT_ROOT/output/prodigal_raw" ;;
        replication) PRODIGAL="$PROJECT_ROOT/output/rep_prodigal_raw" ;;
        *) echo "ERROR: --mode custom|replication is required unless --prodigal is given" >&2; exit 2 ;;
    esac
fi
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
source "$PROJECT_ROOT/scripts/stage_logging.sh"

stop_if_no_candidates "Step3-Assembly" "input" "multi_mapping_rows" csv --path "$INPUT" --required-column mining_protein_id


# S07 requires DIAMOND. Resolve it after entering defense_system so the exact
# executable used by the assembly stage is recorded and passed explicitly.

if [ -n "$DIAMOND_BIN" ]; then
 DIAMOND_BIN="$(command -v "$DIAMOND_BIN" 2>/dev/null || true)"
else
 DIAMOND_BIN="$(command -v diamond 2>/dev/null || true)"
fi
if [ -z "$DIAMOND_BIN" ] || [ ! -x "$DIAMOND_BIN" ]; then
 echo "ERROR: DIAMOND executable not found after activating defense_system; use --diamond_bin /absolute/path/to/diamond" >&2
 exit 1
fi

log() { echo "[Step3-Assembly] $*"; }
log "=== Step3: System Assembly ==="
log "Input: $INPUT"
log "Dict dir: $DICT_DIR"
log "Rules: $RULES"
log "Prodigal: $PRODIGAL"
log "DIAMOND: $DIAMOND_BIN"
log "Output: $OUTPUT"

MERGED_MAPPING="$DICT_DIR/step2_merged_mapping.csv"
mkdir -p "$OUTPUT"
mkdir -p "$DICT_DIR"

log "S00: Merging $INPUT -> mapping ..."
if ls "$INPUT"/*.csv &>/dev/null; then
 run_stage "Step3/S00" --device CPU -- python -u "$SCRIPT_DIR/S00_merge_step2_mappings.py" \
 --input "$INPUT" \
 --output_mapping "$MERGED_MAPPING"
else
 log "ERROR: no CSVs in $INPUT"
 exit 1
fi

stop_if_no_candidates "Step3-Assembly" "S00" "merged_mapping_rows" csv --path "$MERGED_MAPPING" --required-column mining_protein_id

log "S04: Independent novel-hit assembly ..."
run_stage "Step3/S04" --device CPU -- python -u "$SCRIPT_DIR/S04_independent_novel_hit_assembly.py" \
 --input "$MERGED_MAPPING" \
 --rules_dir "$RULES" \
 --output "$OUTPUT/independent_novel_hit_assembly"

log "S05: Partial complementation ..."
REMAIN="$OUTPUT/independent_novel_hit_assembly/remain_for_partial_complementation.csv"
REMAIN_ROWS="$(candidate_count csv --path "$REMAIN" --required-column mining_protein_id)"
echo "[Step3-Assembly] SUMMARY S05_remaining_mapping_rows=$REMAIN_ROWS"
if [ "$REMAIN_ROWS" -gt 0 ]; then
 run_stage "Step3/S05" --device CPU -- python -u "$SCRIPT_DIR/S05_partial_system_complementation.py" \
 --input "$REMAIN" \
 --rules_dir "$RULES" \
 --incomplete_dict "$DICT_DIR/defensefinder_incomplete_dict.pkl" \
 --output "$OUTPUT/partial_complementation"
else
 log "S05 skipped: no remaining candidates; existing assemblies are retained"
 mkdir -p "$OUTPUT/partial_complementation"
fi

log "S06: Complete expansion ..."
REMAIN="$OUTPUT/partial_complementation/remain_for_complete_expansion.csv"
if [ "$REMAIN_ROWS" -gt 0 ]; then
 REMAIN_ROWS="$(candidate_count csv --path "$REMAIN" --required-column mining_protein_id)"
fi
echo "[Step3-Assembly] SUMMARY S06_remaining_mapping_rows=$REMAIN_ROWS"
if [ "$REMAIN_ROWS" -gt 0 ]; then
 run_stage "Step3/S06" --device CPU -- python -u "$SCRIPT_DIR/S06_complete_system_expansion.py" \
 --input "$REMAIN" \
 --rules_dir "$RULES" \
 --systems_dict "$DICT_DIR/defensefinder_systems_dict.pkl" \
 --output "$OUTPUT/complete_expansion"
else
 log "S06 skipped: no remaining candidates; existing assemblies are retained"
 mkdir -p "$OUTPUT/complete_expansion"
fi

log "S07: Query-guided homologous assembly ..."
REMAIN="$OUTPUT/complete_expansion/remain_for_query_guided_assembly.csv"
if [ "$REMAIN_ROWS" -gt 0 ]; then
 REMAIN_ROWS="$(candidate_count csv --path "$REMAIN" --required-column mining_protein_id)"
fi
echo "[Step3-Assembly] SUMMARY S07_remaining_mapping_rows=$REMAIN_ROWS"
if [ "$REMAIN_ROWS" -gt 0 ]; then
 run_stage "Step3/S07" --device CPU --threads "$THREADS" --counted --info "sensitive identity>=50% qcov>=70% tcov>=70% evalue<=1e-5" -- python -u "$SCRIPT_DIR/S07_query_guided_assembly.py" \
 --input "$REMAIN" \
 --query_location "$QUERY_LOCATION" \
 --prodigal "$PRODIGAL" \
 --rules_dir "$RULES" \
 --output "$OUTPUT/query_guided_assembly" \
 --diamond_bin "$DIAMOND_BIN" \
 --threads "$THREADS"
else
 log "S07 skipped: no remaining candidates; existing assemblies are retained"
 mkdir -p "$OUTPUT/query_guided_assembly"
fi

log "=== Assembly complete ==="
log "Output directories: $OUTPUT/independent_novel_hit_assembly/ $OUTPUT/partial_complementation/ $OUTPUT/complete_expansion/ $OUTPUT/query_guided_assembly/"
