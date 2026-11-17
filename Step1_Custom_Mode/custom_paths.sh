#!/usr/bin/env bash
# Shared paths and completion markers for Custom Mode.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

QUERY_PDB="$PROJECT_ROOT/data/query_data/query_pdb"
OUTPUT_BASE="$PROJECT_ROOT/output"
PRODIGAL_OUT="$OUTPUT_BASE/prodigal_raw"
DB_PDB_OUT="$OUTPUT_BASE/custom_pdb"
CUSTOM_RAW="$OUTPUT_BASE/step1_raw"
SEARCH_OUT="$CUSTOM_RAW/foldseek_search_results"
FIDENT_OUT="$CUSTOM_RAW/fident_filter"
COVERAGE_OUT="$CUSTOM_RAW/coverage_filter"
MAPPING_OUT="$OUTPUT_BASE/step1_mapping"
MARKER_DIR="$OUTPUT_BASE/.done/custom/step1"
FOLDSEEK_DIR="$PROJECT_ROOT/foldseek"
DB_PREFIX="custom"
source "$PROJECT_ROOT/scripts/stage_state.sh"
