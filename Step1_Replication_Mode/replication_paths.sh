#!/usr/bin/env bash
# Shared path defaults for Replication Mode stage scripts.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

QUERY_PDB="$PROJECT_ROOT/data/query_data/query_pdb"
QUERY_LOCATION="$PROJECT_ROOT/data/query_data/query_location"
GTDB_GENOMES="$PROJECT_ROOT/data/gtdb/gtdb_genomes_reps_r220"
PROTEIN_FAA_REPS="$PROJECT_ROOT/data/gtdb/protein_faa_reps"
GTDB_DIAMOND_ROOT="$PROJECT_ROOT/data/gtdb/diamond_db"
OUTPUT_BASE="$PROJECT_ROOT/output"
AFDB50_PREFIX="afdb50"

REP_RAW="$OUTPUT_BASE/rep_step1_raw"
REP_MAP="$OUTPUT_BASE/rep_step1_mapping"
MARKER_DIR="$OUTPUT_BASE/.done/replication/step1"
DOWNLOAD_STATUS="$REP_RAW/afdb50_pdb_download_status.csv"

SEARCH_OUT="$REP_RAW/foldseek_search_results"
FIDENT_OUT="$REP_RAW/fident_filter"
COVERAGE_OUT="$REP_RAW/coverage_filter"
AFDB50_IDS_OUT="$REP_RAW/afdb50_isozyme_ids"
AFDB50_PDB_OUT="$REP_RAW/afdb50_isozyme_pdb"
AFDB50_FAA_OUT="$REP_RAW/afdb50_isozyme_faa"
QUERY_DIAMOND_DB="$PROJECT_ROOT/data/query_data/query_diamond_db"
DIAMOND_OUT="$REP_RAW/diamond_filtered_faa"
GTDB_MAPPING_OUT="$REP_RAW/gtdb_mapping"
MAPPING_OUT="$REP_MAP"
CANDIDATE_GENOMES_OUT="$REP_RAW/candidate_genomes"
PRODIGAL_OUT="$OUTPUT_BASE/rep_prodigal_raw"
source "$PROJECT_ROOT/scripts/stage_state.sh"

require_resource() {
    if [[ ! -d "$1" ]]; then
        echo "ERROR: resource directory not found: $1" >&2
        exit 1
    fi
}
