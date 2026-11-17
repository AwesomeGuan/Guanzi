#!/usr/bin/env bash
# Apply the AFDB50 coverage/fident filters and extract candidate IDs (S02-S04).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
MIN_QCOV=0.70
MIN_TCOV=0.70
FIDENT_BLACKLIST_THRESHOLD=0.30
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"
ORIGINAL_ARGS=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --min_qcov) MIN_QCOV="$2"; shift 2 ;;
        --min_tcov) MIN_TCOV="$2"; shift 2 ;;
        --fident_blacklist_threshold) FIDENT_BLACKLIST_THRESHOLD="$2"; shift 2 ;;
        -h|--help)
            cat <<EOF
Usage: bash Step1_Replication_Mode/run_afdb50_filter.sh [OPTIONS]
  --min_qcov FLOAT        Minimum FoldSeek query coverage (default: 0.70)
  --min_tcov FLOAT        Minimum FoldSeek target coverage (default: 0.70)
  --fident_blacklist_threshold FLOAT
                          Remove targets with any coverage-passing hit above
                          this identity (default: 0.30)
EOF
            exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 2 ;;
    esac
done
if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"

mkdir -p "$REP_RAW" "$MARKER_DIR"
clear_filter_dependents() {
    rm -rf -- "$AFDB50_IDS_OUT" "$AFDB50_PDB_OUT" "$AFDB50_FAA_OUT" \
        "$DIAMOND_OUT" "$GTDB_MAPPING_OUT" "$CANDIDATE_GENOMES_OUT" \
        "$PRODIGAL_OUT" "$MAPPING_OUT" "$REP_RAW/.task_markers/diamond"
}
if stage_done 02_coverage_filter && \
        python "$VALIDATOR" coverage --input "$SEARCH_OUT" --output "$COVERAGE_OUT" \
        --min-qcov "$MIN_QCOV" --min-tcov "$MIN_TCOV" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S02" --device CPU
    echo "coverage filter already complete"
else
    invalidate_stages 02_coverage_filter 03_fident_filter 04_extract_afdb50_ids \
        05_download_isozyme_pdb 06_extract_isozyme_faa 07_query_homology_exclusion \
        08_gtdb_mapping 09_extract_candidate_contigs 10_materialize_gtdb_prodigal \
        11_build_replication_mapping
    clear_filter_dependents
    run_stage "Step1-Replication/S02" --device CPU --info "qcov>=${MIN_QCOV:-0.70} tcov>=${MIN_TCOV:-0.70}" -- python -u "$SCRIPT_DIR/S02_coverage_filter.py" --input "$SEARCH_OUT" \
        --output "$COVERAGE_OUT" --min_qcov "$MIN_QCOV" --min_tcov "$MIN_TCOV"
    python "$VALIDATOR" coverage --input "$SEARCH_OUT" --output "$COVERAGE_OUT" \
        --min-qcov "$MIN_QCOV" --min-tcov "$MIN_TCOV"
    mark_stage_done 02_coverage_filter
fi
if stage_done 03_fident_filter && \
        python "$VALIDATOR" fident --input "$COVERAGE_OUT" --output "$FIDENT_OUT" \
        --threshold "$FIDENT_BLACKLIST_THRESHOLD" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S03" --device CPU
    echo "fident filter already complete"
else
    invalidate_stages 03_fident_filter 04_extract_afdb50_ids \
        05_download_isozyme_pdb 06_extract_isozyme_faa 07_query_homology_exclusion \
        08_gtdb_mapping 09_extract_candidate_contigs 10_materialize_gtdb_prodigal \
        11_build_replication_mapping
    clear_filter_dependents
    run_stage "Step1-Replication/S03" --device CPU --info "target_blacklist: any fident>${FIDENT_BLACKLIST_THRESHOLD:-0.30}; then max-TM-score deduplication" -- python -u "$SCRIPT_DIR/S03_fident_filter.py" --input "$COVERAGE_OUT" \
        --output "$FIDENT_OUT" \
        --fident_blacklist_threshold "$FIDENT_BLACKLIST_THRESHOLD"
    python "$VALIDATOR" fident --input "$COVERAGE_OUT" --output "$FIDENT_OUT" \
        --threshold "$FIDENT_BLACKLIST_THRESHOLD"
    mark_stage_done 03_fident_filter
fi
if stage_done 04_extract_afdb50_ids && python "$VALIDATOR" ids \
        --input "$FIDENT_OUT" --output "$AFDB50_IDS_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S04" --device CPU
    echo "AFDB50 ID extraction already complete"
else
    invalidate_stages 04_extract_afdb50_ids 05_download_isozyme_pdb \
        06_extract_isozyme_faa 07_query_homology_exclusion 08_gtdb_mapping \
        09_extract_candidate_contigs 10_materialize_gtdb_prodigal \
        11_build_replication_mapping
    clear_filter_dependents
    run_stage "Step1-Replication/S04" --device CPU -- python -u "$SCRIPT_DIR/S04_extract_afdb50_ids.py" \
        --input "$FIDENT_OUT" --output "$AFDB50_IDS_OUT"
    python "$VALIDATOR" ids --input "$FIDENT_OUT" --output "$AFDB50_IDS_OUT"
    mark_stage_done 04_extract_afdb50_ids
fi
echo "AFDB50 filtering and ID extraction complete"
