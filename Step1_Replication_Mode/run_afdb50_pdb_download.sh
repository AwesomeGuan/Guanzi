#!/usr/bin/env bash
# Download AFDB50 PDBs and derive FAA files with count-based checks (S05-S06).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
THREADS=8
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"
ORIGINAL_ARGS=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash Step1_Replication_Mode/run_afdb50_pdb_download.sh [--threads N]"
            echo "  --threads N  Download workers (default: 8)"
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
if stage_done 05_download_isozyme_pdb && python "$VALIDATOR" pdb \
        --input "$AFDB50_IDS_OUT" --output "$AFDB50_PDB_OUT" \
        --status "$DOWNLOAD_STATUS" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S05" --device CPU --threads "$THREADS"
    echo "AFDB50 PDB download accounting already complete"
else
    python "$PROJECT_ROOT/scripts/preflight.py" network --url "https://alphafold.ebi.ac.uk/"
    invalidate_stages 05_download_isozyme_pdb 06_extract_isozyme_faa
    run_stage "Step1-Replication/S05-reconcile" --device CPU -- python -u "$SCRIPT_DIR/S05_reconcile_pdb_state.py" \
        --input "$AFDB50_IDS_OUT" --pdb "$AFDB50_PDB_OUT" \
        --status "$DOWNLOAD_STATUS"
    run_stage "Step1-Replication/S05" --device CPU --threads "$THREADS" --counted -- python -u "$SCRIPT_DIR/S05_download_isozyme_pdb.py" \
        --input "$AFDB50_IDS_OUT" --output "$AFDB50_PDB_OUT" \
        --status "$DOWNLOAD_STATUS" --threads "$THREADS"
    python "$VALIDATOR" pdb --input "$AFDB50_IDS_OUT" \
        --output "$AFDB50_PDB_OUT" --status "$DOWNLOAD_STATUS"
    mark_stage_done 05_download_isozyme_pdb
fi

stop_if_no_candidates "Step1-Replication" "S05" "successful_pdbs" pdb-status --path "$DOWNLOAD_STATUS" --pdb-root "$AFDB50_PDB_OUT"

if stage_done 06_extract_isozyme_faa && python "$VALIDATOR" faa \
        --input "$AFDB50_PDB_OUT" --output "$AFDB50_FAA_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S06" --device CPU
    echo "AFDB50 FAA extraction already complete"
else
    invalidate_stages 06_extract_isozyme_faa
    run_stage "Step1-Replication/S06-reconcile" --device CPU -- python -u "$SCRIPT_DIR/S06_reconcile_faa_state.py" \
        --pdb "$AFDB50_PDB_OUT" --faa "$AFDB50_FAA_OUT"
    run_stage "Step1-Replication/S06" --device CPU --counted -- python -u "$SCRIPT_DIR/S06_extract_isozyme_faa.py" \
        --input "$AFDB50_PDB_OUT" --output "$AFDB50_FAA_OUT"
    python "$VALIDATOR" faa --input "$AFDB50_PDB_OUT" --output "$AFDB50_FAA_OUT"
    mark_stage_done 06_extract_isozyme_faa
fi
echo "AFDB50 PDB download and FAA extraction complete"
