#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/replication_paths.sh"
VALIDATOR="$PROJECT_ROOT/scripts/validate_step1_stage.py"

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            echo "Usage: bash Step1_Replication_Mode/run_gtdb_prodigal_materialization.sh"
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "defense_system" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" defense_system bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=replication
source "$PROJECT_ROOT/scripts/stage_logging.sh"


require_resource "$PROTEIN_FAA_REPS"
require_resource "$GTDB_GENOMES"


if stage_done 09_extract_candidate_contigs && python "$VALIDATOR" contigs \
        --mapping "$GTDB_MAPPING_OUT" --output "$CANDIDATE_GENOMES_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S09" --device CPU
    echo "Candidate contig extraction already complete"
else
    invalidate_stages 09_extract_candidate_contigs
    run_stage "Step1-Replication/S09" --device CPU --counted -- python -u "$SCRIPT_DIR/S09_extract_candidate_contigs.py" \
        --mapping "$GTDB_MAPPING_OUT" \
        --gtdb_genomes "$GTDB_GENOMES" \
        --output "$CANDIDATE_GENOMES_OUT"
    python "$VALIDATOR" contigs --mapping "$GTDB_MAPPING_OUT" \
        --output "$CANDIDATE_GENOMES_OUT"
    mark_stage_done 09_extract_candidate_contigs
fi

if stage_done 10_materialize_gtdb_prodigal && python "$VALIDATOR" materialized \
        --mapping "$GTDB_MAPPING_OUT" --contigs "$CANDIDATE_GENOMES_OUT" \
        --output "$PRODIGAL_OUT" >/dev/null 2>&1; then
    log_skip "Step1-Replication/S10" --device CPU
    echo "GTDB Prodigal materialization already complete"
    exit 0
fi

invalidate_stages 10_materialize_gtdb_prodigal
run_stage "Step1-Replication/S10" --device CPU --counted -- python -u "$SCRIPT_DIR/S10_materialize_gtdb_prodigal.py" \
    --gtdb_mapping "$GTDB_MAPPING_OUT" \
    --protein_faa_reps "$PROTEIN_FAA_REPS" \
    --candidate_genomes "$CANDIDATE_GENOMES_OUT" \
    --output "$PRODIGAL_OUT"

python "$VALIDATOR" materialized --mapping "$GTDB_MAPPING_OUT" \
    --contigs "$CANDIDATE_GENOMES_OUT" --output "$PRODIGAL_OUT"
mark_stage_done 10_materialize_gtdb_prodigal
