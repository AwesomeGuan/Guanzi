#!/usr/bin/env bash
# Logging adapters and explicit zero-candidate guards; no completion markers are written here.

run_stage() {
    local stage_label="$1"
    shift
    python -u "$PROJECT_ROOT/scripts/logged_stage.py" \
        --label "$stage_label" --mode "${MODE:-unspecified}" "$@"
}

log_skip() {
    local stage_label="$1"
    shift
    python -u "$PROJECT_ROOT/scripts/logged_stage.py" \
        --label "$stage_label" --mode "${MODE:-unspecified}" --event SKIP "$@"
}

log_script_start() {
    SCRIPT_LOG_LABEL="$1"
    SCRIPT_LOG_DEVICE="$2"
    SCRIPT_LOG_GPUS="${3:-}"
    SCRIPT_LOG_STARTED="$(python -c 'import time; print(time.monotonic())')"
    python -u "$PROJECT_ROOT/scripts/logged_stage.py" --event START \
        --label "$SCRIPT_LOG_LABEL" --mode "${MODE:-unspecified}" \
        --device "$SCRIPT_LOG_DEVICE" --gpus "$SCRIPT_LOG_GPUS"
    trap 'log_script_exit "$?"' EXIT
}

log_script_exit() {
    local script_code="$1"
    python -u "$PROJECT_ROOT/scripts/logged_stage.py" --event END \
        --label "$SCRIPT_LOG_LABEL" --mode "${MODE:-unspecified}" \
        --device "$SCRIPT_LOG_DEVICE" --gpus "$SCRIPT_LOG_GPUS" \
        --started-at "$SCRIPT_LOG_STARTED" --exit-code "$script_code" || :
}

# Explicit output counts used by chapter orchestration, never completion markers.
candidate_count() {
    python "$PROJECT_ROOT/scripts/count_candidates.py" "$@"
}

stop_if_no_candidates() {
    local label="$1" after="$2" quantity="$3" count
    shift 3
    count="$(candidate_count "$@")" || return "$?"
    echo "[$label] SUMMARY ${quantity}=$count"
    if [[ "$count" -eq 0 ]]; then
        echo "[$label] STOP mode=${MODE:-unspecified} reason=no_candidates after=$after"
        exit 0
    fi
}
