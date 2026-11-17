#!/usr/bin/env bash
# Shared empty completion-marker helpers. MARKER_DIR must be set by the caller.

stage_marker() {
    printf '%s/%s.done\n' "$MARKER_DIR" "$1"
}

stage_done() {
    [[ -f "$(stage_marker "$1")" ]]
}

mark_stage_done() {
    mkdir -p "$MARKER_DIR"
    : > "$(stage_marker "$1")"
}

invalidate_stages() {
    local stage
    for stage in "$@"; do
        rm -f -- "$(stage_marker "$stage")"
    done
}
