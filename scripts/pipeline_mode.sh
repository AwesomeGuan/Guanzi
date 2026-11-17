#!/usr/bin/env bash
# Guard the shared Step 2-5 output tree against cross-mode reuse.

select_pipeline_mode() {
    local requested="$1"
    local output_root="$2"
    local state_file="$output_root/.active_mode"
    mkdir -p "$output_root"
    if [[ -f "$state_file" ]]; then
        local active
        active="$(tr -d '[:space:]' < "$state_file")"
        if [[ "$active" != "$requested" ]]; then
            echo "ERROR: shared Step 2-5 outputs belong to mode '$active', not '$requested'." >&2
            echo "Archive or remove the shared Step 2-5 outputs and output/.done, then remove output/.active_mode." >&2
            return 1
        fi
    else
        printf '%s\n' "$requested" > "$state_file"
    fi
}
