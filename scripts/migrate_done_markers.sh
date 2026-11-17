#!/usr/bin/env bash
# Copy legacy Step 1 markers into the centralized empty-marker layout.
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

migrate_directory() {
    local source="$1"
    local destination="$2"
    [[ -d "$source" ]] || return 0
    mkdir -p "$destination"
    local marker
    for marker in "$source"/*.done; do
        [[ -f "$marker" ]] || continue
        : > "$destination/$(basename "$marker")"
        echo "Migrated: $(basename "$marker")"
    done
}

migrate_directory "$PROJECT_ROOT/output/step1_raw/.done" \
    "$PROJECT_ROOT/output/.done/custom/step1"
migrate_directory "$PROJECT_ROOT/output/rep_step1_raw/.done" \
    "$PROJECT_ROOT/output/.done/replication/step1"

echo "Legacy Step 1 markers were copied. Stage wrappers will validate outputs before reuse."
