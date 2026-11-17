#!/usr/bin/env bash
# Unpack the per-genome GTDB representative protein files in place.
set -euo pipefail

if [[ $# -gt 0 ]]; then
    if [[ $# -eq 1 && ( "$1" == -h || "$1" == --help ) ]]; then
        echo "Usage: bash scripts/unpack_gtdb_proteins.sh"
        echo "Unpack data/gtdb/protein_faa_reps/{bacteria,archaea}/*_protein.faa.gz."
        exit 0
    fi
    echo "ERROR: this resource operation takes no arguments" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROTEIN_ROOT="$PROJECT_ROOT/data/gtdb/protein_faa_reps"
START_SECONDS=$SECONDS
LAST_REPORT=$SECONDS
TEMPORARY=""
TOTAL=0
PROCESSED=0
UNPACKED=0
REUSED=0
FAILURES=0

cleanup() {
    if [[ -n "$TEMPORARY" && -f "$TEMPORARY" ]]; then
        rm -f -- "$TEMPORARY"
    fi
}

finish() {
    local result=$?
    trap - EXIT
    cleanup
    echo "[GTDB-Proteins] END total_archives=$TOTAL processed=$PROCESSED unpacked=$UNPACKED reused_faa=$REUSED failures=$FAILURES remaining=$((TOTAL - UNPACKED - REUSED)) elapsed=$((SECONDS - START_SECONDS))s exit_code=$result"
    exit "$result"
}

trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

echo "[GTDB-Proteins] START device=CPU operation=unpack root=$PROTEIN_ROOT"
for binary in gzip cmp mktemp mv rm; do
    command -v "$binary" >/dev/null 2>&1 || {
        echo "[GTDB-Proteins] ERROR required executable not found: $binary" >&2
        exit 1
    }
done
for domain in bacteria archaea; do
    [[ -d "$PROTEIN_ROOT/$domain" ]] || {
        echo "[GTDB-Proteins] ERROR protein directory not found: $PROTEIN_ROOT/$domain" >&2
        exit 1
    }
done

shopt -s nullglob
ARCHIVES=("$PROTEIN_ROOT/bacteria/"*_protein.faa.gz "$PROTEIN_ROOT/archaea/"*_protein.faa.gz)
PLAIN_FILES=("$PROTEIN_ROOT/bacteria/"*_protein.faa "$PROTEIN_ROOT/archaea/"*_protein.faa)
TOTAL=${#ARCHIVES[@]}
echo "[GTDB-Proteins] INFO total_archives=$TOTAL existing_faa=${#PLAIN_FILES[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
    if [[ ${#PLAIN_FILES[@]} -eq 0 ]]; then
        echo "[GTDB-Proteins] ERROR no representative protein FAA or FAA.gz files found" >&2
        exit 1
    fi
    echo "[GTDB-Proteins] INFO no compressed protein files remain; existing FAA files retained"
    exit 0
fi

for archive in "${ARCHIVES[@]}"; do
    target="${archive%.gz}"
    PROCESSED=$((PROCESSED + 1))
    TEMPORARY="$(mktemp "${target%/*}/.unpack_gtdb_proteins.XXXXXX")"
    if ! gzip -dc -- "$archive" > "$TEMPORARY"; then
        echo "[GTDB-Proteins] ERROR decompression failed; source retained: $archive" >&2
        FAILURES=$((FAILURES + 1))
        cleanup
        TEMPORARY=""
        continue
    fi
    if [[ ! -s "$TEMPORARY" ]]; then
        echo "[GTDB-Proteins] ERROR empty decoded FAA; source retained: $archive" >&2
        FAILURES=$((FAILURES + 1))
        cleanup
        TEMPORARY=""
        continue
    fi
    if [[ -e "$target" || -L "$target" ]]; then
        if [[ ! -f "$target" ]] || ! cmp -s -- "$TEMPORARY" "$target"; then
            echo "[GTDB-Proteins] ERROR existing FAA differs; both original files retained: $target" >&2
            FAILURES=$((FAILURES + 1))
            cleanup
            TEMPORARY=""
            continue
        fi
        cleanup
        TEMPORARY=""
        if rm -- "$archive"; then
            REUSED=$((REUSED + 1))
        else
            FAILURES=$((FAILURES + 1))
        fi
    else
        # Publish without replacing a file created by another process.
        if ! mv -n -- "$TEMPORARY" "$target" || [[ -e "$TEMPORARY" ]]; then
            echo "[GTDB-Proteins] ERROR FAA publication failed; source retained: $archive" >&2
            FAILURES=$((FAILURES + 1))
            cleanup
            TEMPORARY=""
            continue
        fi
        TEMPORARY=""
        if rm -- "$archive"; then
            UNPACKED=$((UNPACKED + 1))
        else
            FAILURES=$((FAILURES + 1))
        fi
    fi
    if (( SECONDS - LAST_REPORT >= 300 )); then
        echo "[GTDB-Proteins] PROGRESS total_archives=$TOTAL processed=$PROCESSED unpacked=$UNPACKED reused_faa=$REUSED failures=$FAILURES elapsed=$((SECONDS - START_SECONDS))s"
        LAST_REPORT=$SECONDS
    fi
done

if [[ "$FAILURES" -gt 0 ]]; then
    exit 1
fi
