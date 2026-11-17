#!/usr/bin/env bash
# ============================================================
# run_foldseek_db.sh — Build FoldSeek database from PDB files
#
# Runs FoldSeek commands in the dedicated foldseek conda environment.
# Must cd into foldseek/ directory before running foldseek commands.
#
# Usage:
# bash Step1_Custom_Mode/run_foldseek_db.sh \
# --db_pdb output/custom_pdb --db_prefix mydb [options]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---------- defaults ----------
DB_PDB=""
DB_PREFIX=""
FOLDSEEK_DIR="$PROJECT_ROOT/foldseek"
THREADS=16

usage() {
 cat <<EOF
Usage: bash Step1_Custom_Mode/run_foldseek_db.sh [OPTIONS]

Required:
 --db_pdb PATH Directory containing .pdb files (recursive scan)
 --db_prefix NAME Database name (output: foldseek/{NAME}, foldseek/{NAME}_pad)

Optional:
 --threads N CPU threads (default: 16)

Example:
 bash Step1_Custom_Mode/run_foldseek_db.sh \\
 --db_pdb output/custom_pdb --db_prefix mydb
EOF
 exit 1
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
 case "$1" in
 --db_pdb) DB_PDB="$2"; shift 2 ;;
 --db_prefix) DB_PREFIX="$2"; shift 2 ;;
 --threads) THREADS="$2"; shift 2 ;;
 -h|--help) usage ;;
 *) echo "Unknown: $1"; usage ;;
 esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "foldseek" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" foldseek bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=custom
source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step1-Custom/FoldSeek-DB" CPU


# ---------- validate ----------
if [[ -z "$DB_PDB" ]] || [[ -z "$DB_PREFIX" ]]; then
 echo "ERROR: --db_pdb and --db_prefix are required"
 usage
fi

if [[ ! -d "$DB_PDB" ]]; then
 echo "ERROR: --db_pdb directory not found: $DB_PDB"
 exit 1
fi

mkdir -p "$FOLDSEEK_DIR"

# ---------- activate conda ----------
log() { echo "[FoldSeek] $*"; }

log "PDB input: $DB_PDB"
log "Database: $DB_PREFIX"
log "Output dir: $FOLDSEEK_DIR"


# Convert DB_PDB to absolute path before cd
DB_PDB="$(cd "$(dirname "$DB_PDB")" && pwd)/$(basename "$DB_PDB")"

# ---------- run (must cd to foldseek/ directory) ----------
log "Step 1/3: createdb..."
cd "$FOLDSEEK_DIR"

# Calculate relative path from foldseek/ to DB_PDB
DB_PDB_REL=$(python3 -c "import os; print(os.path.relpath('$DB_PDB', '$FOLDSEEK_DIR'))")

foldseek createdb "$DB_PDB_REL" "$DB_PREFIX" --threads "$THREADS"

log "Step 2/3: createindex..."
foldseek createindex "$DB_PREFIX" tmp --threads "$THREADS"

log "Step 3/3: makepaddedseqdb..."
foldseek makepaddedseqdb "$DB_PREFIX" "${DB_PREFIX}_pad" --threads "$THREADS"

log "============================================================"
log "FoldSeek DB build complete!"
log "============================================================"
log "Main DB: $(pwd)/$DB_PREFIX"
log "Pad DB: $(pwd)/${DB_PREFIX}_pad"
log ""
log "============================================================"
