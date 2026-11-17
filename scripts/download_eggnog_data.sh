#!/usr/bin/env bash
# Deploy only the eggNOG 5.0.2 resources required by the DIAMOND route.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
[[ $# -eq 0 ]] || { echo "Usage: bash scripts/download_eggnog_data.sh" >&2; exit 2; }
if [[ "${CONDA_DEFAULT_ENV:-}" != eggnog ]]; then
    exec bash "$SCRIPT_DIR/run_in_env.sh" eggnog bash "$0" "$@"
fi
MODE=setup
source "$SCRIPT_DIR/stage_logging.sh"
log_script_start eggNOG-Resources CPU
DATA_DIR="$PROJECT_ROOT/data/eggnog_5.0.2"
BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
mkdir -p "$DATA_DIR"
sqlite_valid() {
    python -c 'import sqlite3,sys; from pathlib import Path; p=Path(sys.argv[1]); db=sqlite3.connect(p.resolve().as_uri()+"?mode=ro",uri=True); assert db.execute("SELECT name FROM sqlite_master WHERE type=\"table\" LIMIT 1").fetchone(), "No SQLite tables"; db.close()' "$1"
}
fetch() {
    local name="$1"
    wget -c --tries=3 --timeout=60 "$BASE_URL/$name" -O "$DATA_DIR/$name"
}
for NAME in eggnog.db eggnog_proteins.dmnd; do
    VALID=false
    if [[ -s "$DATA_DIR/$NAME" ]]; then
        if [[ "$NAME" == eggnog.db ]]; then
            sqlite_valid "$DATA_DIR/$NAME" && VALID=true
        else
            diamond dbinfo --db "$DATA_DIR/$NAME" && VALID=true
        fi
    fi
    if [[ "$VALID" == true ]]; then
        echo "[eggNOG-Resources] SKIP verified=$NAME"
        continue
    fi
    [[ ! -e "$DATA_DIR/$NAME" ]] || { echo "ERROR: invalid existing $NAME; remove or repair that file manually" >&2; exit 1; }
    [[ -s "$DATA_DIR/$NAME.gz" ]] && gzip -t "$DATA_DIR/$NAME.gz" || fetch "$NAME.gz"
    gzip -t "$DATA_DIR/$NAME.gz"
    TEMPORARY="$DATA_DIR/.unpacking.$NAME"
    gzip -dc "$DATA_DIR/$NAME.gz" > "$TEMPORARY"
    if [[ "$NAME" == eggnog.db ]]; then
        sqlite_valid "$TEMPORARY"
    else
        diamond dbinfo --db "$TEMPORARY"
    fi
    mv -- "$TEMPORARY" "$DATA_DIR/$NAME"
done
if [[ -s "$DATA_DIR/eggnog.taxa.db" ]] && sqlite_valid "$DATA_DIR/eggnog.taxa.db"; then
    echo "[eggNOG-Resources] SKIP verified=eggnog.taxa.db"
else
    [[ ! -e "$DATA_DIR/eggnog.taxa.db" ]] || { echo "ERROR: invalid existing eggnog.taxa.db; repair it manually" >&2; exit 1; }
    [[ -s "$DATA_DIR/eggnog.taxa.tar.gz" ]] && tar -tzf "$DATA_DIR/eggnog.taxa.tar.gz" >/dev/null || fetch eggnog.taxa.tar.gz
    tar -tzf "$DATA_DIR/eggnog.taxa.tar.gz" >/dev/null
    python - "$DATA_DIR" <<'PY'
import sys, tarfile
from pathlib import Path
root = Path(sys.argv[1]).resolve()
with tarfile.open(root / "eggnog.taxa.tar.gz") as archive:
    for member in archive.getmembers():
        target = (root / member.name).resolve()
        if not target.is_relative_to(root) or member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
            raise ValueError(f"Unsafe archive member: {member.name}")
        if member.isfile() and target.exists():
            raise ValueError(f"Refusing to overwrite an existing resource: {target}")
    archive.extractall(root, filter="data")
PY
    sqlite_valid "$DATA_DIR/eggnog.taxa.db"
fi
VERSION="$(emapper.py --version --data_dir "$DATA_DIR")"
printf '%s\n' "$VERSION"
[[ "$VERSION" == *"Installed eggNOG DB version: 5.0.2 /"* ]] || { echo "ERROR: database version is not 5.0.2" >&2; exit 1; }
echo "[eggNOG-Resources] SUMMARY ready=$DATA_DIR"
