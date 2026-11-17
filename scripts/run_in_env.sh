#!/usr/bin/env bash
# Run one pipeline command inside a named conda environment.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ $# -lt 2 ]]; then
    echo "Usage: bash scripts/run_in_env.sh ENV_NAME COMMAND [ARG ...]" >&2
    exit 2
fi

ENV_NAME="$1"
shift

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available in PATH" >&2
    exit 1
fi

CONDA_BASE="$(conda info --base 2>/dev/null)" || {
    echo "ERROR: unable to locate the conda installation" >&2
    exit 1
}
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
if [[ ! -f "$CONDA_SH" ]]; then
    echo "ERROR: conda activation script not found: $CONDA_SH" >&2
    exit 1
fi

set +u
source "$CONDA_SH"
conda activate "$ENV_NAME" || {
    echo "ERROR: failed to activate conda environment: $ENV_NAME" >&2
    exit 1
}
set -u

cd "$PROJECT_ROOT"
exec "$@"
