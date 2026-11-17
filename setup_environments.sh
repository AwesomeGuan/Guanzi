#!/usr/bin/env bash
# Create and verify the six Linux conda environments used by the pipeline.
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DONE_DIR="$PROJECT_ROOT/output/.done/setup_envs"
mkdir -p "$DONE_DIR"

has_env() {
    conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$1"
}

verify_env() {
    local name="$1"
    shift
    conda run --no-capture-output -n "$name" "$@"
}

verify_defense_system() {
    verify_env defense_system python -c 'import Bio, pandas, requests, tqdm' &&
        verify_env defense_system diamond help &&
        verify_env defense_system bash -c 'command -v cd-hit >/dev/null 2>&1' &&
        verify_env defense_system TMalign
}

run_step() {
    local number="$1"
    local name="$2"
    local verifier="$3"
    shift 3
    local marker="$DONE_DIR/${number}_${name}.done"
    local validation_output
    if has_env "$name"; then
        if validation_output="$(eval "$verifier" 2>&1)"; then
            : > "$marker"
            echo "[$number/6] $name: existing environment verified; marker written; skipped"
            return
        fi
        rm -f -- "$marker"
        [[ -z "$validation_output" ]] || printf '%s\n' "$validation_output" >&2
        echo "ERROR: existing conda environment '$name' failed validation." >&2
        echo "Remove or repair it manually, then run setup_environments.sh again:" >&2
        echo "  conda env remove -n $name -y" >&2
        echo "  bash setup_environments.sh" >&2
        exit 1
    fi
    rm -f -- "$marker"
    echo "[$number/6] $name: creating and configuring"
    "$@"
    if ! validation_output="$(eval "$verifier" 2>&1)"; then
        [[ -z "$validation_output" ]] || printf '%s\n' "$validation_output" >&2
        echo "ERROR: newly created conda environment '$name' failed validation." >&2
        echo "No completion marker was written." >&2
        exit 1
    fi
    : > "$marker"
    echo "[$number/6] $name: complete"
}

create_defense_system() {
    conda create -n defense_system -y -c conda-forge -c bioconda \
        python=3.10 pandas biopython requests tqdm diamond cd-hit tmalign
}

create_prodigal() {
    conda create -n prodigal -y -c conda-forge -c bioconda python=3.10 prodigal
}

create_esmfold() {
    conda create -n esmfold -y python=3.9
    conda install -n esmfold -y -c pytorch -c conda-forge \
        pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 \
        cudatoolkit=11.6 numpy=1.24.3 pandas scipy scikit-learn biopython
    conda install -n esmfold -y -c conda-forge -c huggingface \
        transformers==4.27.0 accelerate sentencepiece tqdm requests protobuf
    conda run -n esmfold python -m pip install --force-reinstall safetensors==0.3.3
    conda run -n esmfold python -m pip install accelerate==0.24.1 modelscope
}

create_foldseek() {
    conda create -n foldseek -y -c conda-forge -c bioconda \
        python=3.10 pandas requests biopython foldseek
}

create_defensefinder() {
    conda create -n defensefinder -y -c conda-forge -c bioconda \
        python=3.12 pip hmmer pandas numpy
    conda run -n defensefinder python -m pip install mdmparis-defense-finder
}

run_step 01 defense_system \
    "verify_defense_system" \
    create_defense_system
run_step 02 prodigal \
    "verify_env prodigal python -c 'import sys, time; assert sys.version_info[:2] == (3, 10), sys.version; assert hasattr(time, \"monotonic\")' && verify_env prodigal prodigal -v" \
    create_prodigal
run_step 03 esmfold \
    "verify_env esmfold python -c 'import modelscope; from transformers import EsmForProteinFolding'" \
    create_esmfold
run_step 04 foldseek "verify_env foldseek foldseek -h" create_foldseek
run_step 05 defensefinder \
    "verify_env defensefinder defense-finder run --help" create_defensefinder

create_eggnog() {
    conda create -n eggnog -y -c conda-forge -c bioconda \
        python=3.10 eggnog-mapper=2.1.13 diamond=2.2.8
}

verify_eggnog() {
    verify_env eggnog python -c 'import sys, sqlite3, Bio, psutil, xlsxwriter; from eggnogmapper.version import __VERSION__; assert sys.version_info[:2] == (3, 10), sys.version; assert __VERSION__ == "2.1.13", __VERSION__' &&
        verify_env eggnog emapper.py --help &&
        verify_env eggnog bash -c 'test "$(diamond version)" = "diamond version 2.2.8"'
}

run_step 06 eggnog "verify_eggnog" create_eggnog

echo "All pipeline environments are ready."
