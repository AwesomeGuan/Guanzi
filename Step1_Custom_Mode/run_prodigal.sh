#!/usr/bin/env bash
# ============================================================
# run_prodigal.sh  - Annotate user genome .fna files with Prodigal
#
# Multi-contig .fasta files are auto-split: each contig gets its
# own Prodigal run. Output matches Prodigal-style directory structure.
#
# Usage:
#   bash Step1_Custom_Mode/run_prodigal.sh --genome_dir input/genomes --outdir output/prodigal_raw
#
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

GENOME_DIR=""
OUTDIR=""

usage() {
    cat <<EOF
Usage: bash Step1_Custom_Mode/run_prodigal.sh [OPTIONS]

Required:
  --genome_dir PATH   Directory containing .fna/.fasta genome files
  --outdir PATH       Prodigal output directory

Each input file is treated as one genome (may contain multiple contigs).
Output structure: <outdir>/<genome>/<contig>/{.faa, .ffn, .gff}
EOF
    exit 1
}

ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome_dir) GENOME_DIR="$2"; shift 2 ;;
        --outdir)     OUTDIR="$2";     shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown: $1"; usage ;;
    esac
done


if [[ "${CONDA_DEFAULT_ENV:-}" != "prodigal" ]]; then
    exec bash "$PROJECT_ROOT/scripts/run_in_env.sh" prodigal bash "$0" ${ORIGINAL_ARGS[@]+"${ORIGINAL_ARGS[@]}"}
fi
MODE=custom
source "$PROJECT_ROOT/scripts/stage_logging.sh"
log_script_start "Step1-Custom/Prodigal" CPU


# ---------- validate ----------
if [[ -z "$GENOME_DIR" ]] || [[ -z "$OUTDIR" ]]; then
    echo "ERROR: --genome_dir and --outdir are required"
    usage
fi

mkdir -p "$OUTDIR"

log() { echo "[Prodigal] $*"; }

log "Genome dir: $GENOME_DIR"
log "Output:     $OUTDIR"

# ---------- activate conda ----------

count=0
total_contigs=0

for fna in "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fa; do
    [[ -f "$fna" ]] || continue

    genome_name=$(basename "$fna")
    genome_name="${genome_name%.fna}"
    genome_name="${genome_name%.fasta}"
    genome_name="${genome_name%.fa}"
    n_contigs=$(grep -c "^>" "$fna" || echo 0)

    if [[ "$n_contigs" -eq 1 ]]; then
        contig_name=$(grep "^>" "$fna" | head -1 | sed 's/^>//; s/ .*//')
        out_dir="$OUTDIR/$genome_name/$contig_name"

        if [[ -f "$out_dir/$contig_name.faa" ]]; then
            log "[$contig_name] Already done  - skipping"
        else
            log "[$contig_name] single contig  - running Prodigal ..."
            mkdir -p "$out_dir"

            prodigal -i "$fna" -o "$out_dir/$contig_name.gff" \
                -a "$out_dir/$contig_name.faa" \
                -d "$out_dir/$contig_name.ffn" \
                -f gff -p meta -q

            # Post-process: fix GFF gene IDs  - Prodigal assigns global IDs (e.g. ID=1_1)
            # in a single-contig run; replace with contig-prefixed IDs for consistency.
            sed -i -E "s/ID=1_([0-9]+)/ID=${contig_name}_\1/g" "$out_dir/$contig_name.gff"
            sed -i -E "s/>(1_[0-9]+)/>${contig_name}_\1/" "$out_dir/$contig_name.faa"
            sed -i -E "s/>(1_[0-9]+)/>${contig_name}_\1/" "$out_dir/$contig_name.ffn"

            log "[$contig_name] Done  - $out_dir/"
        fi
        count=$((count + 1))
        total_contigs=$((total_contigs + 1))

    else
        log "Split: $(basename "$fna") ($n_contigs contigs)"
        TMP_DIR=$(mktemp -d)

        python -c "
import sys
with open(sys.argv[1]) as fh:
    content = fh.read()
for record in content.split('\n>'):
    if not record.startswith('>'):
        record = '>' + record
    header_line = record.split('\n', 1)[0]
    contig_name = header_line[1:].split()[0]
    out_text = record if record.endswith('\n') else record + '\n'
    with open(sys.argv[2] + '/' + contig_name + '.fna', 'w') as wf:
        wf.write(out_text)
" "$fna" "$TMP_DIR"

        for split_fna in "$TMP_DIR"/*.fna; do
            contig_name=$(basename "$split_fna" .fna)
            out_dir="$OUTDIR/$genome_name/$contig_name"

            if [[ -f "$out_dir/$contig_name.faa" ]]; then
                log "  [$contig_name] Already done  - skipping"
            else
                log "  [$contig_name] Prodigal ..."
                mkdir -p "$out_dir"

                prodigal -i "$split_fna" -o "$out_dir/$contig_name.gff" \
                    -a "$out_dir/$contig_name.faa" \
                    -d "$out_dir/$contig_name.ffn" \
                    -f gff -p meta -q

                # Fix gene IDs
                sed -i -E "s/ID=1_([0-9]+)/ID=${contig_name}_\1/g" "$out_dir/$contig_name.gff"
                sed -i -E "s/>(1_[0-9]+)/>${contig_name}_\1/" "$out_dir/$contig_name.faa"
                sed -i -E "s/>(1_[0-9]+)/>${contig_name}_\1/" "$out_dir/$contig_name.ffn"

                log "  [$contig_name] Done  - $out_dir/"
            fi
            rm -f "$split_fna"
            total_contigs=$((total_contigs + 1))
        done

        rmdir "$TMP_DIR" 2>/dev/null
        count=$((count + 1))
    fi
done

log "=== $count file(s)  - $total_contigs contig(s) annotated ==="

ABS_OUTPUT="$(cd "$OUTDIR" && pwd)"
echo ""
echo "============================================================"
echo "  Prodigal annotation complete!"
echo "============================================================"
echo ""
echo "  Output: $ABS_OUTPUT"
echo "  Structure: <output>/<genome>/<contig>/{.faa, .ffn, .gff}"
echo ""
echo "  To reuse existing gene predictions, pass their directory with"
echo "  --protein_faa or --prodigal_dir in the relevant downstream wrapper."
echo "  The directory must follow this structure:"
echo ""
echo "    <prodigal_dir>/"
echo "    +-- <genome_name>/"
echo "        +-- <contig_name>/"
echo "            +-- <contig_name>.faa   (protein sequences)"
echo "            +-- <contig_name>.ffn   (CDS nucleotide sequences)"
echo "            +-- <contig_name>.gff   (gene coordinates, CDS rows)"
echo ""
echo "  .gff requirements:"
echo '    - CDS rows only (gene rows and comment lines are ignored)'
echo '    - ID= attribute must follow gene_id = "{contig}_NNNNN" format'
echo "    - Prodigal GFF files with or without locus_tag are supported"
echo "      ID= is read first, with locus_tag= as fallback"
echo "============================================================"
