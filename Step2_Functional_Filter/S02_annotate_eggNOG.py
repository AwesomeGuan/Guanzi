#!/usr/bin/env python3
"""Run eggNOG-mapper v2 locally and publish the canonical Step2 TSV."""

import argparse
import shutil
import subprocess
from pathlib import Path

from eggnog_annotation_utils import (
    fasta_identifiers, file_digest, publish_annotation, read_annotation_rows,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIXED_OPTIONS = [
    "--itype", "proteins", "-m", "diamond", "--sensmode", "sensitive",
    "--dmnd_iterate", "yes", "--dmnd_ignore_warnings", "--matrix", "BLOSUM62",
    "--gapopen", "11", "--gapextend", "1", "--tax_scope", "auto",
    "--evalue", "0.001", "--target_orthologs", "all",
]


def completed_annotation(path, identifiers):
    header, rows = read_annotation_rows(path, identifiers, require_complete=True)
    with Path(path).open(encoding="utf-8-sig") as handle:
        if not any(line.strip() == "## emapper-2.1.13" for line in handle):
            raise ValueError("Native annotations must identify eggNOG-mapper 2.1.13")
    return header, rows


def outputs_match(native, output, source):
    """A completion marker is usable only with complete, matching tables."""
    header, rows = completed_annotation(native, fasta_identifiers(source))
    published_header, published_rows = read_annotation_rows(output)
    return header == published_header and rows == published_rows


def preflight(data_dir):
    for name in ("eggnog.db", "eggnog.taxa.db", "eggnog_proteins.dmnd"):
        path = data_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"eggNOG resource missing or empty: {path}")
    result = subprocess.run(
        ["emapper.py", "--version", "--data_dir", str(data_dir)],
        check=True, text=True, capture_output=True,
    )
    version = result.stdout + result.stderr
    print(version.strip(), flush=True)
    if "emapper-2.1.13 /" not in version or "Installed eggNOG DB version: 5.0.2 /" not in version:
        raise ValueError("This stage requires eggNOG-mapper 2.1.13 and installed eggNOG DB 5.0.2")
    diamond = subprocess.run(["diamond", "version"], check=True, text=True, capture_output=True).stdout.strip()
    if diamond != "diamond version 2.2.8":
        raise ValueError(f"Expected DIAMOND 2.2.8 in the eggnog environment; found {diamond}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(PROJECT_ROOT / "output/step2_raw/eggnog_upload.faa"))
    parser.add_argument("--data_dir", default=str(PROJECT_ROOT / "data/eggnog_5.0.2"))
    parser.add_argument("--output_dir", default=str(PROJECT_ROOT / "output/step2_raw/eggnog_results"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "output/step2_raw/eggnog_results.tsv"))
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--annotations", help="Explicitly import a completed native v2 annotation for this input")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    source = Path(args.input).resolve()
    identifiers = fasta_identifiers(source)
    if not identifiers:
        print("[Step2/S02] STOP reason=no_candidates input_proteins=0")
        return
    data_dir = Path(args.data_dir).resolve()
    print(f"[Step2/S02] INFO input_proteins={len(identifiers)} threads={args.threads}", flush=True)
    preflight(data_dir)
    raw_dir = Path(args.output_dir).resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)
    native = raw_dir / "query.emapper.annotations"
    command = [
        "emapper.py", "-i", str(source), *FIXED_OPTIONS, "--cpu", str(args.threads),
        "--data_dir", str(data_dir), "--output_dir", str(raw_dir), "-o", "query",
    ]
    reused = False
    if args.annotations:
        imported = Path(args.annotations).resolve()
        header, rows = completed_annotation(imported, identifiers)
        if imported != native:
            if native.exists() and file_digest(native) != file_digest(imported):
                raise ValueError(f"A different native annotation already exists: {native}")
            if not native.exists():
                temporary = native.with_name(f".{native.name}.tmp")
                try:
                    shutil.copyfile(imported, temporary)
                    temporary.replace(native)
                finally:
                    temporary.unlink(missing_ok=True)
        reused = True
    elif native.is_file():
        header, rows = completed_annotation(native, identifiers)
        reused = True
    else:
        if list(raw_dir.glob("query.emapper.*")):
            raise ValueError("Incomplete native output exists; recover or clear this annotation run before retrying")
        subprocess.run(command, check=True, cwd=PROJECT_ROOT)
        header, rows = completed_annotation(native, identifiers)
    changed = publish_annotation(args.output, header, rows)
    print(
        f"[Step2/S02] SUMMARY input_proteins={len(identifiers)} annotated_proteins={len(rows)} "
        f"unannotated_proteins={len(identifiers) - len(rows)} reused={int(reused)} "
        f"tsv_changed={int(changed)} output={args.output}"
    )


if __name__ == "__main__":
    main()
