#!/usr/bin/env python3
"""Count explicitly selected candidate outputs without changing files or state."""

import argparse
import csv
import math
import sys
from pathlib import Path


def selected_files(paths, pattern):
    files = set()
    for raw in paths:
        path = Path(raw)
        if path.is_file():
            files.add(path)
        elif path.is_dir():
            files.update(p for p in path.glob(pattern) if p.is_file())
        else:
            raise ValueError(f"Candidate output does not exist: {path}")
    return sorted(files)


def csv_count(path, required_columns):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, strict=True)
        if not reader.fieldnames or not any(reader.fieldnames):
            raise ValueError(f"Missing CSV header: {path}")
        missing = set(required_columns) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing candidate columns in {path}: {sorted(missing)}")
        count = 0
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed CSV record: {path}:{reader.line_num}")
            if any(value.strip() for value in row.values()):
                count += 1
        return count


def fasta_count(path):
    count = 0
    have_header = False
    have_sequence = False
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if have_header:
                    if not have_sequence:
                        raise ValueError(f"FASTA record has no sequence: {path}")
                    count += 1
                if not line[1:].strip():
                    raise ValueError(f"Empty FASTA identifier: {path}")
                have_header, have_sequence = True, False
            else:
                if not have_header:
                    raise ValueError(f"FASTA sequence has no header: {path}")
                have_sequence = True
    if have_header:
        if not have_sequence:
            raise ValueError(f"FASTA record has no sequence: {path}")
        count += 1
    return count


def foldseek_count(path):
    count = 0
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 8:
                raise ValueError(f"Expected eight FoldSeek fields: {path}:{number}")
            values = [float(value) for value in fields[2:]]
            if all(math.isfinite(value) for value in values):
                count += 1
    return count


def pdb_success_count(path, pdb_root):
    if pdb_root is None:
        raise ValueError("--pdb-root is required for pdb-status counts")
    count = 0
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, strict=True)
        if reader.fieldnames != ["defense_system", "protein_id", "status", "message"]:
            raise ValueError(f"Unexpected PDB status columns: {path}")
        for row in reader:
            if row["status"] not in {"success", "failed"}:
                raise ValueError(f"Invalid PDB status: {path}:{reader.line_num}")
            if row["status"] == "success":
                pdb = pdb_root / row["defense_system"] / f"{row['protein_id']}.pdb"
                if not pdb.is_file() or not pdb.stat().st_size:
                    raise ValueError(f"Successful PDB is missing or empty: {pdb}")
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("format", choices=("csv", "foldseek", "ids", "fasta", "pdb", "pdb-status"))
    parser.add_argument("--path", action="append", required=True,
                        help="Existing file or directory; may be repeated")
    parser.add_argument("--pattern", help="Explicit directory glob; may include **")
    parser.add_argument("--required-column", action="append", default=[])
    parser.add_argument("--pdb-root", type=Path)
    args = parser.parse_args()
    defaults = {"csv": "*.csv", "foldseek": "*/*", "ids": "*.txt",
                "fasta": "*/*.faa", "pdb": "**/*.pdb", "pdb-status": "*.csv"}
    try:
        files = selected_files(args.path, args.pattern or defaults[args.format])
        if args.format == "ids":
            count = len({(path.stem, line.strip()) for path in files
                         for line in path.read_text(encoding="utf-8").splitlines()
                         if line.strip()})
        elif args.format == "pdb":
            if any(not path.stat().st_size for path in files):
                raise ValueError("Empty PDB file in validated candidate output")
            count = len(files)
        else:
            counter = {"csv": lambda p: csv_count(p, args.required_column),
                       "foldseek": foldseek_count, "fasta": fasta_count,
                       "pdb-status": lambda p: pdb_success_count(p, args.pdb_root)}[args.format]
            count = sum(counter(path) for path in files)
        print(count)
    except (OSError, ValueError, csv.Error, KeyError, TypeError) as exc:
        print(f"CANDIDATE COUNT FAILED: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
