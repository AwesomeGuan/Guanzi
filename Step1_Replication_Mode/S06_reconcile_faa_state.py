#!/usr/bin/env python3
"""Remove stale or invalid FAA files while retaining valid S06 results."""
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--faa", required=True)
    return parser.parse_args()


def valid_fasta(path):
    header = False
    sequence = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = bool(line[1:].strip())
            elif header:
                sequence = True
    return header and sequence


def main():
    args = parse_args()
    expected = {
        (path.parent.name, path.stem)
        for path in Path(args.pdb).glob("*/*.pdb")
        if path.stat().st_size > 0
    }
    removed = 0
    faa_root = Path(args.faa)
    if faa_root.is_dir():
        for path in faa_root.glob("*/*.faa"):
            if (path.parent.name, path.stem) not in expected or not valid_fasta(path):
                path.unlink()
                removed += 1
    print(f"S06 state reconciled; removed {removed} stale or invalid FAA file(s)")


if __name__ == "__main__":
    main()
