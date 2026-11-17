#!/usr/bin/env python3
"""Reset S05 status while retaining reusable expected PDB files."""
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--status", required=True)
    return parser.parse_args()

def expected_pairs(root):
    pairs = set()
    for path in Path(root).glob("*.txt"):
        pairs.update((path.stem, line.strip())
                     for line in path.read_text(encoding="utf-8").splitlines()
                     if line.strip())
    return pairs

def main():
    args = parse_args()
    expected = expected_pairs(args.input)
    if not expected:
        raise ValueError(f"No AFDB50 IDs found under {args.input}")
    removed = 0
    pdb_root = Path(args.pdb)
    if pdb_root.is_dir():
        for path in pdb_root.glob("*/*.pdb"):
            if (path.parent.name, path.stem) not in expected or path.stat().st_size == 0:
                path.unlink()
                removed += 1
    Path(args.status).unlink(missing_ok=True)
    print(f"S05 state reset; removed {removed} stale or empty PDB file(s)")

if __name__ == "__main__":
    main()
