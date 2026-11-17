#!/usr/bin/env python3
"""Extract protein sequences from downloaded AlphaFold PDB files."""
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    add_progress_argument(parser)
    return parser.parse_args()


def main():
    args = parse_args()
    output_root = Path(args.output)
    parser = PDBParser(QUIET=True)
    systems = [path for path in Path(args.input).iterdir() if path.is_dir()]
    pdb_tasks = {path: list(path.glob("*.pdb")) for path in systems}
    progress = Progress(args.log_label or "PDB-to-FAA", total=sum(map(len, pdb_tasks.values())),
                        processed=0, reused=0, written=0, failed=0, successful=0)
    written = 0
    reused = new = failed = 0

    for system_dir in systems:
        for pdb in pdb_tasks[system_dir]:
            target = output_root / system_dir.name / f"{pdb.stem}.faa"
            if target.exists() and target.stat().st_size > 0:
                written += 1
                reused += 1
                progress.set(processed=written + failed, reused=reused, written=new, failed=failed, successful=written)
                continue
            try:
                structure = parser.get_structure(pdb.stem, pdb)
                sequence = "".join(
                    seq1(residue.resname)
                    for residue in structure.get_residues()
                    if residue.id[0] == " "
                )
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(f">{pdb.stem}\n{sequence}\n", encoding="utf-8")
                written += 1
                new += 1
            except Exception as exc:
                failed += 1
                progress.message(f"{system_dir.name} {pdb.stem}: {type(exc).__name__}: {exc}", "WARN")
                progress.set(processed=written + failed, reused=reused, written=new, failed=failed, successful=written)
                continue
            progress.set(processed=written + failed, reused=reused, written=new, failed=failed, successful=written)

    progress.finish()
    print(f"PDB-to-FAA extraction stage complete: {written} file(s)")


if __name__ == "__main__":
    main()
