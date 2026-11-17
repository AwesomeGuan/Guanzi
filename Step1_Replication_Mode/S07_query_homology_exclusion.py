#!/usr/bin/env python3
"""Exclude AFDB50 proteins homologous to their defense-system Query set."""
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query_location", required=True)
    parser.add_argument("--subject_dir", required=True)
    parser.add_argument("--database_dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--marker_dir", required=True)
    add_progress_argument(parser)
    return parser.parse_args()


def write_query_faa(csv_path, target):
    frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    sequences = {}
    for _, row in frame.iterrows():
        protein_id = str(row.get("protein_id", "")).strip()
        sequence = str(row.get("protein_seq", "")).replace("*", "").strip()
        if not protein_id or not sequence:
            continue
        previous = sequences.setdefault(protein_id, sequence)
        if previous != sequence:
            raise ValueError(
                f"{csv_path}: conflicting sequences for Query protein {protein_id}"
            )
    with target.open("w", encoding="utf-8") as handle:
        for protein_id, sequence in sequences.items():
            handle.write(f">{protein_id}\n{sequence}\n")


def ensure_database(system, query_csv, database_root):
    database_root.mkdir(parents=True, exist_ok=True)
    query_faa = database_root / f"{system}.faa"
    database = database_root / f"{system}.dmnd"
    stale = (
        not query_faa.is_file() or not database.is_file()
        or query_faa.stat().st_mtime < query_csv.stat().st_mtime
        or database.stat().st_mtime < query_csv.stat().st_mtime
    )
    if stale:
        write_query_faa(query_csv, query_faa)
        if not query_faa.stat().st_size:
            raise ValueError(f"No Query sequences found in {query_csv}")
        subprocess.run(
            ["diamond", "makedb", "--in", str(query_faa), "--db", str(database)],
            check=True,
        )
    return database, stale


def is_query_homolog(database, subject_faa, output_tsv):
    command = [
        "diamond", "blastp", "-d", str(database), "-q", str(subject_faa),
        "-o", str(output_tsv), "--id", "70", "--query-cover", "70",
        "--subject-cover", "70",
        "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid",
        "pident", "qcovhsp", "--quiet",
    ]
    subprocess.run(command, check=True)
    return output_tsv.is_file() and output_tsv.stat().st_size > 0


def main():
    args = parse_args()
    output_root = Path(args.output)
    marker_root = Path(args.marker_dir)
    database_root = Path(args.database_dir)
    systems = [path for path in sorted(Path(args.subject_dir).iterdir()) if path.is_dir()]
    subject_tasks = {path: sorted(path.glob("*.faa")) for path in systems}
    progress = Progress(args.log_label or "Query-exclusion", total=sum(map(len, subject_tasks.values())),
                        processed=0, reused=0, searched=0, kept=0, removed=0)
    kept = removed = reused = searched = 0

    for system_dir in systems:
        progress.message(f"system={system_dir.name}")
        progress.set(system=system_dir.name)
        query_csv = Path(args.query_location) / f"{system_dir.name}.csv"
        target_dir = output_root / system_dir.name
        target_dir.mkdir(parents=True, exist_ok=True)
        database = None
        rebuilt = False
        if query_csv.is_file():
            database, rebuilt = ensure_database(system_dir.name, query_csv, database_root)
        if rebuilt:
            shutil.rmtree(marker_root / system_dir.name, ignore_errors=True)
        for subject_faa in subject_tasks[system_dir]:
            marker = marker_root / system_dir.name / f"{subject_faa.stem}.ok"
            target = target_dir / subject_faa.name
            if marker.is_file():
                kept += int(target.is_file())
                removed += int(not target.is_file())
                reused += 1
                progress.set(processed=kept + removed, reused=reused, searched=searched, kept=kept, removed=removed)
                continue
            excluded = False
            if database is not None:
                with tempfile.TemporaryDirectory(prefix="rep_query_exclusion_") as temp:
                    searched += 1
                    progress.set(current=subject_faa.stem, searched=searched)
                    excluded = is_query_homolog(
                        database, subject_faa, Path(temp) / "hits.tsv"
                    )
            if excluded:
                target.unlink(missing_ok=True)
                removed += 1
            else:
                shutil.copy2(subject_faa, target)
                kept += 1
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.write_text("ok\n", encoding="utf-8")
            progress.set(processed=kept + removed, reused=reused, searched=searched, kept=kept, removed=removed)
    progress.finish()
    print(f"Query homology exclusion complete: {kept} kept, {removed} removed")


if __name__ == "__main__":
    main()
