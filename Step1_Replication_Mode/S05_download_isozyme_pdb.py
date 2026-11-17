#!/usr/bin/env python3
"""Download missing AlphaFold PDB files and record every requested outcome."""
import argparse
import csv
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

STATUS_COLUMNS = ["defense_system", "protein_id", "status", "message"]


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--status", required=True)
    parser.add_argument("--threads", type=int, default=8)
    add_progress_argument(parser)
    return parser.parse_args()


def build_url(raw_id, version):
    if "-model_v" in raw_id:
        base = raw_id.split("-model_v", 1)[0]
        return f"https://alphafold.ebi.ac.uk/files/{base}-model_{version}.pdb"
    if raw_id.startswith("AF-"):
        return f"https://alphafold.ebi.ac.uk/files/{raw_id}-model_{version}.pdb"
    return f"https://alphafold.ebi.ac.uk/files/AF-{raw_id}-F1-model_{version}.pdb"


def requested_tasks(input_root):
    tasks = set()
    for txt in Path(input_root).glob("*.txt"):
        for line in txt.read_text(encoding="utf-8").splitlines():
            protein_id = line.strip()
            if protein_id:
                tasks.add((txt.stem, protein_id))
    return sorted(tasks)


def download_one(system, raw_id, output_root):
    target = Path(output_root) / system / f"{raw_id}.pdb"
    if target.is_file() and target.stat().st_size > 0:
        return system, raw_id, "success", "existing file"
    messages = []
    for version in ("v6", "v5", "v4", "v3"):
        url = build_url(raw_id, version)
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200 and response.content:
                target.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile(
                    dir=target.parent, prefix=f".{target.name}.", delete=False
                ) as handle:
                    handle.write(response.content)
                    temporary = Path(handle.name)
                os.replace(temporary, target)
                return system, raw_id, "success", f"downloaded {version}"
            messages.append(f"{version}: HTTP {response.status_code}")
        except requests.RequestException as exc:
            messages.append(f"{version}: {type(exc).__name__}")
    return system, raw_id, "failed", "; ".join(messages)


def write_status(path, rows):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=target.parent,
        prefix=f".{target.name}.", delete=False,
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, target)


def main():
    args = parse_args()
    tasks = requested_tasks(args.input)
    if not tasks:
        raise ValueError(f"No AFDB50 IDs found under {args.input}")
    existing = sum((Path(args.output) / system / f"{protein_id}.pdb").is_file()
                   and (Path(args.output) / system / f"{protein_id}.pdb").stat().st_size > 0
                   for system, protein_id in tasks)
    progress = Progress(args.log_label or "PDB-download", total=len(tasks),
                        processed=0, reused=0, downloaded=0, failed=0, successful=0)
    progress.message(f"existing={existing} pending_downloads={len(tasks) - existing} threads={args.threads}")
    reused = downloaded = failures = 0
    rows = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(download_one, system, protein_id, args.output)
            for system, protein_id in tasks
        ]
        for future in as_completed(futures):
            system, protein_id, status, message = future.result()
            if status == "success":
                if message == "existing file":
                    reused += 1
                else:
                    downloaded += 1
            else:
                failures += 1
                progress.message(f"{system} {protein_id}: {message}", "WARN")
            progress.set(processed=len(rows) + 1, reused=reused, downloaded=downloaded,
                         failed=failures, successful=reused + downloaded)
            rows.append({
                "defense_system": system,
                "protein_id": protein_id,
                "status": status,
                "message": message,
            })
    rows.sort(key=lambda row: (row["defense_system"], row["protein_id"]))
    write_status(args.status, rows)
    progress.finish()
    successful = sum(row["status"] == "success" for row in rows)
    failed = len(rows) - successful
    print(
        f"PDB download pass complete: {successful} successful, {failed} failed; "
        f"status: {args.status}"
    )


if __name__ == "__main__":
    main()
