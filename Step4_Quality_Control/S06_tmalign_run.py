#!/usr/bin/env python3
"""Run TM-align for paired mining/query PDB files."""
import argparse
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


METRIC_SUFFIX = ".tmalign.txt"


def parse_metrics(text):
    """Return the four metrics required downstream, or None if incomplete."""
    tm = re.search(r"TM-score=\s*([\d.]+)", text)
    identity = None
    for pattern in (
        r"Seq_ID\s*=\s*n_identical/n_aligned\s*=\s*([\d.]+)",
        r"Seq_ID\s*=\s*([\d.]+)",
        r"(?<![A-Za-z0-9_])ID\s*=\s*([\d.]+)",
    ):
        identity = re.search(pattern, text)
        if identity:
            break
    aligned = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+)", text)
    if not (tm and identity and aligned):
        return None
    return tm.group(1), identity.group(1), aligned.group(1), aligned.group(2)


def metric_path(output_prefix):
    return Path(str(output_prefix) + METRIC_SUFFIX)


def overlay_exists(output_prefix):
    """Support both current .pdb output and legacy extensionless _all_atm output."""
    candidates = [
        Path(str(output_prefix) + ".pdb"),
        Path(str(output_prefix) + "_all_atm"),
    ]
    return any(path.is_file() and path.stat().st_size > 0 for path in candidates)


def write_atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".partial")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mining_assembly_pdb", required=True)
    parser.add_argument("--query_assembly_pdb", required=True)
    parser.add_argument("--tmalign_results", required=True)
    parser.add_argument("--threads", type=int, default=12)
    add_progress_argument(parser)
    return parser.parse_args()


def tmalign_binary():
    return shutil.which("TMalign") or shutil.which("tmalign") or "TMalign"


def worker_count(threads):
    # TM-align is single-threaded, so its concurrency is exactly the requested
    # worker count.
    return max(1, threads)


def run_pair(binary, mining_pdb, query_pdb, output_prefix):
    command = [binary, str(mining_pdb), str(query_pdb), "-a", "T", "-o", str(output_prefix)]
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except OSError as exc:
        return "execution_error", str(exc)
    combined = completed.stdout
    if completed.stderr:
        combined += ("\n" if combined else "") + completed.stderr
    if completed.returncode:
        return "nonzero_exit", f"returncode={completed.returncode}: {combined[-1000:]}"
    if parse_metrics(combined) is None:
        return "invalid_metrics", "TM-align output lacks TM-score, identity, aligned length or RMSD"
    if not overlay_exists(output_prefix):
        return "missing_overlay", "TM-align produced no non-empty .pdb or legacy _all_atm overlay"
    write_atomic(metric_path(output_prefix), combined)
    return "success", ""


def collect_pairs(mining_root, query_root, output_root):
    """Walk the paired trees and return the pairs that still need a run."""
    pairs = []
    missing = reused = migrated = 0
    for system_dir in sorted(mining_root.iterdir()):
        if not system_dir.is_dir():
            continue
        for contig_dir in sorted(system_dir.iterdir()):
            if not contig_dir.is_dir():
                continue
            for mining_dir in sorted(contig_dir.iterdir()):
                if not mining_dir.is_dir():
                    continue
                query_dir = query_root / system_dir.name / contig_dir.name / mining_dir.name
                if not query_dir.exists():
                    missing += 1
                    continue
                for mining_pdb in sorted(mining_dir.glob("*.pdb")):
                    query_pdb = query_dir / mining_pdb.name
                    if not query_pdb.exists():
                        missing += 1
                        continue
                    output_prefix = (
                        output_root
                        / system_dir.name
                        / contig_dir.name
                        / mining_dir.name
                        / mining_pdb.stem
                    )
                    canonical_metrics = metric_path(output_prefix)
                    if (
                        canonical_metrics.is_file()
                        and canonical_metrics.stat().st_size > 0
                        and parse_metrics(canonical_metrics.read_text(encoding="utf-8", errors="ignore"))
                        and overlay_exists(output_prefix)
                    ):
                        reused += 1
                        continue
                    legacy_metrics = Path(str(output_prefix) + "_all_atm")
                    if legacy_metrics.is_file() and legacy_metrics.stat().st_size > 0:
                        legacy_text = legacy_metrics.read_text(encoding="utf-8", errors="ignore")
                        if parse_metrics(legacy_text) is not None:
                            write_atomic(canonical_metrics, legacy_text)
                            reused += 1
                            migrated += 1
                            continue
                    canonical_metrics.unlink(missing_ok=True)
                    output_prefix.parent.mkdir(parents=True, exist_ok=True)
                    pairs.append(
                        (system_dir.name, mining_dir.name, mining_pdb, query_pdb, output_prefix)
                    )
    return pairs, missing, reused, migrated


def main():
    args = parse_args()
    binary = tmalign_binary()
    workers = worker_count(args.threads)
    progress = Progress(
        args.log_label or "TM-align",
        executed=0,
        reused=0,
        nonzero_exit=0,
        invalid_metrics=0,
        missing_overlay=0,
        execution_error=0,
        missing_pairs=0,
        legacy_migrated=0,
        workers=workers,
    )
    pairs, missing, reused, migrated = collect_pairs(
        Path(args.mining_assembly_pdb), Path(args.query_assembly_pdb), Path(args.tmalign_results)
    )
    progress.set(missing_pairs=missing, reused=reused, legacy_migrated=migrated)
    executed = 0
    failures = {
        "nonzero_exit": 0,
        "invalid_metrics": 0,
        "missing_overlay": 0,
        "execution_error": 0,
    }
    if pairs:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(run_pair, binary, mining_pdb, query_pdb, output_prefix): (
                    system,
                    instance,
                    mining_pdb,
                )
                for system, instance, mining_pdb, query_pdb, output_prefix in pairs
            }
            for future in as_completed(futures):
                system, instance, mining_pdb = futures[future]
                status, message = future.result()
                executed += 1
                if status != "success":
                    failures[status] += 1
                    progress.message(f"{mining_pdb}: {status}: {message}", "WARN")
                progress.set(
                    executed=executed,
                    reused=reused,
                    **failures,
                    current_system=system,
                    current_instance=instance,
                    current_component=mining_pdb.stem,
                )
    progress.finish()
    failure_count = sum(failures.values())
    print(
        f"TM-align stage complete: {executed} executed, {reused} reused "
        f"({migrated} legacy metric file(s) migrated), {failure_count} failed"
    )
    if failure_count:
        raise SystemExit(f"TM-align has {failure_count} failed comparison(s)")


if __name__ == "__main__":
    main()
