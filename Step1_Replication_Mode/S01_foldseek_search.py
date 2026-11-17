#!/usr/bin/env python3
"""FoldSeek structural search with auto-retry."""

import os
import shutil
import subprocess
import sys
import tempfile
import time
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--query_pdb", required=True)
    p.add_argument("--db_prefix", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--num_threads", type=int, default=12)
    p.add_argument("--max_missing", type=int, default=100)
    p.add_argument("--max_rounds", type=int, default=5)
    add_progress_argument(p)
    return p.parse_args()


def count_files(folder):
    if not os.path.exists(folder):
        return 0
    return len([f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))])


def remove_partial_output(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


def process_one(pdb_file, input_dir, output_dir, db_pad):
    input_path = os.path.join(input_dir, pdb_file)
    out_prefix = os.path.join(output_dir, os.path.splitext(pdb_file)[0])
    if os.path.exists(out_prefix):
        return "skip", pdb_file, "", 0.0

    tmp_dir = tempfile.mkdtemp(prefix="foldseek_easy_search_")
    started = time.monotonic()
    try:
        cmd = [
            "foldseek", "easy-search", input_path, db_pad, out_prefix, tmp_dir,
            "--format-output", "query,target,alntmscore,fident,qcov,tcov,evalue,lddt",
            "--tmscore-threshold", "0.7", "--lddt-threshold", "0.7",
            "--gpu", "1", "--prefilter-mode", "1",
        ]
        subprocess.run(
            cmd,
            check=True,
            timeout=600,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        if not os.path.exists(out_prefix):
            return "fail", pdb_file, "FoldSeek exited successfully but created no output file", time.monotonic() - started
        return "ok", pdb_file, "", time.monotonic() - started
    except subprocess.TimeoutExpired as exc:
        remove_partial_output(out_prefix)
        detail = (exc.stderr or "").strip()
        message = "timeout after 600s"
        if detail:
            message += f": {detail[-1000:]}"
        return "fail", pdb_file, message, time.monotonic() - started
    except subprocess.CalledProcessError as exc:
        remove_partial_output(out_prefix)
        detail = (exc.stderr or "").strip()
        message = f"exit code {exc.returncode}"
        if detail:
            message += f": {detail[-1000:]}"
        return "fail", pdb_file, message, time.monotonic() - started
    except Exception as exc:
        remove_partial_output(out_prefix)
        return "fail", pdb_file, str(exc), time.monotonic() - started
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def run_one_round(systems, query_dir, out_dir, db_pad, num_threads, log_label=None):
    skip_total = ok_total = fail_total = 0
    worker_count = max(1, num_threads)

    for sys_name in sorted(systems):
        sys_in = os.path.join(query_dir, sys_name)
        sys_out = os.path.join(out_dir, sys_name)
        os.makedirs(sys_out, exist_ok=True)
        pdbs = [f for f in os.listdir(sys_in) if f.endswith(".pdb")]
        todo = [f for f in pdbs if not os.path.exists(os.path.join(sys_out, os.path.splitext(f)[0]))]
        skipped = len(pdbs) - len(todo)
        results = {"ok": 0, "fail": 0}
        if not todo:
            print(f"  {sys_name:<40} skipped={skipped} ok=0 fail=0")
            skip_total += skipped
            continue
        progress = Progress(log_label or "FoldSeek", total=len(pdbs), system=sys_name,
                            processed=skipped, reused=skipped, searched=0, failed=0,
                            successful=skipped)
        completed = 0
        print(
            f"  {sys_name:<40} skipped={skipped} "
            f"progress={completed}/{len(todo)} ok=0 fail=0",
            end="",
            flush=True,
        )
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(process_one, pdb_file, sys_in, sys_out, db_pad): pdb_file
                for pdb_file in todo
            }
            for future in as_completed(futures):
                status, pdb_file, error, elapsed = future.result()
                completed += 1
                if status == "ok":
                    results["ok"] += 1
                elif status == "fail":
                    results["fail"] += 1
                    print(
                        f"\n    FAIL {pdb_file} ({elapsed:.1f}s): {error}",
                        file=sys.stderr,
                        flush=True,
                    )
                progress.set(processed=skipped + completed, searched=results["ok"],
                             failed=results["fail"], successful=skipped + results["ok"])
                print(
                    f"\r  {sys_name:<40} skipped={skipped} "
                    f"progress={completed}/{len(todo)} "
                    f"ok={results['ok']} fail={results['fail']}",
                    end="",
                    flush=True,
                )
        print()
        progress.finish()
        skip_total += skipped
        ok_total += results["ok"]
        fail_total += results["fail"]
    return skip_total, ok_total, fail_total


def get_stats(query_dir, out_dir, systems):
    grand_in = grand_out = 0
    per_sys = []
    for sys_name in sorted(systems):
        in_c = len([f for f in os.listdir(os.path.join(query_dir, sys_name)) if f.endswith(".pdb")])
        out_c = count_files(os.path.join(out_dir, sys_name))
        grand_in += in_c
        grand_out += out_c
        per_sys.append((sys_name, in_c, out_c))
    return grand_in, grand_out, grand_in - grand_out, per_sys


def main():
    args = parse_args()
    query_dir = args.query_pdb
    out_dir = args.output
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    db_pad = os.path.join(project_root, "foldseek", args.db_prefix + "_pad")
    systems = sorted(d for d in os.listdir(query_dir) if os.path.isdir(os.path.join(query_dir, d)))
    os.makedirs(out_dir, exist_ok=True)

    grand_in, grand_out, grand_missing, per_sys = get_stats(query_dir, out_dir, systems)
    print(f"Systems: {len(systems)}, db={db_pad}")
    print(f"Total: {grand_in} PDBs, done: {grand_out}, missing: {grand_missing}")
    print(
        f"Visible GPUs: {os.environ.get('CUDA_VISIBLE_DEVICES', 'default')}, "
        f"parallel FoldSeek processes: {max(1, args.num_threads)}"
    )

    round_num = 0
    prev_missing = grand_missing
    consecutive_small = 0
    while True:
        round_num += 1
        skip_c, ok_c, fail_c = run_one_round(
            systems, query_dir, out_dir, db_pad, args.num_threads, args.log_label
        )
        _, _, grand_missing, per_sys = get_stats(query_dir, out_dir, systems)
        diff = prev_missing - grand_missing
        prev_missing = grand_missing
        print(f"Round {round_num}: skip={skip_c} ok={ok_c} fail={fail_c} missing={grand_missing} diff={diff}")
        if grand_missing == 0:
            print("All done!")
            break
        if diff <= 10:
            consecutive_small += 1
        else:
            consecutive_small = 0
        if grand_missing <= args.max_missing and diff <= 10:
            break
        if consecutive_small >= 3:
            break
        if args.max_rounds > 0 and round_num >= args.max_rounds:
            break

    _, _, grand_missing, per_sys = get_stats(query_dir, out_dir, systems)
    for sys_name, in_c, out_c in per_sys:
        print(f"{sys_name:<30} in={in_c} out={out_c}")
    print(f"Done: {grand_missing} missing")
    if grand_missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
