#!/usr/bin/env python3
"""Run a pipeline command with consistent resource and elapsed-time messages."""

import argparse
import signal
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.progress import duration


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--mode", choices=("custom", "replication", "unspecified"), default="unspecified")
    parser.add_argument("--device", choices=("CPU", "GPU", "mixed"), default="CPU")
    parser.add_argument("--gpus")
    parser.add_argument("--threads")
    parser.add_argument("--processes-per-gpu")
    parser.add_argument("--info", action="append", default=[])
    parser.add_argument("--counted", action="store_true", help="The child reports its own 300-second snapshots")
    parser.add_argument("--event", choices=("START", "END", "SKIP"))
    parser.add_argument("--started-at", type=float)
    parser.add_argument("--exit-code", type=int, default=0)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return parser.parse_args()


def main():
    args = parse_args()
    command = args.command
    if command[:1] == ["--"]:
        command = command[1:]
    resources = f"device={args.device}"
    if args.gpus:
        resources += f" gpus={args.gpus}"
    if args.threads:
        resources += f" threads={args.threads}"
    if args.processes_per_gpu:
        resources += f" processes_per_gpu={args.processes_per_gpu}"
    prefix = f"[{args.label}]"
    if args.event:
        elapsed = duration(time.monotonic() - args.started_at) if args.started_at is not None else "0s"
        event = args.event
        if event == "END" and args.exit_code:
            event = "PAUSE" if args.label == "Step2" and args.exit_code == 3 else "ERROR"
        detail = " execution=skipped output_check=passed" if event == "SKIP" else ""
        suffix = "" if event == "START" else f" elapsed={elapsed}{detail}"
        if args.event == "END":
            suffix += f" exit_code={args.exit_code}"
        print(f"{prefix} {event} mode={args.mode} {resources}{suffix}",
              file=sys.stderr if event == "ERROR" else sys.stdout, flush=True)
        return 0
    if not command:
        raise SystemExit("A command is required after --")
    print(f"{prefix} START mode={args.mode} {resources}", flush=True)
    for information in args.info:
        print(f"{prefix} INFO {information}", flush=True)
    if args.counted:
        command += ["--log_label", args.label]
    started = time.monotonic()
    child = None
    old_handlers = {}

    def forward(signum, frame):
        if child is not None and child.poll() is None:
            child.send_signal(signum)

    try:
        # Inherit both streams: tool output remains visible and is never captured.
        child = subprocess.Popen(command)
        for signum in (signal.SIGINT, signal.SIGTERM):
            old_handlers[signum] = signal.signal(signum, forward)
        while True:
            try:
                code = child.wait(timeout=300)
                break
            except subprocess.TimeoutExpired:
                if not args.counted:
                    print(f"{prefix} RUNNING {resources} elapsed={duration(time.monotonic() - started)}", flush=True)
    except OSError as exc:
        print(f"{prefix} ERROR {exc} {resources} elapsed={duration(time.monotonic() - started)}", file=sys.stderr, flush=True)
        return 1
    finally:
        for signum, handler in old_handlers.items():
            signal.signal(signum, handler)
    event = "END" if code == 0 else ("PAUSE" if args.label == "Step2" and code == 3 else "ERROR")
    print(f"{prefix} {event} {resources} elapsed={duration(time.monotonic() - started)} exit_code={code}",
          file=sys.stderr if event == "ERROR" else sys.stdout, flush=True)
    return code if code >= 0 else 128 - code


if __name__ == "__main__":
    sys.exit(main())
