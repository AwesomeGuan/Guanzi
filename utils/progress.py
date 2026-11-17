"""Read-only, time-based progress reporting for pipeline tasks."""

import atexit
import sys
import threading
import time


def duration(seconds):
    seconds = max(0, seconds)
    if seconds < 60:
        return f"{seconds:.1f}s"
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{seconds:02d}s"
    if minutes:
        return f"{minutes}m{seconds:02d}s"
    return f"{seconds}s"


def add_progress_argument(parser):
    parser.add_argument("--log_label", help="Optional progress label supplied by the Shell wrapper")


class Progress:
    """Print in-memory snapshots without changing task scheduling or outputs."""

    def __init__(self, label, total=None, interval=300, **fields):
        self.label = label
        self.started = time.monotonic()
        self.fields = ({"total": total} if total is not None else {})
        self.fields.update(fields)
        self.lock = threading.Lock()
        self.stopped = threading.Event()
        self.interval = interval
        self.thread = threading.Thread(target=self._report, daemon=True)
        self.thread.start()
        atexit.register(self.close)

    def set(self, **fields):
        with self.lock:
            self.fields.update(fields)

    def advance(self, **increments):
        with self.lock:
            for name, value in increments.items():
                self.fields[name] = self.fields.get(name, 0) + value

    def message(self, text, event="INFO"):
        stream = sys.stderr if event in ("WARN", "ERROR") else sys.stdout
        print(f"[{self.label}] {event} {text}", file=stream, flush=True)

    def snapshot(self, event="PROGRESS"):
        with self.lock:
            fields = dict(self.fields)
        if "total" in fields and "processed" in fields:
            fields["remaining"] = max(0, fields["total"] - fields["processed"])
        text = " ".join(f"{name}={value}" for name, value in fields.items())
        print(f"[{self.label}] {event} {text} elapsed={duration(time.monotonic() - self.started)}", flush=True)

    def _report(self):
        while not self.stopped.wait(self.interval):
            try:
                self.snapshot()
            except (OSError, ValueError):
                return

    def close(self):
        self.stopped.set()
        if threading.current_thread() is not self.thread:
            self.thread.join(timeout=1)
        atexit.unregister(self.close)

    def finish(self):
        self.close()
        self.snapshot("SUMMARY")
