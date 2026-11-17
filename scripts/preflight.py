#!/usr/bin/env python3
"""Fail-fast checks for network and NVIDIA GPU stages."""
import argparse
import shutil
import subprocess
import urllib.request
from pathlib import Path


def check_network(args):
    request = urllib.request.Request(args.url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=args.timeout) as response:
            if response.status >= 500:
                raise RuntimeError(f"server returned HTTP {response.status}")
    except Exception as exc:
        raise RuntimeError(f"network endpoint is unavailable: {args.url}: {exc}") from exc
    print(f"Network preflight passed: {args.url}")


def parse_gpu_ids(value):
    fields = [item.strip() for item in value.split(",")]
    if not fields or any(not item.isdigit() for item in fields):
        raise ValueError(f"invalid --gpus value: {value}")
    ids = [int(item) for item in fields]
    if len(ids) != len(set(ids)):
        raise ValueError(f"duplicate GPU ID in --gpus: {value}")
    return ids


def check_gpu(args):
    requested = parse_gpu_ids(args.gpus)
    executable = shutil.which("nvidia-smi")
    if executable is None:
        raise RuntimeError("nvidia-smi is not available")
    result = subprocess.run(
        [executable, "--query-gpu=index", "--format=csv,noheader,nounits"],
        check=True, capture_output=True, text=True,
    )
    available = {int(line.strip()) for line in result.stdout.splitlines() if line.strip()}
    missing = sorted(set(requested) - available)
    if missing:
        raise RuntimeError(f"requested GPU IDs are unavailable: {missing}")
    if args.torch:
        import torch
        if not torch.cuda.is_available():
            raise RuntimeError("PyTorch cannot access CUDA")
    if args.database:
        prefix = Path(args.database)
        if not Path(str(prefix) + ".dbtype").is_file():
            raise RuntimeError(f"FoldSeek database is missing: {prefix}.dbtype")
    print(f"GPU preflight passed: {','.join(map(str, requested))}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="check", required=True)
    network = subparsers.add_parser("network")
    network.add_argument("--url", required=True)
    network.add_argument("--timeout", type=float, default=10.0)
    network.set_defaults(function=check_network)
    gpu = subparsers.add_parser("gpu")
    gpu.add_argument("--gpus", required=True)
    gpu.add_argument("--torch", action="store_true")
    gpu.add_argument("--database")
    gpu.set_defaults(function=check_gpu)
    args = parser.parse_args()
    try:
        args.function(args)
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as exc:
        print(f"PREFLIGHT FAILED: {exc}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
