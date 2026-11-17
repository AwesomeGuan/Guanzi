#!/usr/bin/env python3
"""Coverage filter for FoldSeek hits."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

COLUMNS = [
    "query_id", "subject_id", "alntmscore", "fident",
    "qcov", "tcov", "evalue", "lddt",
]
NUMERIC_COLUMNS = [
    "alntmscore", "fident", "qcov", "tcov", "evalue", "lddt",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min_qcov", type=float, default=0.70)
    p.add_argument("--min_tcov", type=float, default=0.70)
    return p.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    for stale in output_dir.glob("*.csv"):
        stale.unlink()

    if not 0 <= args.min_qcov <= 1 or not 0 <= args.min_tcov <= 1:
        raise ValueError("--min_qcov and --min_tcov must be between 0 and 1")

    invalid_total = 0
    for system_dir in sorted(path for path in input_dir.iterdir() if path.is_dir()):
        frames = []
        for hit_file in sorted(path for path in system_dir.iterdir() if path.is_file()):
            try:
                frame = pd.read_csv(hit_file, sep="\t", header=None, names=COLUMNS)
            except Exception:
                continue
            if not frame.empty:
                frames.append(frame)
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        for column in NUMERIC_COLUMNS:
            df[column] = pd.to_numeric(df[column], errors="coerce")
        finite = np.isfinite(df[NUMERIC_COLUMNS].to_numpy(dtype=float)).all(axis=1)
        invalid_count = int((~finite).sum())
        invalid_total += invalid_count
        if invalid_count:
            print(
                f"WARNING: skipped {invalid_count} non-finite or non-numeric "
                f"FoldSeek row(s) for {system_dir.name}"
            )
        df = df.loc[finite]
        if df.empty:
            continue
        df = df[(df["qcov"] >= args.min_qcov) & (df["tcov"] >= args.min_tcov)]
        df.to_csv(output_dir / f"{system_dir.name}.csv", index=False)
    print(
        "coverage filter complete: "
        f"qcov >= {args.min_qcov}, tcov >= {args.min_tcov}; "
        f"invalid rows skipped={invalid_total}"
    )


if __name__ == "__main__":
    main()
