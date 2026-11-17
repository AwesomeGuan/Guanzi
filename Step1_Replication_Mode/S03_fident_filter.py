#!/usr/bin/env python3
"""Blacklist sequence-similar targets and retain one best Query per target."""
import argparse
from pathlib import Path

import pandas as pd


NUMERIC_COLUMNS = ["alntmscore", "fident", "qcov", "tcov", "evalue", "lddt"]


def filter_and_select_best(frame, blacklist_threshold):
    """Exclude targets above the identity boundary, then keep their best hit."""
    frame = frame.copy()
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["subject_id", "fident", "alntmscore"])
    blacklist = set(
        frame.loc[frame["fident"] > blacklist_threshold, "subject_id"].astype(str)
    )
    frame = frame.loc[~frame["subject_id"].astype(str).isin(blacklist)].copy()
    frame["_input_order"] = range(len(frame))
    frame = frame.sort_values(
        ["alntmscore", "_input_order"],
        ascending=[False, True],
        na_position="last",
        kind="mergesort",
    )
    return frame.drop_duplicates(subset=["subject_id"], keep="first").drop(
        columns="_input_order"
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fident_blacklist_threshold", type=float, default=0.30)
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    for stale in output_dir.glob("*.csv"):
        stale.unlink()
    if not 0 <= args.fident_blacklist_threshold <= 1:
        raise ValueError("--fident_blacklist_threshold must be between 0 and 1")

    for csv_path in sorted(input_dir.glob("*.csv")):
        try:
            frame = pd.read_csv(csv_path)
        except Exception:
            continue
        if frame.empty:
            continue
        raw_rows = len(frame)
        raw_targets = frame["subject_id"].astype(str).nunique()
        frame = filter_and_select_best(frame, args.fident_blacklist_threshold)
        frame.to_csv(output_dir / csv_path.name, index=False)
        print(
            f"{csv_path.stem}: {raw_rows} rows/{raw_targets} targets -> "
            f"{len(frame)} non-homologous unique targets"
        )
    print(
        "fident filter complete: targets with any fident > "
        f"{args.fident_blacklist_threshold} were removed"
    )


if __name__ == "__main__":
    main()
