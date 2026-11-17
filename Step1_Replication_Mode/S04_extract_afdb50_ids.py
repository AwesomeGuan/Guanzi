#!/usr/bin/env python3
"""Extract AFDB50 hit IDs from the final filtered FoldSeek output."""
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    for stale in out.glob("*.txt"):
        stale.unlink()
    for csv_path in sorted(Path(args.input).glob("*.csv")):
        df = pd.read_csv(csv_path)
        if "subject_id" not in df.columns:
            continue
        ids = df["subject_id"].astype(str).unique()
        (out / f"{csv_path.stem}.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")
    print("AFDB50 IDs extracted")


if __name__ == "__main__":
    main()
