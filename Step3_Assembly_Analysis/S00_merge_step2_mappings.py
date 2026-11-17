#!/usr/bin/env python3
"""Merge Step2 multi-gene mapping CSVs into one unified mapping table."""
import argparse
from pathlib import Path
import pandas as pd
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from utils.pipeline_schema import canonicalize, order_frame, step2_columns


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output_mapping", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    frames = []
    for csv in sorted(Path(args.input).glob("*.csv")):
        df = canonicalize(pd.read_csv(csv, dtype=str, keep_default_na=False))
        if "defense_system" not in df.columns:
            df.insert(0, "defense_system", csv.stem)
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not merged.empty:
        merged = order_frame(merged, step2_columns(merged))
    Path(args.output_mapping).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output_mapping, index=False)
    print(f"merged {len(frames)} mapping CSVs")


if __name__ == "__main__":
    main()
