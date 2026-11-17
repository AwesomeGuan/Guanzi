#!/usr/bin/env python3
"""Add genome_id and contig_id to multi-gene mapping CSVs."""
import argparse
from pathlib import Path
import pandas as pd

from step2_candidate_utils import (
    build_contig_genome_index,
    enrich_candidates,
    order_candidate_columns,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--prodigal", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    contig_genome_index = build_contig_genome_index(args.prodigal)
    for csv_path in sorted(Path(args.input).glob("*.csv")):
        df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        system = (
            str(df.iloc[0].get("defense_system", "")).strip()
            if not df.empty else csv_path.stem
        ) or csv_path.stem
        df = enrich_candidates(df, system, contig_genome_index)
        order_candidate_columns(df).to_csv(csv_path, index=False)
    print("mining genome/contig supplement complete")


if __name__ == "__main__":
    main()
