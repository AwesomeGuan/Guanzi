#!/usr/bin/env python3
"""Generate query-to-mining mappings while retaining FoldSeek evidence."""
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--isozyme_dir", required=True)
    p.add_argument("--output_dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.isozyme_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for stale in output_dir.glob("*.csv"):
        stale.unlink()

    for csv_path in sorted(input_dir.glob("*.csv")):
        df = pd.read_csv(csv_path)
        if df.empty:
            continue
        columns = [
            "query_id", "subject_id", "alntmscore", "fident",
            "qcov", "tcov", "evalue", "lddt",
        ]
        missing = [column for column in columns if column not in df.columns]
        if missing:
            raise ValueError(f"{csv_path} is missing FoldSeek columns: {missing}")
        mapping = df[columns].rename(columns={
            "query_id": "query_protein_id",
            "subject_id": "mining_protein_id",
            "alntmscore": "foldseek_alntmscore",
            "fident": "foldseek_sequence_identity",
            "qcov": "foldseek_query_coverage",
            "tcov": "foldseek_target_coverage",
            "evalue": "foldseek_evalue",
            "lddt": "foldseek_lddt",
        })
        mapping.insert(0, "defense_system", csv_path.stem)
        mapping = mapping.drop_duplicates(
            subset=["defense_system", "query_protein_id", "mining_protein_id"]
        )
        mapping.to_csv(output_dir / csv_path.name, index=False)
    print("custom mapping complete")


if __name__ == "__main__":
    main()
