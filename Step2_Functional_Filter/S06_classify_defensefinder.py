#!/usr/bin/env python3
"""Classify novel candidates into single/multi/complex/cas groups."""
import argparse
from pathlib import Path
import pandas as pd
from step2_candidate_utils import remove_stale_csvs


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--style_file", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.output)
    for d in ("step2_mapping_single", "step2_mapping_multi", "step2_mapping_complex", "step2_mapping_cas"):
        (out_dir / d).mkdir(parents=True, exist_ok=True)
    style_file = Path(args.style_file)
    style = {}
    if not style_file.is_file():
        raise FileNotFoundError(f"Defense-system style table not found: {style_file}")
    sdf = pd.read_csv(style_file)
    required = {"defense_system", "style", "minimum_genes"}
    if not required.issubset(sdf.columns):
        missing = sorted(required.difference(sdf.columns))
        raise ValueError(f"Style table is missing columns: {missing}")
    for _, row in sdf.iterrows():
        style[str(row["defense_system"])] = (
            str(row["style"]),
            int(row.get("minimum_genes", 0) or 0),
        )
    summaries = {k: [] for k in ("single", "multi", "complex", "cas")}
    expected = {k: set() for k in summaries}
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"Novel candidate table not found: {input_path}")
    all_candidates = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    if "defense_system" not in all_candidates.columns:
        raise ValueError("Novel candidate table is missing defense_system")
    for system, df in all_candidates.groupby("defense_system", sort=True):
        system = str(system).strip()
        if not system or df.empty:
            continue
        if system.startswith("CAS_"):
            group = "cas"
        else:
            group = style.get(system, ("multi", 2))[0]
        df.to_csv(out_dir / f"step2_mapping_{group}" / f"{system}.csv", index=False)
        expected[group].add(f"{system}.csv")
        summaries[group].append({"defense_system": system, "candidates": len(df)})
    for group, rows in summaries.items():
        remove_stale_csvs(out_dir / f"step2_mapping_{group}", expected[group])
        pd.DataFrame(rows, columns=["defense_system", "candidates"]).to_csv(
            out_dir / f"step2_summary_{group}.csv", index=False
        )
    print("classification complete")


if __name__ == "__main__":
    main()
