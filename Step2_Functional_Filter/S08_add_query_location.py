#!/usr/bin/env python3
"""Add one exact normalized Query-system context to multi-gene mappings."""

import argparse
import sys
from pathlib import Path

import pandas as pd
from step2_candidate_utils import order_candidate_columns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.query_location import read_query_location, select_reference_instance
from utils.rules_loader import load_rules


CONTEXT_COLUMNS = [
    "mining_component",
    "query_genome_id",
    "query_contig_id",
    "query_contig_system_index",
    "query_component_count",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--query_location", required=True)
    parser.add_argument("--rules_dir", required=True)
    return parser.parse_args()


def load_tables(query_root):
    tables = {}
    for path in sorted(Path(query_root).glob("*.csv")):
        tables[path.stem] = read_query_location(path, path.stem)
    return tables


def main():
    args = parse_args()
    tables = load_tables(args.query_location)
    rules = load_rules(args.rules_dir)
    matched = missing = ambiguous = rule_incomplete = 0

    for csv_path in sorted(Path(args.input).glob("*.csv")):
        frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        if "query_protein_id" not in frame.columns:
            continue
        for column in CONTEXT_COLUMNS:
            frame[column] = ""

        for index, row in frame.iterrows():
            system = str(row.get("defense_system", csv_path.stem)).strip()
            query_id = str(row.get("query_protein_id", "")).strip()
            table = tables.get(system)
            if table is None or not query_id:
                missing += 1
                continue

            selected, candidate_count = select_reference_instance(
                table, query_id, system, rules
            )
            if selected is None:
                missing += 1
                continue
            ambiguous += int(candidate_count > 1)
            rule_incomplete += int(not selected["complete"])
            query_row = selected["query_row"]
            key = selected["key"]
            values = {
                "mining_component": str(query_row["component"]),
                "query_genome_id": key[0],
                "query_contig_id": key[1],
                "query_contig_system_index": key[2],
                "query_component_count": str(selected["component_count"]),
            }
            for column, value in values.items():
                frame.at[index, column] = value
            matched += 1

        order_candidate_columns(frame).to_csv(csv_path, index=False)

    print(
        "query location supplement complete: "
        f"{matched} matched, {missing} missing, "
        f"{ambiguous} selected from multiple unique instances, "
        f"{rule_incomplete} selected instance(s) not rule-complete"
    )


if __name__ == "__main__":
    main()
