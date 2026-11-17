#!/usr/bin/env python3
"""Build Replication mapping without grafting AFDB50 FoldSeek metrics."""
import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


FINAL_COLUMNS = [
    "defense_system", "query_protein_id", "afdb50_id", "mining_domain",
    "mining_genome_id", "mining_contig_id", "mining_protein_id",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtdb_mapping", required=True)
    parser.add_argument("--foldseek_results", required=True)
    parser.add_argument("--prodigal", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def load_gtdb_mapping(mapping_root):
    frames = []
    for path in sorted(Path(mapping_root).glob("*/*.csv")):
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        raise ValueError("No GTDB mappings were found")
    mapping = pd.concat(frames, ignore_index=True)
    required = {
        "defense_system", "domain", "genome_id", "contig_id",
        "afdb50_id", "gtdb_protein_id",
    }
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"GTDB mapping is missing columns: {sorted(missing)}")
    return mapping


def build_protein_index(prodigal_root):
    index = set()
    for faa in Path(prodigal_root).rglob("*.faa"):
        genome_id = faa.parent.parent.name
        contig_id = faa.parent.name
        with faa.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    protein_id = line[1:].split()[0]
                    if protein_id:
                        index.add((genome_id, contig_id, protein_id))
    return index


def load_query_links(path):
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"query_id", "subject_id"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing FoldSeek columns: {sorted(missing)}")
    lookup = defaultdict(set)
    for _, row in frame.iterrows():
        lookup[str(row["subject_id"])].add(str(row["query_id"]))
    return lookup


def main():
    args = parse_args()
    gtdb_mapping = load_gtdb_mapping(args.gtdb_mapping)
    protein_index = build_protein_index(args.prodigal)
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    for stale in output_root.glob("*.csv"):
        stale.unlink()
    written = 0
    skipped_unmaterialized = 0

    for system, system_rows in gtdb_mapping.groupby("defense_system", sort=True):
        evidence_path = Path(args.foldseek_results) / f"{system}.csv"
        query_links = load_query_links(evidence_path)
        rows = []
        for _, mapping_row in system_rows.iterrows():
            genome_id = str(mapping_row["genome_id"])
            contig_id = str(mapping_row["contig_id"])
            mining_id = str(mapping_row["gtdb_protein_id"])
            if (genome_id, contig_id, mining_id) not in protein_index:
                skipped_unmaterialized += 1
                continue
            afdb50_id = str(mapping_row["afdb50_id"])
            for query_id in sorted(query_links.get(afdb50_id, set())):
                record = {
                    "defense_system": system,
                    "query_protein_id": query_id,
                    "afdb50_id": afdb50_id,
                    "mining_domain": str(mapping_row["domain"]),
                    "mining_genome_id": genome_id,
                    "mining_contig_id": contig_id,
                    "mining_protein_id": mining_id,
                }
                rows.append(record)
        result = pd.DataFrame(rows, columns=FINAL_COLUMNS)
        output_path = output_root / f"{system}.csv"
        if result.empty:
            output_path.unlink(missing_ok=True)
            continue
        result = result.drop_duplicates(
            FINAL_COLUMNS
        )
        result.to_csv(output_path, index=False)
        written += len(result)
    print(
        f"Replication mapping complete: {written} row(s), "
        f"{skipped_unmaterialized} unmaterialized GTDB mapping(s) skipped"
    )


if __name__ == "__main__":
    main()
