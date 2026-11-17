#!/usr/bin/env python3
"""Import precomputed GTDB DefenseFinder tables for candidate contigs."""
import argparse
from pathlib import Path

import pandas as pd
import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path: sys.path.insert(0, str(PROJECT_ROOT))
from utils.pipeline_schema import canonicalize


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--mapping_dir", required=True)
    parser.add_argument("--prodigal", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_table(path):
    options = {
        "sep": "\t",
        "dtype": str,
        "keep_default_na": False,
    }
    try:
        return pd.read_csv(path, encoding="utf-8", **options)
    except UnicodeDecodeError:
        if hasattr(path, "seek"):
            path.seek(0)
        return pd.read_csv(
            path, encoding="cp1252", encoding_errors="replace", **options
        )


def filter_hit_table(frame, allowed_ids):
    if "hit_id" not in frame.columns:
        return frame.iloc[0:0].copy()
    return frame[frame["hit_id"].isin(allowed_ids)].copy()


def filter_system_table(frame, allowed_ids):
    if "protein_in_syst" not in frame.columns:
        return frame.iloc[0:0].copy()
    rows = []
    for _, row in frame.iterrows():
        source_ids = [value.strip() for value in str(row["protein_in_syst"]).split(",") if value.strip()]
        retained = [value for value in source_ids if value in allowed_ids]
        if not retained:
            continue
        record = row.to_dict()
        record["protein_in_syst"] = ",".join(retained)
        rows.append(record)
    return pd.DataFrame(rows, columns=frame.columns)


def fasta_ids(path):
    identifiers = set()
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.add(line[1:].split()[0])
    return identifiers


def candidate_contigs(mapping_dir):
    rows = []
    for path in sorted(Path(mapping_dir).glob("*.csv")):
        frame = canonicalize(pd.read_csv(path, dtype=str, keep_default_na=False))
        required = {"mining_domain", "mining_genome_id", "mining_contig_id"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        rows.append(frame[["mining_domain", "mining_genome_id", "mining_contig_id"]])
    if not rows:
        raise ValueError(f"No Replication mapping CSV files found under {mapping_dir}")
    return pd.concat(rows, ignore_index=True).drop_duplicates()


def main():
    args = parse_args()
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    written = 0
    for row in candidate_contigs(args.mapping_dir).itertuples(index=False):
        domain, genome, contig = row
        faa = Path(args.prodigal) / genome / contig / f"{contig}.faa"
        if not faa.is_file():
            raise FileNotFoundError(f"Materialized GTDB FAA not found: {faa}")
        allowed_ids = fasta_ids(faa)
        if not allowed_ids:
            raise ValueError(f"No protein IDs found in {faa}")
        source_dir = Path(args.source) / domain / f"{genome}_protein"
        if not source_dir.is_dir():
            raise FileNotFoundError(f"GTDB DefenseFinder directory not found: {source_dir}")
        output_dir = output_root / genome / contig
        output_dir.mkdir(parents=True, exist_ok=True)
        patterns = {
            "hmmer": "*_defense_finder_hmmer.tsv",
            "genes": "*_defense_finder_genes.tsv",
            "systems": "*_defense_finder_systems.tsv",
        }
        current_files = set()
        for kind, pattern in patterns.items():
            matches = sorted(source_dir.glob(pattern))
            if not matches:
                continue
            frame = read_table(matches[0])
            required_column = "protein_in_syst" if kind == "systems" else "hit_id"
            if required_column not in frame.columns:
                raise ValueError(f"DefenseFinder source lacks {required_column}: {matches[0]}")
            filtered = (
                filter_system_table(frame, allowed_ids)
                if kind == "systems" else filter_hit_table(frame, allowed_ids)
            )
            filtered.to_csv(
                output_dir / f"{contig}_defense_finder_{kind}.tsv",
                sep="\t", index=False,
            )
            written += 1
            current_files.add(f"{contig}_defense_finder_{kind}.tsv")
        if not current_files:
            raise ValueError(f"No DefenseFinder tables found: {source_dir}")
        for obsolete in output_dir.glob(f"{contig}_defense_finder_*.tsv"):
            if obsolete.name not in current_files:
                obsolete.unlink()
    print(f"Replication DefenseFinder import complete: {written} table(s), IDs unchanged")


def import_complete(source, mapping_dir, output):
    for domain, genome, contig in candidate_contigs(mapping_dir).itertuples(index=False):
        source_dir = Path(source) / domain / f"{genome}_protein"
        target = Path(output) / genome / contig
        found = False
        for kind, column in (("hmmer", "hit_id"), ("genes", "hit_id"), ("systems", "protein_in_syst")):
            if not list(source_dir.glob(f"*_defense_finder_{kind}.tsv")):
                continue
            found = True
            table = target / f"{contig}_defense_finder_{kind}.tsv"
            if not table.is_file() or column not in read_table(table).columns:
                return False
        if not found:
            return False
    return True


if __name__ == "__main__":
    main()
