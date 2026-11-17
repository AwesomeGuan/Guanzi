#!/usr/bin/env python3
"""Build the DefenseFinder genes dictionary from hmmer results.

Input:
    output/step2_raw/defensefinder_results/
        <genome_id>/<contig_id>/<contig_id>_defense_finder_hmmer.tsv

Output:
    output/step3_raw/defensefinder_genes_dict.pkl

Dictionary layout:
    {
        genome_id: {
            contig_id: (
                gene_count,
                {hit_id: gene_name}
            )
        }
    }
"""
import argparse
import pickle
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--defensefinder", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def clean_columns(frame):
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def read_tsv(path):
    return clean_columns(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False))


def extract_contig_from_gene_id(gene_id):
    parts = str(gene_id).rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return str(gene_id)


def main():
    args = parse_args()
    genes_dict = {}
    processed_files = 0
    skipped_files = []

    for hmmer_file in sorted(Path(args.defensefinder).rglob("*_defense_finder_hmmer.tsv")):
        genome_id = hmmer_file.parent.parent.name
        contig_id = hmmer_file.parent.name

        try:
            frame = read_tsv(hmmer_file)
        except Exception as exc:
            skipped_files.append(f"{hmmer_file}: {exc}")
            continue

        if "hit_id" not in frame.columns or "gene_name" not in frame.columns:
            skipped_files.append(f"{hmmer_file}: missing hit_id or gene_name column")
            continue

        contig_genes = {}
        for _, row in frame.iterrows():
            hit_id = str(row["hit_id"]).strip()
            gene_name = str(row["gene_name"]).strip()
            if not hit_id or not gene_name or gene_name.lower() == "nan":
                continue
            target_contig = extract_contig_from_gene_id(hit_id) or contig_id
            contig_genes.setdefault(target_contig, {})
            contig_genes[target_contig].setdefault(hit_id, gene_name)

        genome_entry = genes_dict.setdefault(genome_id, {})
        for target_contig, genes in contig_genes.items():
            genome_entry[target_contig] = (len(genes), genes)
        processed_files += 1

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as handle:
        pickle.dump(genes_dict, handle)

    total_genes = sum(
        count for contigs in genes_dict.values() for count, _ in contigs.values()
    )
    print(
        f"genes dict built: {len(genes_dict)} genome(s), "
        f"{processed_files} hmmer file(s), {total_genes} gene(s)"
    )
    for message in skipped_files[:10]:
        print(f"  skipped: {message}")


if __name__ == "__main__":
    main()
