#!/usr/bin/env python3
"""Build the DefenseFinder complete-systems dictionary.

Input:
    output/step2_raw/defensefinder_results/
        <genome_id>/<contig_id>/<contig_id>_defense_finder_systems.tsv
        <genome_id>/<contig_id>/<contig_id>_defense_finder_genes.tsv

Output:
    output/step3_raw/defensefinder_systems_dict.pkl

Dictionary layout:
    {
        genome_id: {
            contig_id: {
                subtype: [
                    {hit_id: gene_name},
                    ...
                ]
            }
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


def split_values(value):
    return [item.strip() for item in str(value).split(",") if item.strip()]


def extract_contig_from_gene_id(gene_id):
    parts = str(gene_id).rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return str(gene_id)


def build_gene_name_map(genes_file):
    if genes_file is None:
        return {}
    try:
        frame = read_tsv(genes_file)
    except Exception:
        return {}
    if "hit_id" not in frame.columns or "gene_name" not in frame.columns:
        return {}
    gene_map = {}
    for _, row in frame.iterrows():
        hit_id = str(row["hit_id"]).strip()
        gene_name = str(row["gene_name"]).strip()
        if hit_id and gene_name and gene_name.lower() != "nan":
            gene_map.setdefault(hit_id, gene_name)
    return gene_map


def main():
    args = parse_args()
    systems_dict = {}
    processed_files = 0
    skipped_files = []

    for systems_file in sorted(Path(args.defensefinder).rglob("*_defense_finder_systems.tsv")):
        genome_id = systems_file.parent.parent.name
        contig_id = systems_file.parent.name
        genes_file = systems_file.with_name(
            systems_file.name.replace("_defense_finder_systems.tsv", "_defense_finder_genes.tsv")
        )
        if not genes_file.exists():
            genes_file = None

        try:
            systems_frame = read_tsv(systems_file)
        except Exception as exc:
            skipped_files.append(f"{systems_file}: {exc}")
            continue

        required = {"sys_id", "subtype", "sys_beg", "protein_in_syst"}
        if not required.issubset(systems_frame.columns):
            missing = sorted(required.difference(systems_frame.columns))
            skipped_files.append(f"{systems_file}: missing columns {missing}")
            continue

        gene_name_map = build_gene_name_map(genes_file)
        genome_entry = systems_dict.setdefault(genome_id, {})

        for _, row in systems_frame.iterrows():
            sys_id = str(row["sys_id"]).strip()
            subtype = str(row["subtype"]).strip()
            sys_beg = str(row["sys_beg"]).strip()
            protein_ids = split_values(row["protein_in_syst"])
            if not sys_id or not subtype or not sys_beg or not protein_ids:
                continue

            target_contig = extract_contig_from_gene_id(sys_beg) or contig_id
            assembly_genes = {}

            profile_names = split_values(row.get("name_of_profiles_in_sys", ""))
            for index, protein_id in enumerate(protein_ids):
                if index < len(profile_names) and profile_names[index]:
                    assembly_genes[protein_id] = profile_names[index]

            for protein_id in protein_ids:
                if protein_id in assembly_genes:
                    continue
                if protein_id in gene_name_map:
                    assembly_genes[protein_id] = gene_name_map[protein_id]
                else:
                    assembly_genes[protein_id] = protein_id

            contig_entry = genome_entry.setdefault(target_contig, {})
            contig_entry.setdefault(subtype, []).append(assembly_genes)

        processed_files += 1

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as handle:
        pickle.dump(systems_dict, handle)

    total_assemblies = sum(
        len(assemblies)
        for contigs in systems_dict.values()
        for subtypes in contigs.values()
        for assemblies in subtypes.values()
    )
    print(
        f"systems dict built: {len(systems_dict)} genome(s), "
        f"{processed_files} systems file(s), {total_assemblies} assembly record(s)"
    )
    for message in skipped_files[:10]:
        print(f"  skipped: {message}")


if __name__ == "__main__":
    main()
