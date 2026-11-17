#!/usr/bin/env python3
"""Build the incomplete DefenseFinder genes dictionary.

This retains DefenseFinder hits that appear in the genes dictionary but are not
part of any complete system assembly in the systems dictionary. These genes are
used later for partial-system complementation.

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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes_dict", required=True)
    parser.add_argument("--systems_dict", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def load_pickle(path):
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def main():
    args = parse_args()
    genes_dict = load_pickle(args.genes_dict)
    systems_dict = load_pickle(args.systems_dict)
    incomplete_dict = {}

    for genome_id, contigs in genes_dict.items():
        system_contigs = systems_dict.get(genome_id, {})
        genome_incomplete = {}

        for contig_id, (_, contig_genes) in contigs.items():
            system_gene_ids = set()
            for assemblies in system_contigs.get(contig_id, {}).values():
                for assembly_genes in assemblies:
                    system_gene_ids.update(assembly_genes.keys())

            remaining_genes = {
                hit_id: gene_name
                for hit_id, gene_name in contig_genes.items()
                if hit_id not in system_gene_ids
            }
            if remaining_genes:
                genome_incomplete[contig_id] = (len(remaining_genes), remaining_genes)

        if genome_incomplete:
            incomplete_dict[genome_id] = genome_incomplete

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as handle:
        pickle.dump(incomplete_dict, handle)

    total_genes = sum(
        count for contigs in incomplete_dict.values() for count, _ in contigs.values()
    )
    print(
        f"incomplete genes dict built: {len(incomplete_dict)} genome(s), "
        f"{total_genes} incomplete gene(s)"
    )


if __name__ == "__main__":
    main()
