#!/usr/bin/env python3
"""Shared helpers for the four canonical Step3 assembly routes."""

import fnmatch
import sys
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from utils.pipeline_schema import (
    assembly_columns, canonicalize, order_frame, step3_mapping_columns,
)
from utils.rules_loader import load_rules as load_rules_from_dir


def gene_matches_pattern(name, patterns):
    return any(fnmatch.fnmatch(str(name), pattern) for pattern in patterns)


def check_completeness(component_types, defense_system, rules):
    rule = rules.get(defense_system)
    if not rule:
        return False
    if any(any(gene_matches_pattern(c, group) for c in component_types)
           for group in rule["forbidden_patterns"]):
        return False
    mandatory = sum(any(gene_matches_pattern(c, group) for c in component_types)
                    for group in rule["mandatory_patterns"])
    accessory = sum(any(gene_matches_pattern(c, group) for c in component_types)
                    for group in rule["accessory_patterns"])
    return mandatory >= rule["min_mandatory"] and mandatory + accessory >= rule["min_genes"]


def extract_position(protein_id):
    parts = str(protein_id).rsplit("_", 1)
    return int(parts[1]) if len(parts) == 2 and parts[1].isdigit() else None


def sorted_components(primary_protein_id, primary_component, other_genes):
    primary = (str(primary_protein_id), str(primary_component))
    others = [(str(pid), str(comp)) for pid, comp in other_genes.items()
              if str(pid) and str(pid) != primary[0]]
    others.sort(key=lambda item: (
        extract_position(item[0]) is None,
        extract_position(item[0]) if extract_position(item[0]) is not None else 10**18,
        item[0], item[1],
    ))
    return [primary, *others]


def make_assembly_record(system, domain, genome_id, contig_id, instance_index,
                         primary_protein_id, primary_component, other_genes,
                         assembly_type):
    components = sorted_components(primary_protein_id, primary_component, other_genes)
    record = {"defense_system": system}
    if domain != "" and domain is not None:
        record["mining_domain"] = domain
    record.update({
        "mining_genome_id": genome_id,
        "mining_contig_id": contig_id,
        "mining_contig_system_index": str(instance_index),
        "mining_component_count": str(len(components)),
        "assembly_type": assembly_type,
    })
    for index, (protein_id, component) in enumerate(components, start=1):
        record[f"mining_component_protein_id_{index}"] = protein_id
        record[f"mining_component_{index}"] = component
    return record


def assembly_records_to_dataframe(records):
    frame = pd.DataFrame(records)
    if frame.empty:
        frame = pd.DataFrame(columns=[
            "defense_system", "mining_genome_id", "mining_contig_id",
            "mining_contig_system_index", "mining_component_count", "assembly_type",
        ])
    return order_frame(frame, assembly_columns(frame))


def make_mapping_rows(group, instance_index, assembly_type):
    frame = canonicalize(group)
    frame["mining_contig_system_index"] = str(instance_index)
    frame["assembly_type"] = assembly_type
    return order_frame(frame, step3_mapping_columns(frame))


def read_csv_or_empty(path, columns=None):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns or [])
    return canonicalize(pd.read_csv(path, dtype=str, keep_default_na=False))


def write_csv_with_columns(frame, path, columns):
    order_frame(frame, columns).to_csv(Path(path), index=False)


def genes_matching_system(genes, system, rules):
    rule = rules.get(system)
    if not rule:
        return {}
    groups = rule["mandatory_patterns"] + rule["accessory_patterns"] + rule["neutral_patterns"]
    return {gid: name for gid, name in genes.items()
            if any(gene_matches_pattern(name, group) for group in groups)}


def get_complete_assemblies(systems_dict, genome_id, contig_id, system):
    return [dict(item) for item in systems_dict.get(genome_id, {}).get(contig_id, {}).get(system, [])
            if isinstance(item, dict) and item]


def is_duplicate_component(component, existing, system, rules):
    if component in existing:
        return True
    rule = rules.get(system)
    if not rule:
        return False
    groups = rule["mandatory_patterns"] + rule["accessory_patterns"] + rule["neutral_patterns"]
    matched = [group for group in groups if gene_matches_pattern(component, group)]
    return any(gene_matches_pattern(old, group) for group in matched for old in existing)


def next_instance_indices(frames):
    counters = {}
    result = []
    for frame in frames:
        system = str(frame.iloc[0]["defense_system"])
        genome = str(frame.iloc[0]["mining_genome_id"])
        contig = str(frame.iloc[0]["mining_contig_id"])
        key = (system, genome, contig)
        counters[key] = counters.get(key, 0) + 1
        result.append(counters[key])
    return result
