#!/usr/bin/env python3
"""Canonical table schemas shared by Step1--Step5."""

from pathlib import Path
import re

import pandas as pd


LEGACY_COLUMNS = {
    "query_id": "query_protein_id",
    "mining_id": "mining_protein_id",
    "domain": "mining_domain",
    "genome_id": "mining_genome_id",
    "contig_id": "mining_contig_id",
    "mining_subunit_name": "mining_component",
    "query_genome": "query_genome_id",
    "query_contig": "query_contig_id",
    "query_subunits_count": "query_component_count",
}

QUERY_COLUMNS = [
    "query_genome_id", "query_contig_id", "query_contig_system_index",
    "query_protein_id", "query_component_count",
]
MINING_COLUMNS = [
    "mining_genome_id", "mining_contig_id", "mining_protein_id",
    "mining_component",
]
STEP2_MAPPING_COLUMNS = ["defense_system", *QUERY_COLUMNS]
STEP3_MAPPING_COLUMNS = ["defense_system", *QUERY_COLUMNS]
MINING_INSTANCE_KEY = [
    "defense_system", "mining_genome_id", "mining_contig_id",
    "mining_contig_system_index",
]
QUERY_INSTANCE_KEY = [
    "defense_system", "query_genome_id", "query_contig_id",
    "query_contig_system_index",
]


def canonicalize(frame):
    """Rename legacy columns without overwriting an already-canonical column."""
    frame = frame.copy()
    rename = {
        old: new for old, new in LEGACY_COLUMNS.items()
        if old in frame.columns and new not in frame.columns
    }
    return frame.rename(columns=rename)


def with_optional_domain(columns, position):
    columns = list(columns)
    if "mining_domain" not in columns:
        columns.insert(position, "mining_domain")
    return columns


def step2_columns(frame):
    columns = ["defense_system", *QUERY_COLUMNS]
    if "mining_domain" in frame.columns:
        columns.append("mining_domain")
    columns.extend(MINING_COLUMNS)
    return columns


def step3_mapping_columns(frame):
    columns = ["defense_system", *QUERY_COLUMNS]
    if "mining_domain" in frame.columns:
        columns.append("mining_domain")
    columns.extend([
        "mining_genome_id", "mining_contig_id", "mining_contig_system_index",
        "mining_protein_id", "mining_component", "assembly_type",
    ])
    return columns


def order_frame(frame, columns):
    frame = frame.copy()
    for column in columns:
        if column not in frame.columns:
            frame[column] = ""
    return frame[list(columns)]


def component_indices(columns):
    result = []
    for column in columns:
        match = re.fullmatch(r"mining_component_protein_id_(\d+)", str(column))
        if match:
            result.append(int(match.group(1)))
    return sorted(set(result))


def assembly_columns(frame):
    columns = ["defense_system"]
    if "mining_domain" in frame.columns:
        columns.append("mining_domain")
    columns.extend([
        "mining_genome_id", "mining_contig_id", "mining_contig_system_index",
        "mining_component_count", "assembly_type",
    ])
    for index in component_indices(frame.columns):
        columns.extend([
            f"mining_component_protein_id_{index}",
            f"mining_component_{index}",
        ])
    return columns


def iter_assembly_components(row):
    for index in component_indices(row.index):
        protein_id = str(row.get(f"mining_component_protein_id_{index}", "")).strip()
        component = str(row.get(f"mining_component_{index}", "")).strip()
        if protein_id:
            yield index, protein_id, component


def instance_key(row):
    return tuple(str(row.get(column, "")).strip() for column in MINING_INSTANCE_KEY)


def read_csv(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    return canonicalize(pd.read_csv(path, dtype=str, keep_default_na=False))
