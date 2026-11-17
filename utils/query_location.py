#!/usr/bin/env python3
"""Shared accessors for the normalized Query-location resource."""

import fnmatch

import pandas as pd


INSTANCE_COLUMNS = ["genome_id", "contig_id", "contig_system_index"]
REQUIRED_COLUMNS = {
    "defense_system",
    "genome_id",
    "contig_id",
    "contig_system_index",
    "component",
    "protein_id",
    "specific_location",
    "protein_nun",
    "raw_instance_frequency",
    "protein_seq",
}


def read_query_location(path, defense_system=None):
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing normalized query columns: {sorted(missing)}")
    if defense_system is not None:
        values = set(frame["defense_system"].astype(str).str.strip())
        if values != {str(defense_system)}:
            raise ValueError(
                f"{path}: defense_system values {sorted(values)} do not match "
                f"file/system name {defense_system}"
            )
    return frame


def instance_mask(frame, genome_id, contig_id, contig_system_index):
    return (
        frame["genome_id"].astype(str).eq(str(genome_id))
        & frame["contig_id"].astype(str).eq(str(contig_id))
        & frame["contig_system_index"].astype(str).eq(str(contig_system_index))
    )


def get_instance(frame, genome_id, contig_id, contig_system_index):
    return frame[
        instance_mask(frame, genome_id, contig_id, contig_system_index)
    ].copy()


def _matches(name, patterns):
    return any(fnmatch.fnmatch(str(name), pattern) for pattern in patterns)


def instance_is_rule_complete(components, defense_system, rules):
    rule = rules.get(defense_system)
    if not rule:
        return False

    components = {str(component) for component in components}
    for forbidden in rule["forbidden_patterns"]:
        if any(_matches(component, forbidden) for component in components):
            return False

    mandatory = sum(
        any(_matches(component, group) for component in components)
        for group in rule["mandatory_patterns"]
    )
    accessory = sum(
        any(_matches(component, group) for component in components)
        for group in rule["accessory_patterns"]
    )
    return (
        mandatory >= rule["min_mandatory"]
        and mandatory + accessory >= rule["min_genes"]
    )


def select_reference_instance(frame, query_id, defense_system, rules=None):
    """Choose one normalized instance containing query_id deterministically.

    Rule-complete instances rank first when rules are supplied, followed by
    distinct component richness, member count, genomic compactness and the
    normalized composite instance key. raw_instance_frequency is deliberately
    excluded because it measures source-dataset redundancy.
    """
    matches = frame[frame["protein_id"].astype(str).eq(str(query_id))]
    if matches.empty:
        return None, None

    candidates = []
    for key in matches[INSTANCE_COLUMNS].drop_duplicates().itertuples(
        index=False, name=None
    ):
        group = get_instance(frame, *key)
        query_rows = group[group["protein_id"].astype(str).eq(str(query_id))]
        if len(query_rows) != 1:
            raise ValueError(
                f"{defense_system}/{query_id}: expected one row in instance "
                f"{key}, found {len(query_rows)}"
            )
        locations = pd.to_numeric(group["specific_location"], errors="raise")
        components = set(group["component"].astype(str))
        complete = (
            instance_is_rule_complete(components, defense_system, rules)
            if rules is not None
            else False
        )
        candidates.append(
            {
                "key": tuple(str(value) for value in key),
                "group": group,
                "query_row": query_rows.iloc[0],
                "complete": complete,
                "component_count": len(components),
                "member_count": len(group),
                "span": int(locations.max() - locations.min()),
            }
        )

    candidates.sort(
        key=lambda item: (
            -int(item["complete"]),
            -item["component_count"],
            -item["member_count"],
            item["span"],
            item["key"],
        )
    )
    selected = candidates[0]
    return selected, len(candidates)
