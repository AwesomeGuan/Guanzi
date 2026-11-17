#!/usr/bin/env python3
"""
Shared DefenseFinder rules loader.

Always reads DefenseFinder_rules_special.tsv first (higher priority),
then DefenseFinder_rules.tsv. Special entries override standard entries
for the same Subsystem key.

Returns one canonical rule dictionary per subsystem. Pattern fields are lists
of pattern groups: each outer item is one component type and the strings in an
inner list are interchangeable profile patterns for that type.

Usage:
  from utils.rules_loader import load_rules
  rules = load_rules("data/DefenseFinder_rules")

  # S05 classification: check gene count
  min_genes = rules["RM_Type_II"]["min_genes"]

  # Assembly scripts: get full patterns
  entry = rules["RM_Type_II"]
  entry["mandatory_patterns"]  # [["RM_Type_II__Type_II_MTase*"], ...]
"""

import csv
from pathlib import Path


def _parse_patterns(value):
    """Parse comma-separated component patterns, including ``{a,b}`` groups."""
    text = "" if value is None else str(value).strip()
    if not text:
        return []

    groups = []
    index = 0
    while index < len(text):
        if text[index] == "{":
            end = index + 1
            depth = 1
            while end < len(text) and depth:
                if text[end] == "{":
                    depth += 1
                elif text[end] == "}":
                    depth -= 1
                end += 1
            variants = [
                item.strip()
                for item in text[index + 1 : end - 1].split(",")
                if item.strip()
            ]
            if variants:
                groups.append(variants)
            index = end
        else:
            end = index
            while end < len(text) and text[end] != ",":
                end += 1
            token = text[index:end].strip()
            if token:
                groups.append([token])
            index = end

        if index < len(text) and text[index] == ",":
            index += 1
        while index < len(text) and text[index].isspace():
            index += 1
    return groups


def _parse_file(path, is_special=False):
    """Parse one rules TSV, return {subsystem: entry_dict}."""
    result = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            subsys = row.get('Subsystem', '').strip()
            if not subsys:
                continue
            try:
                min_mandatory = int(float(row.get('Min_mandatory', 0)))
            except (ValueError, TypeError):
                min_mandatory = 0
            try:
                min_genes = int(float(row.get('Minimum_genes', 0)))
            except (ValueError, TypeError):
                min_genes = 0

            result[subsys] = {
                'system': row.get('System', '').strip(),
                'subsystem': subsys,
                'min_mandatory': min_mandatory,
                'min_genes': min_genes,
                'mandatory_patterns': _parse_patterns(row.get('Mandatory', '')),
                'accessory_patterns': _parse_patterns(row.get('Accessory', '')),
                'forbidden_patterns': _parse_patterns(row.get('Forbidden', '')),
                'neutral_patterns': _parse_patterns(row.get('Neutral', '')),
                'is_special': is_special,
            }
    return result


def load_rules(rules_dir):
    """
    Load DefenseFinder rules from a directory.

    Special rules (DefenseFinder_rules_special.tsv) are loaded first
    and take priority over standard rules (DefenseFinder_rules.tsv)
    for the same Subsystem.

    Args:
        rules_dir: path to directory containing the rules TSV files

    Returns:
        dict: {subsystem: {system, subsystem, min_genes,
               mandatory_patterns, accessory_patterns, ...}}
    """
    base = Path(rules_dir)
    if not base.is_dir():
        print(f"WARNING: rules directory not found: {rules_dir}")
        return {}

    rules = {}

    # Load special rules first (higher priority)
    spc_path = base / "DefenseFinder_rules_special.tsv"
    if spc_path.exists():
        special = _parse_file(spc_path, is_special=True)
        rules.update(special)

    # Load standard rules (special overrides already in place)
    std_path = base / "DefenseFinder_rules.tsv"
    if std_path.exists():
        standard = _parse_file(std_path, is_special=False)
        for key, val in standard.items():
            if key not in rules:
                rules[key] = val

    return rules
