#!/usr/bin/env python3
"""Shared Step4 instance and component helpers."""
import fnmatch, re, sys
from pathlib import Path
import pandas as pd

PROJECT_ROOT=Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path: sys.path.insert(0,str(PROJECT_ROOT))
from utils.pipeline_schema import *
from utils.rules_loader import load_rules

ROUTES=[
    "independent_novel_hit_assembly","partial_complementation",
    "complete_expansion","query_guided_assembly",
]


def component_name_pattern(pattern):
    """Match rule patterns against the sanitized QMA filenames."""
    return re.sub(r"[^A-Za-z0-9_.*?\[\]-]+", "_", str(pattern))


def name_distance(left, right):
    """Use name distance only to order alternatives within one rule role."""
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        current = [i]
        for j, b in enumerate(right, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (a != b)))
        previous = current
    return previous[-1]


def pair_component_names(system, mining_names, query_names, rules):
    """Pair exactly first, then one-to-one within rule-defined roles.

    This follows the historical Code/601TM-align.py policy. Name distance is
    not a sequence/structure threshold and never permits cross-role matching.
    Original Mining and Query names remain separate in the returned mapping.
    """
    mining = set(mining_names)
    query = set(query_names)
    pairs = {name: name for name in sorted(mining & query)}
    mining.difference_update(pairs)
    query.difference_update(pairs.values())
    rule = rules.get(system, {})
    groups = [group for category in ("mandatory_patterns", "accessory_patterns", "neutral_patterns")
              for group in rule.get(category, [])]
    roles = lambda name: {i for i, group in enumerate(groups)
                          if any(fnmatch.fnmatchcase(name, component_name_pattern(p)) for p in group)}
    mining_roles = {name: roles(name) for name in mining}
    query_roles = {name: roles(name) for name in query}
    alternatives = sorted(
        (name_distance(m, q), m, q) for m in mining for q in query
        if mining_roles[m] & query_roles[q]
    )
    for _, m, q in alternatives:
        if m in mining and q in query:
            pairs[m] = q
            mining.remove(m)
            query.remove(q)
    return pairs

def load_step3(root):
    mappings=[]; assemblies=[]
    for route in ROUTES:
        directory=Path(root)/route
        m=read_csv(directory/"mapping.csv"); a=read_csv(directory/"assembly.csv")
        if not m.empty: mappings.append(m)
        if not a.empty: assemblies.append(a)
    mapping=pd.concat(mappings,ignore_index=True) if mappings else pd.DataFrame()
    assembly=pd.concat(assemblies,ignore_index=True) if assemblies else pd.DataFrame()
    return mapping,assembly

def key_mask(frame,key):
    if frame.empty: return pd.Series([],dtype=bool)
    mask=pd.Series(True,index=frame.index)
    for column,value in zip(MINING_INSTANCE_KEY,key):
        mask &= frame[column].astype(str).eq(str(value))
    return mask

def slug(row):
    values=[str(row.get(c,"")) for c in MINING_INSTANCE_KEY]
    return "__".join(re.sub(r"[^A-Za-z0-9_.-]+","_",v) for v in values)

def component_map(row):
    return {pid:comp for _,pid,comp in iter_assembly_components(row)}

def primary_component(row):
    return (
        str(row.get("mining_component_protein_id_1","")).strip(),
        str(row.get("mining_component_1","")).strip(),
    )

def table_keys(frame):
    if frame.empty: return set()
    return set(map(tuple,frame[MINING_INSTANCE_KEY].astype(str).itertuples(index=False,name=None)))

def subset_keys(frame,keys):
    if frame.empty: return frame.copy()
    values=list(map(tuple,frame[MINING_INSTANCE_KEY].astype(str).itertuples(index=False,name=None)))
    return frame[[value in keys for value in values]].copy()

def write_pair(mapping,assembly,directory,mapping_name="mapping.csv",assembly_name="assembly.csv"):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    order_frame(mapping,step3_mapping_columns(mapping)).to_csv(directory/mapping_name,index=False)
    order_frame(assembly,assembly_columns(assembly)).to_csv(directory/assembly_name,index=False)
