#!/usr/bin/env python3
"""Materialize paired Mining/Query component PDB trees from the QMA manifest."""
import argparse,shutil
from pathlib import Path
from Bio import SeqIO
from step4_utils import PROJECT_ROOT, load_rules, pair_component_names

def args():
    p=argparse.ArgumentParser(); p.add_argument("--qma",required=True); p.add_argument("--mining_pdb",required=True); p.add_argument("--query_pdb",required=True); p.add_argument("--mining_output",required=True); p.add_argument("--query_output",required=True); return p.parse_args()
def index(root):
    result={}
    for path in Path(root).rglob("*.pdb"):result.setdefault(path.stem,path)
    return result
def pid(faa):
    try:return next(SeqIO.parse(faa,"fasta")).id
    except StopIteration:return ""
def prepare_pairs(qma, mining_pdbs, query_pdbs, mining_output, query_output, rules):
    counts = dict(mining=0, query=0, rule_matched=0, unmatched_components=0,
                  missing_mining_pdb=0, missing_query_pdb=0)
    qma = Path(qma)
    for mining_dir in sorted(qma.glob("*/*/*/mining_assembly")):
        leaf = mining_dir.parent
        relative = leaf.relative_to(qma)
        mining = {p.stem: p for p in sorted(mining_dir.glob("*.faa"))}
        query = {p.stem: p for p in sorted((leaf / "query_assembly").glob("*.faa"))}
        pairs = pair_component_names(relative.parts[0], mining, query, rules)
        for component, faa in mining.items():
            mining_source = mining_pdbs.get(pid(faa))
            if mining_source is None:
                counts["missing_mining_pdb"] += 1
            else:
                target = Path(mining_output) / relative / f"{component}.pdb"
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(mining_source, target)
                counts["mining"] += 1
            query_component = pairs.get(component)
            if query_component is None:
                counts["unmatched_components"] += 1
                print(f"[S05] WARN no rule-compatible Query component: {relative}/{component}")
                continue
            if query_component != component:
                counts["rule_matched"] += 1
                print(f"[S05] INFO rule pair: {relative}: {component} -> {query_component}")
            query_source = query_pdbs.get(pid(query[query_component]))
            if query_source is None:
                counts["missing_query_pdb"] += 1
                print(f"[S05] WARN missing Query PDB: {relative}/{query_component}")
                continue
            # The paired PDB alias uses the Mining stem; the original Query FAA
            # and its component name remain unchanged for S07 provenance.
            target = Path(query_output) / relative / f"{component}.pdb"
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(query_source, target)
            counts["query"] += 1
    return counts
def main():
    a=args()
    counts=prepare_pairs(a.qma,index(a.mining_pdb),index(a.query_pdb),a.mining_output,a.query_output,
                         load_rules(PROJECT_ROOT/"data"/"DefenseFinder_rules"))
    print("PDB pairs prepared: "+", ".join(f"{name}={count}" for name,count in counts.items()))
if __name__=="__main__":main()
