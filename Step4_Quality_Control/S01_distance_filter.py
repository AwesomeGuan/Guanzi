#!/usr/bin/env python3
"""Keep anchor-connected component genes and recheck system-rule completeness.

Edges join gene-index suffixes differing by at most max_gene_gap; connectivity
may be transitive. Same-role duplicates are reduced before rule evaluation.
"""
import argparse,fnmatch,sys
from pathlib import Path
import pandas as pd
from step4_utils import *
from utils.rules_loader import load_rules

def args():
    p=argparse.ArgumentParser(); p.add_argument("--input",default="output/step3_output"); p.add_argument("--rules",default="data/DefenseFinder_rules"); p.add_argument("--output",default="output/step4_raw/distance_filter"); p.add_argument("--max_gene_gap",type=int,default=10); return p.parse_args()
def pos(pid):
    x=str(pid).rsplit("_",1); return int(x[1]) if len(x)==2 and x[1].isdigit() else None
def match(name,patterns): return any(fnmatch.fnmatch(str(name),p) for p in patterns)
def complete(components,system,rules):
    rule=rules.get(system)
    if not rule:return False
    if any(any(match(c,g) for c in components) for g in rule["forbidden_patterns"]):return False
    mandatory=sum(any(match(c,g) for c in components) for g in rule["mandatory_patterns"]); accessory=sum(any(match(c,g) for c in components) for g in rule["accessory_patterns"])
    return mandatory>=rule["min_mandatory"] and mandatory+accessory>=rule["min_genes"]
def ctype(component,system,rules):
    rule=rules.get(system,{})
    for category in ("mandatory_patterns","accessory_patterns","neutral_patterns"):
        for index,patterns in enumerate(rule.get(category,[])):
            if match(component,patterns):return f"{category}:{index}"
    return str(component)
def clean(row,gap,rules):
    members=[(pid,comp,pos(pid)) for _,pid,comp in iter_assembly_components(row)]
    if not members or members[0][2] is None:return []
    valid=[m for m in members if m[2] is not None]; kept={valid[0][0]}; changed=True
    while changed:
        changed=False
        for pid,comp,p in valid:
            if pid not in kept and any(abs(p-q)<=gap for qid,_,q in valid if qid in kept): kept.add(pid); changed=True
    connected=[m for m in valid if m[0] in kept]; anchor=connected[0]; selected=[anchor]; used={ctype(anchor[1],str(row["defense_system"]),rules)}
    for member in sorted(connected[1:],key=lambda m:(abs(m[2]-anchor[2]),m[2],m[0])):
        kind=ctype(member[1],str(row["defense_system"]),rules)
        if kind not in used:selected.append(member);used.add(kind)
    selected.sort(key=lambda m:(m[2],m[0]))
    if anchor in selected:
        selected.remove(anchor);selected.insert(0,anchor)
    return selected if complete({m[1] for m in selected},str(row["defense_system"]),rules) else []
def rebuild(row,members):
    fixed={c:str(row.get(c,"")) for c in ["defense_system","mining_domain","mining_genome_id","mining_contig_id","mining_contig_system_index"] if c in row.index}
    fixed["mining_component_count"]=str(len(members)); fixed["assembly_type"]=str(row.get("assembly_type",""))
    for i,(pid,comp,_) in enumerate(members,1): fixed[f"mining_component_protein_id_{i}"]=pid; fixed[f"mining_component_{i}"]=comp
    return fixed
def main():
    a=args(); mapping,assembly=load_step3(a.input); rules=load_rules(a.rules); kept=[]; rejected=[]
    for _,row in assembly.iterrows():
        members=clean(row,a.max_gene_gap,rules)
        if members: kept.append(rebuild(row,members))
        else: rejected.append("\t".join(instance_key(row)))
    out_assembly=assembly_records_to_dataframe_local(kept); out_assembly=order_frame(out_assembly,assembly_columns(assembly)); keys=table_keys(out_assembly); out_mapping=subset_keys(mapping,keys)
    member_ids={instance_key(row):set(component_map(row)) for _,row in out_assembly.iterrows()}
    if not out_mapping.empty:
        mask=[str(row["mining_protein_id"]) in member_ids.get(instance_key(row),set()) for _,row in out_mapping.iterrows()]; out_mapping=out_mapping[mask].copy()
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True); order_frame(out_mapping,step3_mapping_columns(out_mapping)).to_csv(out/"step4_distance_filtered_mapping.csv",index=False); order_frame(out_assembly,assembly_columns(out_assembly)).to_csv(out/"step4_distance_filtered_assembly.csv",index=False); (out/"step4_cleaning_log.txt").write_text("\n".join(rejected)+("\n" if rejected else ""),encoding="utf-8")
    print(f"distance filter complete: {len(out_assembly)} assembly instance(s), {len(out_mapping)} mapping row(s), {len(rejected)} rejected")
def assembly_records_to_dataframe_local(records):
    frame=pd.DataFrame(records)
    if frame.empty: frame=pd.DataFrame(columns=["defense_system","mining_genome_id","mining_contig_id","mining_contig_system_index","mining_component_count","assembly_type"])
    return order_frame(frame,assembly_columns(frame))
if __name__=="__main__":main()
