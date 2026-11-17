#!/usr/bin/env python3
"""Expand a complete DefenseFinder system with one novel candidate."""
import argparse,pickle
from pathlib import Path
import pandas as pd
from step3_assembly_utils import *

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--rules_dir",required=True); p.add_argument("--systems_dict",required=True); p.add_argument("--output",required=True); return p.parse_args()
def main():
    a=parse_args(); out=Path(a.output); out.mkdir(parents=True,exist_ok=True); rules=load_rules_from_dir(a.rules_dir)
    with Path(a.systems_dict).open("rb") as h: systems=pickle.load(h)
    frame=read_csv_or_empty(a.input); successes=[]; remains=[]; records=[]; counters={}; group_cols=["defense_system"]+(["mining_domain"] if "mining_domain" in frame.columns else [])+["mining_genome_id","mining_contig_id","mining_protein_id","mining_component"]
    for _,group in frame.groupby(group_cols,sort=False):
        row=group.iloc[0]
        system=str(row.get("defense_system","")); genome=str(row.get("mining_genome_id","")); contig=str(row.get("mining_contig_id","")); protein=str(row.get("mining_protein_id","")); component=str(row.get("mining_component","")); chosen=None
        for genes in get_complete_assemblies(systems,genome,contig,system):
            existing=set(genes.values())
            if check_completeness(existing,system,rules) and not is_duplicate_component(component,existing,system,rules): chosen=genes; break
        if chosen is None or not protein or not component: remains.append(group); continue
        key=(system,genome,contig); counters[key]=counters.get(key,0)+1; index=200000+counters[key]; successes.append(make_mapping_rows(group,index,"complete_expansion")); records.append(make_assembly_record(system,str(row.get("mining_domain","")),genome,contig,index,protein,component,chosen,"complete_expansion"))
    success=pd.concat(successes,ignore_index=True) if successes else pd.DataFrame(columns=step3_mapping_columns(frame)); remain=pd.concat(remains,ignore_index=True) if remains else frame.iloc[0:0].copy()
    success.to_csv(out/"mapping.csv",index=False); assembly_records_to_dataframe(records).to_csv(out/"assembly.csv",index=False); remain.to_csv(out/"remain_for_query_guided_assembly.csv",index=False)
    print(f"complete expansion complete: {len(records)} system(s), {len(remain)} remain")
if __name__=="__main__": main()
