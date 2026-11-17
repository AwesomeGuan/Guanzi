#!/usr/bin/env python3
"""Complement incomplete DefenseFinder systems with novel candidates."""
import argparse,pickle
from pathlib import Path
import pandas as pd
from step3_assembly_utils import *

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--rules_dir",required=True); p.add_argument("--incomplete_dict",required=True); p.add_argument("--output",required=True); return p.parse_args()
def main():
    a=parse_args(); out=Path(a.output); out.mkdir(parents=True,exist_ok=True); rules=load_rules_from_dir(a.rules_dir)
    with Path(a.incomplete_dict).open("rb") as h: incomplete=pickle.load(h)
    frame=read_csv_or_empty(a.input); successes=[]; remains=[]; records=[]; counters={}; group_cols=["defense_system"]+(["mining_domain"] if "mining_domain" in frame.columns else [])+["mining_genome_id","mining_contig_id","mining_protein_id","mining_component"]
    for _,group in frame.groupby(group_cols,sort=False):
        row=group.iloc[0]
        system=str(row.get("defense_system","")); genome=str(row.get("mining_genome_id","")); contig=str(row.get("mining_contig_id","")); protein=str(row.get("mining_protein_id","")); component=str(row.get("mining_component","")); genes=genes_matching_system(incomplete.get(genome,{}).get(contig,(0,{}))[1],system,rules)
        if not protein or not component or check_completeness(set(genes.values()),system,rules) or not check_completeness(set(genes.values())|{component},system,rules): remains.append(group); continue
        key=(system,genome,contig); counters[key]=counters.get(key,0)+1; index=100000+counters[key]; successes.append(make_mapping_rows(group,index,"partial_complementation")); records.append(make_assembly_record(system,str(row.get("mining_domain","")),genome,contig,index,protein,component,genes,"partial_complementation"))
    success=pd.concat(successes,ignore_index=True) if successes else pd.DataFrame(columns=step3_mapping_columns(frame)); remain=pd.concat(remains,ignore_index=True) if remains else frame.iloc[0:0].copy()
    success.to_csv(out/"mapping.csv",index=False); assembly_records_to_dataframe(records).to_csv(out/"assembly.csv",index=False); remain.to_csv(out/"remain_for_complete_expansion.csv",index=False)
    print(f"partial complementation complete: {len(records)} system(s), {len(remain)} remain")
if __name__=="__main__": main()
