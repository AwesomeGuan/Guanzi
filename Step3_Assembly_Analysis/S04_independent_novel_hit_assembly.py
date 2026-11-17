#!/usr/bin/env python3
"""Assemble complete systems from independent same-contig novel hits."""
import argparse
from pathlib import Path
import pandas as pd
from step3_assembly_utils import *

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--rules_dir",required=True); p.add_argument("--output",required=True); return p.parse_args()

def main():
    a=parse_args(); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    rules=load_rules_from_dir(a.rules_dir); frame=read_csv_or_empty(a.input)
    required={"defense_system","mining_genome_id","mining_contig_id","mining_protein_id","mining_component"}
    if not frame.empty and not required.issubset(frame.columns): raise ValueError(f"missing mapping columns: {sorted(required-set(frame.columns))}")
    group_cols=["defense_system"] + (["mining_domain"] if "mining_domain" in frame.columns else []) + ["mining_genome_id","mining_contig_id"]
    successes=[]; remains=[]; records=[]; counters={}
    for _,group in frame.groupby(group_cols,sort=False):
        if not check_completeness(set(group["mining_component"].astype(str)),str(group.iloc[0]["defense_system"]),rules): remains.append(group); continue
        ordered=group[["mining_protein_id","mining_component"]].drop_duplicates().values.tolist(); ordered.sort(key=lambda x:(extract_position(x[0]) is None,extract_position(x[0]) or 10**18,x[0]))
        primary_id,primary_component=ordered[0]; key=(str(group.iloc[0]["defense_system"]),str(group.iloc[0]["mining_genome_id"]),str(group.iloc[0]["mining_contig_id"])); counters[key]=counters.get(key,0)+1; index=counters[key]
        successes.append(make_mapping_rows(group,index,"independent_novel_hit")); records.append(make_assembly_record(key[0],str(group.iloc[0].get("mining_domain","")),key[1],key[2],index,primary_id,primary_component,dict(ordered[1:]),"independent_novel_hit"))
    success=pd.concat(successes,ignore_index=True) if successes else pd.DataFrame(columns=step3_mapping_columns(frame)); remain=pd.concat(remains,ignore_index=True) if remains else frame.iloc[0:0].copy()
    success.to_csv(out/"mapping.csv",index=False); assembly_records_to_dataframe(records).to_csv(out/"assembly.csv",index=False); remain.to_csv(out/"remain_for_partial_complementation.csv",index=False)
    print(f"independent novel-hit assembly complete: {len(records)} system(s), {len(success)} mapping row(s), {len(remain)} remain")
if __name__=="__main__": main()
