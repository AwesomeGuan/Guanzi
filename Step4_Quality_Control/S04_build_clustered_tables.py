#!/usr/bin/env python3
"""Select CD-HIT assembly-state representatives and preserve formal schemas."""
import argparse,json,shutil
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from step4_utils import *

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument("--cdhit_output",required=True); p.add_argument("--qma_faa",required=True); p.add_argument("--mapping_csv",required=True); p.add_argument("--assembly_csv",required=True); p.add_argument("--output",required=True); return p.parse_args()
def first(path):
    try:return next(SeqIO.parse(path,"fasta"))
    except StopIteration:return None
def main():
    a=parse_args(); qma=Path(a.qma_faa); out=Path(a.output); out.mkdir(parents=True,exist_ok=True); mapping=read_csv(a.mapping_csv); assembly=read_csv(a.assembly_csv); manifest=read_csv(qma/"instance_manifest.csv")
    clusters=json.loads((Path(a.cdhit_output)/"all_clusters.json").read_text(encoding="utf-8")) if (Path(a.cdhit_output)/"all_clusters.json").is_file() else {}
    candidates=[]
    for _,item in manifest.iterrows():
        key=instance_key(item); system=key[0]; context=f"{key[1]}__{key[2]}"; leaf=qma/system/context/key[3]; state=[]; length=0
        for faa in sorted((leaf/"mining_assembly").glob("*.faa")):
            rec=first(faa)
            if rec is None:continue
            length+=len(rec.seq); cid=clusters.get(system,{}).get(faa.stem,{}).get(rec.id,rec.id); state.append(f"{faa.stem}:{cid}")
        candidates.append((system,"|".join(sorted(state)),length,key,leaf))
    reps={}
    for item in candidates:
        group=item[:2]
        if group not in reps or (item[2],item[3])<(reps[group][2],reps[group][3]):reps[group]=item
    keys={item[3] for item in reps.values()}; cm=subset_keys(mapping,keys); ca=subset_keys(assembly,keys)
    order_frame(cm,step3_mapping_columns(cm)).to_csv(out/"clustered_mapping.csv",index=False); order_frame(ca,assembly_columns(ca)).to_csv(out/"clustered_assembly.csv",index=False)
    target=out/"qma_faa_clustered"
    for item in reps.values():
        key=item[3]; destination=target/key[0]/f"{key[1]}__{key[2]}"/key[3]; destination.parent.mkdir(parents=True,exist_ok=True); shutil.copytree(item[4],destination,dirs_exist_ok=True)
    subset_keys(manifest,keys).to_csv(target/"instance_manifest.csv",index=False)
    print(f"CD-HIT representative tables complete: {len(ca)} assembly instance(s), {len(cm)} mapping row(s)")
if __name__=="__main__":main()
