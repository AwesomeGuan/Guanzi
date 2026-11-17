#!/usr/bin/env python3
"""Build the canonical QMA FAA tree and its instance manifest."""
import argparse,re,sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from step4_utils import *
PROJECT_ROOT=Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:sys.path.insert(0,str(PROJECT_ROOT))
from utils.query_location import get_instance,read_query_location

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument("--prodigal_dir",required=True); p.add_argument("--assembly_csv",required=True); p.add_argument("--mapping_csv",required=True); p.add_argument("--query_location",required=True); p.add_argument("--qma_faa",required=True); return p.parse_args()
def safe(v):return re.sub(r"[^A-Za-z0-9_.-]+","_",str(v))
def fasta_index(root):
    result={}
    for path in Path(root).rglob("*.faa"):
        for rec in SeqIO.parse(path,"fasta"):result.setdefault(rec.id,str(rec.seq).replace("*",""))
    return result
def write(path,pid,seq):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(f">{pid}\n{seq}\n",encoding="utf-8")
def query_components(instance, query_protein_id):
    own=instance[instance["protein_id"].astype(str).eq(str(query_protein_id))]
    if own.empty:return instance.iloc[0:0]
    anchor=own.iloc[0]; location=pd.to_numeric(anchor.get("specific_location",""),errors="coerce"); selected=[]
    for component,group in instance.groupby("component",sort=False):
        exact=group[group["protein_id"].astype(str).eq(str(query_protein_id))]
        if not exact.empty: selected.append(exact.iloc[0]);continue
        ranked=group.copy(); ranked["_distance"]=(pd.to_numeric(ranked["specific_location"],errors="coerce")-location).abs(); ranked=ranked.sort_values(["_distance","specific_location","protein_id"],kind="stable");selected.append(ranked.iloc[0])
    return pd.DataFrame(selected)
def main():
    a=parse_args(); assembly=read_csv(a.assembly_csv); mapping=read_csv(a.mapping_csv); proteins=fasta_index(a.prodigal_dir); out=Path(a.qma_faa); out.mkdir(parents=True,exist_ok=True)
    tables={p.stem:read_query_location(p,p.stem) for p in Path(a.query_location).glob("*.csv")}; manifest=[]; missing=0
    for _,row in assembly.iterrows():
        key=instance_key(row); members=list(iter_assembly_components(row)); primary=members[0][1] if members else ""; related=mapping[key_mask(mapping,key)].copy()
        preferred=related[related["mining_protein_id"].astype(str).eq(primary)] if not related.empty else related
        selected=(preferred if not preferred.empty else related).iloc[0] if not related.empty else None
        context=safe(f"{row.get('mining_genome_id','')}__{row.get('mining_contig_id','')}"); leaf=out/safe(row.get("defense_system",""))/context/safe(row.get("mining_contig_system_index","")); mining_dir=leaf/"mining_assembly"; query_dir=leaf/"query_assembly"
        for _,pid,component in members:
            seq=proteins.get(pid,"")
            if seq:write(mining_dir/f"{safe(component)}.faa",pid,seq)
            else:missing+=1
        record={c:str(row.get(c,"")) for c in MINING_INSTANCE_KEY}; record["instance_path_context"]=context; record["primary_mining_protein_id"]=primary
        if selected is not None:
            for c in QUERY_INSTANCE_KEY[1:]+["query_protein_id"]:record[c]=str(selected.get(c,""))
            table=tables.get(str(row.get("defense_system",""))); instance=get_instance(table,record.get("query_genome_id",""),record.get("query_contig_id",""),record.get("query_contig_system_index","")) if table is not None else pd.DataFrame()
            instance=query_components(instance,record.get("query_protein_id",""))
            for _,q in instance.iterrows():
                seq=str(q.get("protein_seq","")).replace("*","")
                if seq:write(query_dir/f"{safe(q.get('component',''))}.faa",str(q.get("protein_id","")),seq)
        manifest.append(record)
    pd.DataFrame(manifest).to_csv(out/"instance_manifest.csv",index=False)
    print(f"QMA FAA complete: {len(manifest)} instance(s), {missing} missing mining sequence(s)")
if __name__=="__main__":main()
