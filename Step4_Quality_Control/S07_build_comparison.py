#!/usr/bin/env python3
"""Build the canonical component Comparison table from TM-align outputs."""
import argparse,re
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from step4_utils import *

COMPARISON_BASE=["defense_system","query_genome_id","query_contig_id","query_contig_system_index","query_component_protein_id","query_component"]
COMPARISON_MINING=["mining_genome_id","mining_contig_id","mining_contig_system_index","mining_component_protein_id","mining_component","is_candidate_mining"]
METRICS=["tm_score","sequence_identity","aligned_length","rmsd"]
METRIC_SUFFIX=".tmalign.txt"
def args():
    p=argparse.ArgumentParser(); p.add_argument("--qma",required=True); p.add_argument("--assembly",required=True); p.add_argument("--tmalign_results",required=True); p.add_argument("--output",required=True); return p.parse_args()
def parse(path):
    text=path.read_text(encoding="utf-8",errors="ignore")
    tm=re.search(r"TM-score=\s*([\d.]+)",text)
    identity=None
    for pattern in (
        r"Seq_ID\s*=\s*n_identical/n_aligned\s*=\s*([\d.]+)",
        r"Seq_ID\s*=\s*([\d.]+)",
        r"(?<![A-Za-z0-9_])ID\s*=\s*([\d.]+)",
    ):
        identity=re.search(pattern,text)
        if identity:break
    align=re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+)",text)
    return (tm.group(1),identity.group(1),align.group(1),align.group(2)) if tm and identity and align else None
def fid(path):
    if not path.is_file():return ""
    try:return next(SeqIO.parse(path,"fasta")).id
    except StopIteration:return ""
def main():
    a=args(); qma=Path(a.qma); manifest=read_csv(qma/"instance_manifest.csv"); assembly=read_csv(a.assembly); rows=[]
    rules=load_rules(PROJECT_ROOT/"data"/"DefenseFinder_rules"); pairings={}
    counters={"metric_files":0,"parse_failed":0,"manifest_missing":0,"assembly_missing":0,"faa_missing":0}
    metric_files=sorted(Path(a.tmalign_results).glob(f"*/*/*/*{METRIC_SUFFIX}"))
    for result in metric_files:
        counters["metric_files"]+=1
        metrics=parse(result)
        if metrics is None:
            counters["parse_failed"]+=1;continue
        system,context,index,_=result.relative_to(a.tmalign_results).parts; match=manifest[(manifest["defense_system"]==system)&(manifest["instance_path_context"]==context)&(manifest["mining_contig_system_index"]==index)]
        if match.empty:
            counters["manifest_missing"]+=1;continue
        meta=match.iloc[0]; key=instance_key(meta); ar=assembly[key_mask(assembly,key)]
        if ar.empty:
            counters["assembly_missing"]+=1;continue
        ar=ar.iloc[0]; component=result.name[:-len(METRIC_SUFFIX)]; leaf=qma/system/context/index
        pairing_key=(system,context,index)
        if pairing_key not in pairings:
            pairings[pairing_key]=pair_component_names(system,
                [p.stem for p in (leaf/"mining_assembly").glob("*.faa")],
                [p.stem for p in (leaf/"query_assembly").glob("*.faa")],rules)
        query_component=pairings[pairing_key].get(component)
        if query_component is None:
            counters["faa_missing"]+=1;continue
        mining_faa=leaf/"mining_assembly"/f"{component}.faa"; query_faa=leaf/"query_assembly"/f"{query_component}.faa"; mpid=fid(mining_faa); qpid=fid(query_faa)
        if not mpid or not qpid:
            counters["faa_missing"]+=1;continue
        row={"defense_system":system,"query_genome_id":meta.get("query_genome_id",""),"query_contig_id":meta.get("query_contig_id",""),"query_contig_system_index":meta.get("query_contig_system_index",""),"query_component_protein_id":qpid,"query_component":query_component}
        if "mining_domain" in ar.index:row["mining_domain"]=ar["mining_domain"]
        row.update({"mining_genome_id":ar["mining_genome_id"],"mining_contig_id":ar["mining_contig_id"],"mining_contig_system_index":ar["mining_contig_system_index"],"mining_component_protein_id":mpid,"mining_component":component,"is_candidate_mining":"T" if mpid==str(ar.get("mining_component_protein_id_1","")) else "F","tm_score":metrics[0],"sequence_identity":metrics[1],"aligned_length":metrics[2],"rmsd":metrics[3]}); rows.append(row)
    frame=pd.DataFrame(rows); columns=COMPARISON_BASE+(["mining_domain"] if "mining_domain" in frame.columns else [])+COMPARISON_MINING+METRICS; frame=order_frame(frame,columns).drop_duplicates(); frame.to_csv(a.output,index=False)
    candidate_rows=int(frame["is_candidate_mining"].eq("T").sum()) if not frame.empty else 0
    print("component Comparison audit: "+", ".join(f"{key}={value}" for key,value in counters.items())+f", rows={len(frame)}, candidate_rows={candidate_rows}")
    if not assembly.empty and not metric_files:
        raise SystemExit("No canonical TM-align metric files were found for non-empty Assembly input")
    if not assembly.empty and frame.empty:
        raise SystemExit("TM-align metrics produced zero Comparison rows; see audit counters")
    if not assembly.empty and candidate_rows==0:
        raise SystemExit("Comparison contains no candidate Mining component rows")
    print(f"component Comparison complete: {len(frame)} row(s)")
if __name__=="__main__":main()
