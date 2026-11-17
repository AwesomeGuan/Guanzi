#!/usr/bin/env python3
"""Build an internal component-coordinate table from Prodigal GFF files."""
import argparse,os,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
KEY=["mining_genome_id","mining_contig_id","mining_component_protein_id"]
def args():
    p=argparse.ArgumentParser(); p.add_argument("--step4",required=True); p.add_argument("--prodigal_dir",required=True); p.add_argument("--output",required=True)
    p.add_argument("--check",action="store_true",help="Only verify that the existing table covers the current Step 4 components and exit 0/1 without writing")
    return p.parse_args()
def attrs(text):
    return {x.split("=",1)[0]:x.split("=",1)[1] for x in text.split(";") if "=" in x}
def worker_count():
    # GFF parsing is file-bound, so independent genomes are parsed concurrently
    # (the pre-rewrite pipeline used 8 threads); capped to avoid oversubscribing.
    return max(1,min(os.cpu_count() or 8,8))
def gff_rows(job):
    gff,genome,contig,ids=job; rows=[]
    for line in gff.read_text(encoding="utf-8",errors="ignore").splitlines():
        if not line or line.startswith("#"):continue
        fields=line.split("\t")
        if len(fields)<9 or fields[2]!="CDS":continue
        values=attrs(fields[8]); pid=values.get("ID") or values.get("locus_tag") or values.get("gene_id")
        if pid in ids:rows.append({"mining_genome_id":genome,"mining_contig_id":contig,"mining_component_protein_id":pid,"start":fields[3],"end":fields[4],"strand":fields[6]})
    return rows
def required_keys(step4):
    keys=set()
    for status in ("qualified","disqualified"):
        path=Path(step4)/f"comparison_{status}.csv"
        if not path.is_file():continue
        frame=pd.read_csv(path,dtype=str,keep_default_na=False)
        for _,r in frame.iterrows():keys.add(tuple(str(r[c]) for c in KEY))
    return keys
def existing_keys(path):
    try:frame=pd.read_csv(path,dtype=str,keep_default_na=False)
    except Exception:return None
    if any(c not in frame.columns for c in KEY):return None
    return {tuple(str(r[c]) for c in KEY) for _,r in frame.iterrows()}
def main():
    a=args(); keys=required_keys(a.step4); have=existing_keys(Path(a.output))
    if a.check:
        covered=0 if have is None else len(keys&have)
        print(f"coverage check: {covered} of {len(keys)} component key(s) already located")
        sys.exit(0 if have is not None and keys.issubset(have) else 1)
    wanted={}
    for genome,contig,pid in keys:wanted.setdefault((genome,contig),set()).add(pid)
    jobs=[]
    for gff in sorted(Path(a.prodigal_dir).rglob("*.gff")):
        genome=gff.parent.parent.name; contig=gff.parent.name; ids=wanted.get((genome,contig))
        if ids:jobs.append((gff,genome,contig,ids))
    rows=[]
    if jobs:
        with ThreadPoolExecutor(max_workers=worker_count()) as pool:
            for chunk in pool.map(gff_rows,jobs):rows.extend(chunk)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); pd.DataFrame(rows,columns=[*KEY,"start","end","strand"]).to_csv(out,index=False)
    located={(r["mining_genome_id"],r["mining_contig_id"],r["mining_component_protein_id"]) for r in rows}; missing=keys-located
    print(f"component locations complete: {len(rows)} row(s)")
    if missing:print(f"[S01] WARN unlocated_component={len(missing)} example={sorted(missing)[:3]}",flush=True)
if __name__=="__main__":main()
