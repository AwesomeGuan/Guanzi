#!/usr/bin/env python3
"""Copy Mapping/Assembly and add transcription order plus syntax to Comparison.

This operation owns all six final tables: they are cleared before the new ones
are written, so an interrupted run cannot leave a complete-looking stale set.

An instance is kept only when all located components share a coding strand. It is
dropped whole, from all six tables, when either

  * any component has no Prodigal coordinate row, or
  * its located components span more than one strand (mixed strand).

The second rule keeps a single `[Chain]` value and coding-direction order
well defined. Same-strand coordinates do not establish co-transcription.
Mapping, Assembly and Comparison always describe the same retained instance set.

The syntax follows the nomenclature rules of the paper. Subunits are lettered by
their transcriptional order (forward-strand subunits ascending by start, then
reverse-strand subunits descending), the mining subunit that guided discovery is
the start point, and subunits meeting the structural standard (CDP-normalized
TM-score >= 0.7 and
sequence identity <= 0.3) are marked with an asterisk:

    (Pre-group)*-Start-(Post-group)*[Chain]

Empty groups are omitted; a group carries one asterisk per conforming subunit.
The start point has its own asterisk when conforming. Group stars count
conforming members; individual metrics remain in Comparison. The coding-strand
direction is taken from the start point.
"""
import argparse
from pathlib import Path
import pandas as pd

LETTERS="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
KEY=["defense_system","mining_genome_id","mining_contig_id","mining_contig_system_index"]
LOCATION=["mining_genome_id","mining_contig_id","mining_component_protein_id"]
KINDS=("mapping","assembly","comparison")
TM_SCORE=0.7
MAX_IDENTITY=0.3


def args():
    p=argparse.ArgumentParser()
    p.add_argument("--step4",required=True); p.add_argument("--locations",required=True); p.add_argument("--output",required=True)
    return p.parse_args()


def label(position):
    return LETTERS[position] if position<len(LETTERS) else f"{LETTERS[position//26-1]}{LETTERS[position%26]}"


def instance_keys(frame):
    return pd.Series([tuple(row) for row in frame[KEY].to_numpy()],index=frame.index)


def order_index(group):
    forward=group[group["_strand"].ne("-")].sort_values(["_start","mining_component_protein_id"],kind="stable")
    reverse=group[group["_strand"].eq("-")].sort_values(["_start","mining_component_protein_id"],kind="stable",ascending=[False,True])
    return [*forward.index,*reverse.index]


def group_syntax(group,order):
    if not order:return ""
    letters={index:label(position) for position,index in enumerate(order)}
    conforming={index:bool(group.at[index,"_conforming"]) for index in order}
    start=next((index for index in order if group.at[index,"is_candidate_mining"]=="T"),order[0])
    cut=order.index(start); front=order[:cut]; back=order[cut+1:]
    parts=[]
    if front:parts.append("("+"".join(letters[index] for index in front)+")"+"*"*sum(conforming[index] for index in front))
    parts.append(letters[start]+("*" if conforming[start] else ""))
    if back:parts.append("("+"".join(letters[index] for index in back)+")"+"*"*sum(conforming[index] for index in back))
    return "-".join(parts)+f"[{group.at[start,'_strand']}]"


def main():
    a=args(); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    for status in ("qualified","disqualified"):
        for kind in KINDS:(out/f"{kind}_{status}.csv").unlink(missing_ok=True)
    locations=pd.read_csv(a.locations,dtype=str,keep_default_na=False).drop_duplicates(subset=LOCATION)
    for status in ("qualified","disqualified"):
        comparison=pd.read_csv(Path(a.step4)/f"comparison_{status}.csv",dtype=str,keep_default_na=False)
        merged=comparison.merge(locations,on=LOCATION,how="left")
        merged["_start"]=pd.to_numeric(merged["start"],errors="coerce")
        merged["_strand"]=merged["strand"].astype(str).where(merged["_start"].notna(),"")
        merged["_key"]=instance_keys(merged)
        missing={key for key,group in merged.groupby("_key",sort=False) if group["_start"].isna().any()}
        mixed={key for key,group in merged.groupby("_key",sort=False) if group.loc[group["_start"].notna(),"_strand"].nunique()>1}
        dropped=missing|mixed
        kept=merged[~merged["_key"].isin(dropped)].copy()
        kept["_conforming"]=(pd.to_numeric(kept["tm_score"],errors="coerce")>=TM_SCORE)&(pd.to_numeric(kept["sequence_identity"],errors="coerce")<=MAX_IDENTITY)
        kept["transcription_order"]=""; kept["syntax_representation"]=""
        instances=0
        for _,group in kept.groupby("_key",sort=False):
            order=order_index(group); syntax=group_syntax(group,order); instances+=1
            for position,index in enumerate(order,1):kept.at[index,"transcription_order"]=str(position)
            for index in group.index:kept.at[index,"syntax_representation"]=syntax
        drop=[c for c in ["start","end","strand","_start","_strand","_conforming","_key"] if c in kept.columns]
        kept.drop(columns=drop).to_csv(out/f"comparison_{status}.csv",index=False)
        for kind in ("mapping","assembly"):
            frame=pd.read_csv(Path(a.step4)/f"{kind}_{status}.csv",dtype=str,keep_default_na=False)
            absent=[column for column in KEY if column not in frame.columns]
            if absent:raise ValueError(f"{kind}_{status}.csv is missing instance key columns: {absent}")
            frame["_key"]=instance_keys(frame)
            frame[~frame["_key"].isin(dropped)].drop(columns="_key").to_csv(out/f"{kind}_{status}.csv",index=False)
        print(f"comparison_{status}: {len(kept)} row(s), {instances} instance(s), {len(dropped)} instance(s) dropped "
              f"(missing_coordinates={len(missing)}, mixed_strand={len(mixed)})")
        if dropped:print(f"[S02] WARN dropped_instances={len(dropped)} missing_coordinates={len(missing)} mixed_strand={len(mixed)} example={sorted(dropped)[:3]}",flush=True)
    print("Step5 six-file output complete")


if __name__=="__main__":main()
