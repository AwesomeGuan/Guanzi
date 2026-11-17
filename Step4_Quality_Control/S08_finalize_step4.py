#!/usr/bin/env python3
"""Split all three formal file classes by candidate-component TM-align QC."""
import argparse
from pathlib import Path
import pandas as pd
from step4_utils import *
def args():
    p=argparse.ArgumentParser(); p.add_argument("--mapping",required=True); p.add_argument("--assembly",required=True); p.add_argument("--comparison",required=True); p.add_argument("--output",required=True); p.add_argument("--tm_score",type=float,default=.7); p.add_argument("--max_identity",type=float,default=.3); return p.parse_args()
def main():
    a=args(); mapping=read_csv(a.mapping); assembly=read_csv(a.assembly); comparison=read_csv(a.comparison)
    if not assembly.empty and comparison.empty:
        raise SystemExit("Comparison is empty for non-empty Assembly input")
    candidate=comparison[comparison["is_candidate_mining"].eq("T")].copy()
    if not assembly.empty and candidate.empty:
        raise SystemExit("Comparison contains no candidate Mining component rows")
    tm=pd.to_numeric(candidate["tm_score"],errors="coerce"); identity=pd.to_numeric(candidate["sequence_identity"],errors="coerce"); qualified=table_keys(candidate[(tm>=a.tm_score)&(identity<=a.max_identity)]); evaluated=table_keys(candidate); groups={"qualified":qualified,"disqualified":evaluated-qualified}; out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    for status,keys in groups.items():
        subset_keys(mapping,keys).to_csv(out/f"mapping_{status}.csv",index=False); subset_keys(assembly,keys).to_csv(out/f"assembly_{status}.csv",index=False); subset_keys(comparison,keys).to_csv(out/f"comparison_{status}.csv",index=False)
    print(f"Step4 formal outputs: {len(qualified)} qualified, {len(evaluated-qualified)} disqualified, {len(table_keys(assembly)-evaluated)} structurally unassessed instance(s)")
if __name__=="__main__":main()
