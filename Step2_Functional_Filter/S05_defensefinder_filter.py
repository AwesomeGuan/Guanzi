#!/usr/bin/env python3
"""Filter DefenseFinder annotations and export annotated/novel candidates."""
import argparse
from pathlib import Path
import pandas as pd

from step2_candidate_utils import (
    enrich_candidates,
    finalize_candidate_table,
    join_unique,
    order_candidate_columns,
    remove_stale_csvs,
)


ANNOTATION_COLUMNS = ["gene_name_annotated", "subtype_annotated"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--defensefinder", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--final_output", required=True)
    return p.parse_args()


def read_tsv(path):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def add_annotation(annotation, hit_id, gene_name="", subtype=""):
    hit_id = str(hit_id).strip()
    if not hit_id:
        return
    entry = annotation.setdefault(hit_id, {"gene_name": [], "subtype": []})
    if str(gene_name).strip():
        entry["gene_name"].append(str(gene_name).strip())
    if str(subtype).strip():
        entry["subtype"].append(str(subtype).strip())


def load_defensefinder_annotations(root):
    annotation = {}
    root = Path(root)
    for pattern in ("*_defense_finder_hmmer.tsv", "*_defense_finder_genes.tsv"):
        for path in root.rglob(pattern):
            try:
                frame = read_tsv(path)
            except Exception:
                continue
            if "hit_id" not in frame.columns:
                continue
            for _, row in frame.iterrows():
                add_annotation(
                    annotation,
                    row.get("hit_id", ""),
                    row.get("gene_name", ""),
                    row.get("subtype", ""),
                )

    for path in root.rglob("*_defense_finder_systems.tsv"):
        try:
            frame = read_tsv(path)
        except Exception:
            continue
        if "protein_in_syst" not in frame.columns:
            continue
        for _, row in frame.iterrows():
            subtype = str(row.get("subtype", "")).strip()
            protein_ids = [
                item.strip()
                for item in str(row.get("protein_in_syst", "")).split(",")
                if item.strip()
            ]
            for hit_id in protein_ids:
                add_annotation(annotation, hit_id, subtype=subtype)
    return annotation


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    final_output = Path(args.final_output)
    final_output.mkdir(parents=True, exist_ok=True)
    annotation = load_defensefinder_annotations(args.defensefinder)
    annotated_ids = set(annotation)
    annotated_frames = []
    novel_frames = []
    for csv_path in sorted(Path(args.input).glob("*.csv")):
        frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        frame = enrich_candidates(frame, csv_path.stem, {})
        is_annotated = frame["mining_protein_id"].astype(str).isin(annotated_ids)
        keep = order_candidate_columns(frame.loc[~is_annotated])
        annotated = frame.loc[is_annotated].copy()
        annotated["gene_name_annotated"] = annotated["mining_protein_id"].astype(str).map(
            lambda hit_id: join_unique(annotation.get(hit_id, {}).get("gene_name", []))
        )
        annotated["subtype_annotated"] = annotated["mining_protein_id"].astype(str).map(
            lambda hit_id: join_unique(annotation.get(hit_id, {}).get("subtype", []))
        )
        annotated_frames.append(order_candidate_columns(annotated, ANNOTATION_COLUMNS))
        novel_frames.append(keep)
        keep.to_csv(out_dir / csv_path.name, index=False)

    final_annotated = finalize_candidate_table(annotated_frames, ANNOTATION_COLUMNS)
    final_novel = finalize_candidate_table(novel_frames)
    final_annotated.to_csv(
        final_output / "isozyme_candidates_defensefinder_annotated.csv", index=False
    )
    final_novel.to_csv(final_output / "isozyme_candidates_novel.csv", index=False)
    remove_stale_csvs(out_dir, {path.name for path in Path(args.input).glob("*.csv")})
    print("DefenseFinder novelty filter complete")


if __name__ == "__main__":
    main()
