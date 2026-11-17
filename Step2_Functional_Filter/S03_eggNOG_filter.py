#!/usr/bin/env python3
"""Filter eggNOG defense annotations and export annotated candidates."""
import argparse
from pathlib import Path
import pandas as pd
from eggnog_annotation_utils import read_annotation_rows

from step2_candidate_utils import (
    build_contig_genome_index,
    enrich_candidates,
    finalize_candidate_table,
    join_unique,
    order_candidate_columns,
    remove_stale_csvs,
)


ANNOTATION_COLUMNS = [
    "Description_annotated",
    "PFAMs_annotated",
    "COG_category_annotated",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--eggNOG_tsv", dest="egg_nog_tsv", required=True)
    p.add_argument("--description_keywords", required=True)
    p.add_argument("--pfam_keywords", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--prodigal", required=True)
    p.add_argument("--final_output", required=True)
    return p.parse_args()


def read_eggnog_table(path):
    header, rows = read_annotation_rows(path)
    return pd.DataFrame(rows, columns=header, dtype=str)


def load_keyword_rules(path):
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "defense_system" not in frame.columns:
        raise ValueError(f"Keyword table is missing defense_system: {path}")
    rules = {}
    keyword_columns = [column for column in frame.columns if column != "defense_system"]
    for _, row in frame.iterrows():
        system = str(row.get("defense_system", "")).strip()
        if not system:
            continue
        keywords = []
        for column in keyword_columns:
            keyword = str(row.get(column, "")).strip().strip('"').strip()
            if keyword and keyword not in keywords:
                keywords.append(keyword)
        rules[system] = keywords
    return rules


def contains_any(value, keywords):
    text = str(value).casefold()
    return any(keyword.casefold() in text for keyword in keywords)


def build_annotation_table(egg, query_col):
    rename = {
        "Description": "Description_annotated",
        "PFAMs": "PFAMs_annotated",
        "COG_category": "COG_category_annotated",
    }
    missing = set(rename) - set(egg.columns)
    if missing:
        raise ValueError(f"eggNOG table lacks required annotation columns: {sorted(missing)}")
    annotated = egg.loc[:, [query_col, *rename]].copy()
    annotated = annotated.rename(columns={query_col: "mining_protein_id", **rename})
    if annotated.empty:
        return annotated[["mining_protein_id", *ANNOTATION_COLUMNS]]
    return annotated.groupby("mining_protein_id", as_index=False)[ANNOTATION_COLUMNS].agg(join_unique)


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    final_output = Path(args.final_output)
    final_output.mkdir(parents=True, exist_ok=True)
    contig_genome_index = build_contig_genome_index(args.prodigal)
    egg = read_eggnog_table(args.egg_nog_tsv)
    query_col = "#query" if "#query" in egg.columns else "query"
    if query_col not in egg.columns:
        raise ValueError("eggNOG table is missing #query/query column")
    annotation = build_annotation_table(egg, query_col)
    description_rules = load_keyword_rules(args.description_keywords)
    pfam_rules = load_keyword_rules(args.pfam_keywords)
    missing_description_rules = set()
    missing_pfam_rules = set()
    annotated_frames = []
    for csv_path in sorted(Path(args.input).glob("*.csv")):
        frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        system = csv_path.stem
        frame = enrich_candidates(frame, system, contig_genome_index)
        joined = frame[["mining_protein_id"]].merge(annotation, on="mining_protein_id", how="left")
        for column in ANNOTATION_COLUMNS:
            joined[column] = joined[column].fillna("")

        description_keywords = description_rules.get(system, [])
        pfam_keywords = pfam_rules.get(system, [])
        if system not in description_rules:
            missing_description_rules.add(system)
        if system not in pfam_rules:
            missing_pfam_rules.add(system)

        description_hit = joined["Description_annotated"].map(
            lambda value: contains_any(value, description_keywords)
        )
        pfam_hit = joined["PFAMs_annotated"].map(
            lambda value: contains_any(value, pfam_keywords)
        )
        cog_hit = joined["COG_category_annotated"].astype(str).str.contains(
            "V", case=False, regex=False, na=False
        )
        is_annotated = description_hit | pfam_hit | cog_hit
        keep = order_candidate_columns(frame.loc[~is_annotated])
        annotated = frame.loc[is_annotated].copy()
        for column in ANNOTATION_COLUMNS:
            annotated[column] = joined.loc[is_annotated, column].to_numpy()
        annotated_frames.append(order_candidate_columns(annotated, ANNOTATION_COLUMNS))
        keep.to_csv(out_dir / csv_path.name, index=False)

    final = finalize_candidate_table(annotated_frames, ANNOTATION_COLUMNS)
    final.to_csv(final_output / "isozyme_candidates_eggnog_annotated.csv", index=False)
    remove_stale_csvs(out_dir, {path.name for path in Path(args.input).glob("*.csv")})
    if missing_description_rules:
        print("WARNING: no Description keywords for: " + ", ".join(sorted(missing_description_rules)))
    if missing_pfam_rules:
        print("WARNING: no PFAM keywords for: " + ", ".join(sorted(missing_pfam_rules)))
    print("eggNOG filter complete")


if __name__ == "__main__":
    main()
