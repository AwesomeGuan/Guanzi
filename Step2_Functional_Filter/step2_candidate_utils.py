#!/usr/bin/env python3
"""Shared helpers for canonical Step2 protein-candidate tables."""

import sys
import argparse
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from utils.pipeline_schema import canonicalize, order_frame, step2_columns


def extract_contig(protein_id):
    parts = str(protein_id).rsplit("_", 1)
    return parts[0] if len(parts) == 2 and parts[1].isdigit() else str(protein_id)


def build_contig_genome_index(prodigal_root):
    index = {}
    for faa in Path(prodigal_root).rglob("*.faa"):
        index.setdefault(faa.stem, faa.parent.parent.name)
    return index


def enrich_candidates(frame, defense_system, contig_genome_index):
    frame = canonicalize(frame)
    frame["defense_system"] = str(defense_system)
    if "mining_protein_id" not in frame.columns:
        raise ValueError("Candidate table is missing mining_protein_id")
    inferred_contigs = frame["mining_protein_id"].astype(str).map(extract_contig)
    if "mining_contig_id" not in frame.columns:
        frame["mining_contig_id"] = inferred_contigs
    else:
        empty = frame["mining_contig_id"].astype(str).str.strip().eq("")
        frame.loc[empty, "mining_contig_id"] = inferred_contigs[empty]
    inferred_genomes = frame["mining_contig_id"].astype(str).map(contig_genome_index).fillna("Unknown")
    if "mining_genome_id" not in frame.columns:
        frame["mining_genome_id"] = inferred_genomes
    else:
        empty = frame["mining_genome_id"].astype(str).str.strip().isin(["", "Unknown"])
        frame.loc[empty, "mining_genome_id"] = inferred_genomes[empty]
    return frame


def order_candidate_columns(frame, extra_columns=None):
    frame = canonicalize(frame)
    columns = step2_columns(frame)
    for column in list(extra_columns or []):
        if column not in columns:
            columns.append(column)
    return order_frame(frame, columns)


def finalize_candidate_table(frames, extra_columns=None):
    if not frames:
        return pd.DataFrame()
    merged = pd.concat([canonicalize(frame) for frame in frames], ignore_index=True)
    merged = order_candidate_columns(merged, extra_columns)
    key = [column for column in [
        "defense_system", "query_protein_id", "mining_protein_id",
        "mining_genome_id", "mining_contig_id",
    ] if column in merged.columns]
    return merged.drop_duplicates(subset=key, keep="first")


def join_unique(values):
    unique = []
    for value in values:
        text = str(value).strip()
        if text and text.lower() != "nan" and text not in unique:
            unique.append(text)
    return ";".join(unique)


def remove_stale_csvs(directory, expected):
    """Remove only obsolete CSV outputs owned by this operation."""
    for path in Path(directory).glob("*.csv"):
        if path.name not in expected:
            path.unlink()


def clear_classification(root):
    root = Path(root)
    for group in ("single", "multi", "complex", "cas"):
        remove_stale_csvs(root / "step2_output" / f"step2_mapping_{group}", set())
        (root / "step2_output" / f"step2_summary_{group}.csv").unlink(missing_ok=True)


def cache_valid(marker, inputs, outputs, mirrors):
    """Check concrete outputs and dependency modification times, not just .done."""
    marker = Path(marker)
    if not marker.is_file():
        return False
    for value in inputs:
        path = Path(value)
        if not path.exists():
            return False
        # Directories catch removed files; their children catch edited files.
        paths = [path]
        if path.is_dir():
            paths.extend(path.rglob("*"))
        if any(item.stat().st_mtime_ns > marker.stat().st_mtime_ns for item in paths):
            return False
    for value in outputs:
        path = Path(value)
        if not path.exists():
            return False
        csv_paths = [path] if path.is_file() and path.suffix == ".csv" else list(path.glob("*.csv"))
        for csv_path in csv_paths:
            frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
            if "mining_protein_id" not in frame.columns:
                return False
    for source, target in mirrors:
        if {p.name for p in Path(source).glob("*.csv")} != {p.name for p in Path(target).glob("*.csv")}:
            return False
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step2 cache checks and scoped generated-table cleanup")
    parser.add_argument("action", choices=("check", "clear-classification", "clear-novelty", "clear-eggnog"))
    parser.add_argument("--marker")
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--output", action="append", default=[])
    parser.add_argument("--mirror", nargs=2, action="append", default=[])
    parser.add_argument("--root")
    args = parser.parse_args()
    if args.action == "check":
        if not args.marker:
            parser.error("--marker is required for check")
        try:
            valid = cache_valid(args.marker, args.input, args.output, args.mirror)
        except (OSError, ValueError, pd.errors.ParserError):
            valid = False
        sys.exit(0 if valid else 1)
    if not args.root:
        parser.error("--root is required for cleanup")
    clear_classification(args.root)
    if args.action in ("clear-novelty", "clear-eggnog"):
        remove_stale_csvs(Path(args.root) / "step2_raw/defensefinder_filtered", set())
        for name in ("isozyme_candidates_novel.csv", "isozyme_candidates_defensefinder_annotated.csv"):
            (Path(args.root) / name).unlink(missing_ok=True)
    if args.action == "clear-eggnog":
        remove_stale_csvs(Path(args.root) / "step2_raw/eggnog_filtered", set())
        (Path(args.root) / "isozyme_candidates_eggnog_annotated.csv").unlink(missing_ok=True)
