#!/usr/bin/env python3
"""Run DefenseFinder on existing Prodigal protein predictions."""
import argparse
import subprocess
from pathlib import Path
import pandas as pd
import sys
import json
import tempfile
from eggnog_annotation_utils import file_digest
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path: sys.path.insert(0, str(PROJECT_ROOT))
from utils.pipeline_schema import canonicalize


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--models_dir")
    return p.parse_args()


def extract_contig(gene_id):
    parts = str(gene_id).rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return str(gene_id)


def candidate_contigs(reference_root):
    candidates = set()
    for csv_path in Path(reference_root).glob("*.csv"):
        frame = canonicalize(pd.read_csv(csv_path, dtype=str, keep_default_na=False))
        if "mining_protein_id" not in frame.columns:
            continue
        for _, row in frame.iterrows():
            genome_id = str(row.get("mining_genome_id", "")).strip()
            contig_id = str(row.get("mining_contig_id", "")).strip()
            if not contig_id:
                contig_id = extract_contig(row.get("mining_protein_id", ""))
            if not genome_id or genome_id == "Unknown":
                raise ValueError(
                    f"{csv_path} lacks a valid genome_id for contig {contig_id}"
                )
            if contig_id:
                candidates.add((genome_id, contig_id))
    return candidates


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    prodigal_root = Path(args.input)
    contigs = candidate_contigs(args.reference)
    reused = 0
    new = 0
    for genome_id, contig_id in sorted(contigs):
        faa = prodigal_root / genome_id / contig_id / f"{contig_id}.faa"
        if not faa.is_file():
            raise FileNotFoundError(f"Prodigal FAA not found: {faa}")
        target = out_dir / genome_id / contig_id
        target.mkdir(parents=True, exist_ok=True)
        record = target / ".step2_defensefinder.json"
        context = {"faa_sha256": file_digest(faa), "models_dir": args.models_dir or "bundled"}
        try:
            previous = json.loads(record.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            previous = None
        tables = {}
        for kind, column in (("hmmer", "hit_id"), ("genes", "hit_id"), ("systems", "protein_in_syst")):
            matches = list(target.glob(f"*_defense_finder_{kind}.tsv"))
            if len(matches) != 1:
                break
            try:
                frame = pd.read_csv(matches[0], sep="\t", dtype=str, keep_default_na=False)
            except (ValueError, pd.errors.ParserError):
                break
            if column not in frame.columns:
                break
            tables[matches[0].name] = file_digest(matches[0])
        if previous and previous.get("context") == context and len(tables) == 3 and previous.get("tables") == tables:
            reused += 1
            continue
        # Stage a fresh tool run, avoiding conflicts with partial or legacy tool directories.
        with tempfile.TemporaryDirectory(prefix=".step2_df_", dir=target.parent) as temporary:
            work = Path(temporary)
            cmd = ["defense-finder", "run", str(faa), "-o", str(work), "--db-type", "ordered_replicon"]
            if args.models_dir:
                cmd += ["--models-dir", args.models_dir]
            subprocess.run(cmd, check=True)
            tables = {}
            for kind, column in (("hmmer", "hit_id"), ("genes", "hit_id"), ("systems", "protein_in_syst")):
                matches = list(work.glob(f"*_defense_finder_{kind}.tsv"))
                if len(matches) != 1 or column not in pd.read_csv(matches[0], sep="\t", nrows=0).columns:
                    raise ValueError(f"DefenseFinder did not create a valid {kind} table: {work}")
                tables[matches[0].name] = file_digest(matches[0])
            for obsolete in target.glob("*_defense_finder_*.tsv"):
                if obsolete.name not in tables:
                    obsolete.unlink()
            for name in tables:
                (work / name).replace(target / name)
        record.write_text(json.dumps({"context": context, "tables": tables}, indent=2) + "\n", encoding="utf-8")
        new += 1
    print(f"DefenseFinder done: total_contigs={len(contigs)} reused={reused} new={new}")


if __name__ == "__main__":
    main()
