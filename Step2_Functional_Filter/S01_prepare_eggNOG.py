#!/usr/bin/env python3
"""Prepare mining FAA for eggNOG upload with one-pass contig extraction."""
import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

import sys
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from utils.pipeline_schema import canonicalize


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--protein_faa", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def extract_wanted_sequences(faa_path, wanted_ids):
    """Read one contig FAA once and return only requested protein sequences."""
    found = {}
    current_id = None
    sequence = []

    def store_current():
        if current_id in wanted_ids and current_id not in found:
            found[current_id] = "".join(sequence).replace("*", "")

    with Path(faa_path).open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store_current()
                if len(found) == len(wanted_ids):
                    break
                current_id = line[1:].split()[0]
                sequence = []
            elif current_id is not None:
                sequence.append(line)
        else:
            store_current()
    return found


def load_targets(mapping_root):
    """Group unique mining proteins by their explicit genome/contig location."""
    targets = defaultdict(set)
    failures = []
    csv_paths = sorted(Path(mapping_root).glob("*.csv"))
    if not csv_paths:
        raise FileNotFoundError(f"No mapping CSV files found under {mapping_root}")

    required = {"mining_protein_id"}
    unlocated = set()
    for csv_path in csv_paths:
        frame = canonicalize(pd.read_csv(csv_path, dtype=str, keep_default_na=False))
        missing_columns = required - set(frame.columns)
        if missing_columns:
            raise ValueError(
                f"{csv_path} is missing columns: {sorted(missing_columns)}"
            )
        for _, row in frame.iterrows():
            genome_id = str(row.get("mining_genome_id", "")).strip()
            contig_id = str(row.get("mining_contig_id", "")).strip()
            mining_id = str(row["mining_protein_id"]).strip()
            if not mining_id:
                failures.append((
                    genome_id,
                    contig_id,
                    mining_id,
                    "mapping row lacks mining_protein_id",
                ))
                continue
            if genome_id and contig_id:
                targets[(genome_id, contig_id)].add(mining_id)
            else:
                unlocated.add(mining_id)
    return targets, unlocated, failures


def resolve_unlocated(prodigal_root, wanted):
    """Resolve Custom-Mode IDs in one pass over FAA headers."""
    resolved = defaultdict(set)
    remaining = set(wanted)
    for faa_path in Path(prodigal_root).rglob("*.faa"):
        if not remaining:
            break
        present = set()
        with faa_path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    protein_id = line[1:].split()[0]
                    if protein_id in remaining:
                        present.add(protein_id)
        if present:
            resolved[(faa_path.parent.parent.name, faa_path.parent.name)].update(present)
            remaining.difference_update(present)
    return resolved, remaining


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_faa = out_dir / "eggnog_upload.faa"
    temporary_faa = out_dir / ".eggnog_upload.faa.tmp"
    failure_path = out_dir / "eggnog_extraction_failures.csv"
    prodigal_root = Path(args.protein_faa)

    targets, unlocated, failures = load_targets(args.input)
    resolved, unresolved = resolve_unlocated(prodigal_root, unlocated)
    for key, protein_ids in resolved.items():
        targets[key].update(protein_ids)
    failures.extend(("", "", protein_id, "protein location not found in Prodigal FAA") for protein_id in sorted(unresolved))
    seen = set()
    total_contigs = len(targets)
    with temporary_faa.open("w", encoding="utf-8") as handle:
        for index, ((genome_id, contig_id), mining_ids) in enumerate(
            sorted(targets.items()), start=1
        ):
            faa_path = (
                prodigal_root / genome_id / contig_id / f"{contig_id}.faa"
            )
            if not faa_path.is_file():
                failures.extend(
                    (genome_id, contig_id, mining_id, f"FAA file not found: {faa_path}")
                    for mining_id in sorted(mining_ids)
                )
                continue
            try:
                found = extract_wanted_sequences(faa_path, mining_ids)
            except Exception as exc:
                failures.extend(
                    (
                        genome_id,
                        contig_id,
                        mining_id,
                        f"unexpected {type(exc).__name__}: {exc}",
                    )
                    for mining_id in sorted(mining_ids)
                )
                continue

            for mining_id in sorted(mining_ids):
                sequence = found.get(mining_id)
                if sequence is None:
                    failures.append((
                        genome_id,
                        contig_id,
                        mining_id,
                        f"protein ID absent from {faa_path}",
                    ))
                elif mining_id not in seen:
                    handle.write(f">{mining_id}\n{sequence}\n")
                    seen.add(mining_id)

            if index % 500 == 0 or index == total_contigs:
                print(
                    f"eggNOG extraction: {index}/{total_contigs} contigs, "
                    f"{len(seen)} unique proteins"
                )

    temporary_faa.replace(out_faa)
    if failures:
        pd.DataFrame(
            failures,
            columns=["mining_genome_id", "mining_contig_id", "mining_protein_id", "reason"],
        ).to_csv(failure_path, index=False)
        print(
            f"WARNING: skipped {len(failures)} extraction failure(s); "
            f"see {failure_path}"
        )
    else:
        failure_path.unlink(missing_ok=True)
    print(f"eggNOG upload: {out_faa} ({len(seen)} sequences)")


if __name__ == "__main__":
    main()
