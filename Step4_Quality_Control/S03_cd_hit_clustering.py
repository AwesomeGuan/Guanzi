#!/usr/bin/env python3
"""Run CD-HIT clustering per defense-system subunit type."""
import argparse
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cdhit_output", required=True)
    parser.add_argument("--qma_faa", required=True)
    parser.add_argument("--threads", type=int, default=12)
    add_progress_argument(parser)
    return parser.parse_args()


def sequence_id_from_record(record):
    return record.id.strip()


def parse_cluster_file(path):
    seq_to_cluster = {}
    current_cluster = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">Cluster"):
            current_cluster = line.split()[1]
        elif ">" in line:
            seq_id = line.split(">")[-1].split("...")[0].strip()
            if seq_id:
                seq_to_cluster[seq_id] = current_cluster or "0"
    return seq_to_cluster


def main():
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    output_dir = Path(args.cdhit_output)
    output_dir.mkdir(parents=True, exist_ok=True)

    progress = Progress(
        args.log_label or "CD-HIT",
        phase="collect-sequences",
        processed=0,
        clustered=0,
        fallback=0,
        threads=args.threads,
    )
    processed = clustered = fallback = 0
    qma_root = Path(args.qma_faa)
    mapping_rows = []
    subunit_sequences = defaultdict(lambda: defaultdict(list))

    for system_dir in qma_root.iterdir():
        if not system_dir.is_dir():
            continue
        for contig_dir in system_dir.iterdir():
            if not contig_dir.is_dir():
                continue
            for leaf_dir in contig_dir.iterdir():
                if not leaf_dir.is_dir():
                    continue
                mining_assembly = leaf_dir / "mining_assembly"
                if not mining_assembly.exists():
                    continue
                for faa_path in mining_assembly.glob("*.faa"):
                    try:
                        records = list(SeqIO.parse(faa_path, "fasta"))
                    except Exception as exc:
                        progress.message(f"{faa_path}: {type(exc).__name__}: {exc}", "WARN")
                        continue
                    for record in records:
                        seq_id = sequence_id_from_record(record)
                        subunit_name = faa_path.stem
                        sequence = str(record.seq)
                        mapping_rows.append(
                            {
                                "defense_system": system_dir.name,
                                "contig_id": contig_dir.name,
                                "mining_id": leaf_dir.name,
                                "subunit_name": subunit_name,
                                "sequence_id": seq_id,
                                "sequence_length": len(sequence),
                            }
                        )
                        subunit_sequences[system_dir.name][subunit_name].append(
                            (seq_id, sequence)
                        )

    all_clusters = defaultdict(lambda: defaultdict(dict))
    cdhit_available = shutil.which("cd-hit") is not None
    progress.set(total=sum(len(subunits) for subunits in subunit_sequences.values()))

    for system, subunits in subunit_sequences.items():
        system_temp = output_dir / system
        system_temp.mkdir(parents=True, exist_ok=True)
        for subunit_name, sequences in subunits.items():
            progress.set(phase="clustering", system=system, component=subunit_name, sequences=len(sequences))
            progress.message(f"system={system} component={subunit_name} sequences={len(sequences)}")
            merged_faa = system_temp / f"{subunit_name}_input.faa"
            with merged_faa.open("w", encoding="utf-8") as handle:
                for seq_id, sequence in sequences:
                    handle.write(f">{seq_id}\n{sequence}\n")

            if len(sequences) < 2 or not cdhit_available:
                for seq_id, _ in sequences:
                    all_clusters[system][subunit_name][seq_id] = "0"
                processed += 1
                fallback += 1
                reason = "single sequence" if len(sequences) < 2 else "cd-hit unavailable"
                progress.message(f"{reason}; applying existing single-cluster handling",
                                 "INFO" if len(sequences) < 2 else "WARN")
                progress.set(processed=processed, clustered=clustered, fallback=fallback)
                continue

            output_prefix = system_temp / f"{subunit_name}_clusters"
            command = [
                "cd-hit",
                "-i",
                str(merged_faa),
                "-o",
                str(output_prefix),
                "-c",
                "0.5",
                "-n",
                "3",
                "-d",
                "0",
                "-M",
                "16000",
                "-T",
                str(args.threads),
            ]
            result = subprocess.run(command, check=False, capture_output=True, text=True)
            cluster_file = output_prefix.with_suffix(".clstr")
            if result.returncode == 0 and cluster_file.exists():
                parsed = parse_cluster_file(cluster_file)
                if parsed:
                    all_clusters[system][subunit_name] = parsed
                    processed += 1
                    clustered += 1
                    progress.set(processed=processed, clustered=clustered, fallback=fallback)
                    progress.message(f"system={system} component={subunit_name} clustering finished")
                    continue

            progress.message(f"system={system} component={subunit_name} clustering failed or produced no usable clusters; "
                             f"returncode={result.returncode}; applying existing fallback; "
                             f"detail={(result.stderr or result.stdout or '')[-1000:].strip()}", "WARN")
            processed += 1
            fallback += 1
            progress.set(processed=processed, clustered=clustered, fallback=fallback)
            for seq_id, _ in sequences:
                all_clusters[system][subunit_name][seq_id] = "0"

    progress.finish()
    with (output_dir / "all_clusters.json").open("w", encoding="utf-8") as handle:
        json.dump(all_clusters, handle, indent=2)

    mapping_frame = pd.DataFrame(mapping_rows)
    mapping_frame.to_csv(output_dir / "sequence_id_mappings.csv", index=False)

    total_systems = len(all_clusters)
    total_sequences = sum(
        len(cluster_map)
        for subunits in all_clusters.values()
        for cluster_map in subunits.values()
    )
    print(
        f"CD-HIT stage complete: {total_systems} system(s), "
        f"{total_sequences} sequence(s), cd-hit={cdhit_available}"
    )


if __name__ == "__main__":
    main()
