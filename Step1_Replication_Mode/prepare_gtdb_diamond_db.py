#!/usr/bin/env python3
"""Build domain-specific GTDB DIAMOND databases and protein-location indexes."""
import argparse
import gzip
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein_faa_reps", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def contig_from_protein(protein_id):
    parts = str(protein_id).rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    raise ValueError(f"Cannot derive contig ID from GTDB protein ID: {protein_id}")


def build_domain(domain, source_root, output_root, threads, force):
    source_dir = source_root / domain
    if not source_dir.is_dir():
        raise FileNotFoundError(f"GTDB {domain} FAA directory not found: {source_dir}")

    db_prefix = output_root / f"gtdb_r220_{domain}"
    database = db_prefix.with_suffix(".dmnd")
    index_path = output_root / f"gtdb_r220_{domain}_protein_locations.tsv.gz"
    if database.is_file() and index_path.is_file() and not force:
        print(f"{domain}: database and location index already exist, skipped")
        return

    faa_files = sorted(source_dir.glob("*_protein.faa"))
    if not faa_files:
        raise FileNotFoundError(f"No *_protein.faa files found under {source_dir}")

    merged = output_root / f".gtdb_r220_{domain}.faa"
    index_tmp = output_root / f".gtdb_r220_{domain}_protein_locations.tsv.gz"
    records = 0
    with merged.open("wb") as merged_handle, gzip.open(index_tmp, "wt", encoding="utf-8") as index:
        index.write("gtdb_protein_id\tgenome_id\tcontig_id\n")
        for number, faa in enumerate(faa_files, start=1):
            genome_id = faa.name.removesuffix("_protein.faa")
            last_line = b""
            with faa.open("rb") as source:
                for line in source:
                    merged_handle.write(line)
                    last_line = line
                    if not line.startswith(b">"):
                        continue
                    protein_id = line[1:].split()[0].decode("utf-8", errors="replace")
                    contig_id = contig_from_protein(protein_id)
                    index.write(f"{protein_id}\t{genome_id}\t{contig_id}\n")
                    records += 1
            if last_line and not last_line.endswith(b"\n"):
                merged_handle.write(b"\n")
            if number % 10000 == 0:
                print(f"{domain}: indexed {number}/{len(faa_files)} FAA files")

    command = [
        "diamond", "makedb", "--in", str(merged), "--db", str(db_prefix),
        "--threads", str(threads),
    ]
    completed = subprocess.run(command, check=False)
    merged.unlink(missing_ok=True)
    if completed.returncode != 0 or not database.is_file():
        index_tmp.unlink(missing_ok=True)
        raise RuntimeError(f"DIAMOND database creation failed for {domain}")
    index_tmp.replace(index_path)
    print(f"{domain}: {records} proteins indexed -> {database}")


def main():
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    source_root = Path(args.protein_faa_reps)
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    for domain in ("bacteria", "archaea"):
        build_domain(domain, source_root, output_root, args.threads, args.force)


if __name__ == "__main__":
    main()
