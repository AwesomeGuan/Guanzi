#!/usr/bin/env python3
"""Incrementally map Query-excluded AFDB50 proteins to GTDB homologs."""
import argparse
import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

COLUMNS = [
    "defense_system", "domain", "genome_id", "contig_id", "afdb50_id",
    "gtdb_protein_id", "gtdb_pident", "gtdb_qcov", "gtdb_tcov",
    "gtdb_evalue", "gtdb_bitscore",
]
PROCESSED_COLUMNS = ["defense_system", "afdb50_id"]


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_fasta(path):
    record_id = None
    sequence = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(sequence)
                record_id = line[1:].split()[0]
                sequence = []
            elif record_id is not None:
                sequence.append(line)
        if record_id is not None:
            yield record_id, "".join(sequence)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query_dir", required=True)
    parser.add_argument("--diamond_db_dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--min_identity", type=float, default=90.0)
    parser.add_argument("--min_qcov", type=float, default=90.0)
    parser.add_argument("--min_tcov", type=float, default=90.0)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--max_target_seqs", type=int, default=0)
    parser.add_argument("--reset", action="store_true")
    add_progress_argument(parser)
    return parser.parse_args()


def collect_queries(query_root):
    records = {}
    for system_dir in sorted(Path(query_root).iterdir()):
        if not system_dir.is_dir():
            continue
        for faa in sorted(system_dir.glob("*.faa")):
            for record_id, sequence in parse_fasta(faa):
                key = (system_dir.name, record_id)
                if key in records and records[key] != sequence:
                    raise ValueError(f"Conflicting AFDB50 sequence for {key}")
                records[key] = sequence.replace("*", "")
    return records


def load_processed(path):
    path = Path(path)
    if not path.is_file():
        return set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != PROCESSED_COLUMNS:
            return set()
        return {
            (row["defense_system"], row["afdb50_id"])
            for row in reader if row["defense_system"] and row["afdb50_id"]
        }


def write_csv_atomic(path, columns, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent,
        prefix=f".{path.name}.", delete=False,
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def write_queries(records, keys, target):
    lookup = {}
    with Path(target).open("w", encoding="utf-8") as handle:
        for counter, key in enumerate(sorted(keys)):
            token = f"Q{counter:09d}"
            lookup[token] = key
            handle.write(f">{token}\n{records[key]}\n")
    return lookup


def load_locations(path, needed):
    locations = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein_id = row.get("gtdb_protein_id", "")
            if protein_id in needed:
                locations[protein_id].append(
                    (row.get("genome_id", ""), row.get("contig_id", ""))
                )
    return locations


def existing_rows(path):
    path = Path(path)
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != COLUMNS:
            raise ValueError(f"Unexpected GTDB columns in {path}")
        return list(reader)


def append_domain_rows(domain_output, systems, additions):
    domain_output.mkdir(parents=True, exist_ok=True)
    for system in sorted(systems):
        path = domain_output / f"{system}.csv"
        rows = existing_rows(path)
        rows.extend(additions.get(system, []))
        unique = {
            tuple(row[column] for column in COLUMNS): row
            for row in rows
        }
        ordered = [unique[key] for key in sorted(unique)]
        write_csv_atomic(path, COLUMNS, ordered)


def search_domain(domain, query_faa, query_lookup, db_root, output_root, args):
    db = db_root / f"gtdb_r220_{domain}"
    index = db_root / f"gtdb_r220_{domain}_protein_locations.tsv.gz"
    if not db.with_suffix(".dmnd").is_file() or not index.is_file():
        raise FileNotFoundError(f"Missing {domain} GTDB DIAMOND resources under {db_root}")
    with tempfile.TemporaryDirectory(prefix=f"rep_{domain}_") as temp_dir:
        hits_path = Path(temp_dir) / "hits.tsv"
        command = [
            "diamond", "blastp", "-d", str(db), "-q", str(query_faa),
            "-o", str(hits_path), "--id", str(args.min_identity),
            "--query-cover", str(args.min_qcov), "--subject-cover", str(args.min_tcov),
            "--evalue", str(args.evalue), "--max-target-seqs", str(args.max_target_seqs),
            "--threads", str(args.threads), "--outfmt", "6", "qseqid", "sseqid",
            "pident", "qcovhsp", "scovhsp", "evalue", "bitscore", "--quiet",
        ]
        args.progress.message(f"domain={domain} phase=DIAMOND-search queries={len(query_lookup)}")
        args.progress.set(domain=domain, phase="DIAMOND-search")
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"DIAMOND GTDB mapping failed for {domain}")
        raw_hits = []
        if hits_path.is_file() and hits_path.stat().st_size:
            with hits_path.open(encoding="utf-8") as handle:
                raw_hits = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    args.progress.message(f"domain={domain} phase=location-index")
    args.progress.set(domain=domain, phase="location-index")
    locations = load_locations(index, {row[1] for row in raw_hits})
    rows_by_system = defaultdict(list)
    unresolved = set()
    for row in raw_hits:
        token, gtdb_id, pident, qcov, scov, evalue, bitscore = row
        system, afdb50_id = query_lookup[token]
        matched_locations = locations.get(gtdb_id, [])
        if not matched_locations:
            unresolved.add(gtdb_id)
            continue
        for genome_id, contig_id in matched_locations:
            rows_by_system[system].append({
                "defense_system": system, "domain": domain,
                "genome_id": genome_id, "contig_id": contig_id,
                "afdb50_id": afdb50_id, "gtdb_protein_id": gtdb_id,
                "gtdb_pident": pident, "gtdb_qcov": qcov,
                "gtdb_tcov": scov, "gtdb_evalue": evalue,
                "gtdb_bitscore": bitscore,
            })
    if unresolved:
        raise ValueError(f"{domain}: {len(unresolved)} DIAMOND IDs absent from location index")
    args.progress.message(f"domain={domain} phase=write-results")
    args.progress.set(domain=domain, phase="write-results")
    append_domain_rows(
        output_root / domain,
        {system for system, _ in query_lookup.values()},
        rows_by_system,
    )
    print(f"{domain}: added {sum(map(len, rows_by_system.values()))} GTDB locations")


def main():
    args = parse_args()
    for name in ("min_identity", "min_qcov", "min_tcov"):
        value = getattr(args, name)
        if not 0 <= value <= 100:
            raise ValueError(f"--{name} must be between 0 and 100: {value}")
    if args.evalue < 0 or args.max_target_seqs < 0:
        raise ValueError("E-value and max_target_seqs must be non-negative")
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = output_root / "processed_queries.csv"
    args.progress = Progress(args.log_label or "GTDB-mapping", phase="collect-queries")
    records = collect_queries(args.query_dir)
    if not records:
        raise ValueError(f"No AFDB50 FAA records found under {args.query_dir}")
    current = set(records)
    processed = load_processed(manifest)
    processed_systems = {system for system, _ in processed}
    outputs_complete = True
    for domain in ("bacteria", "archaea"):
        for system in processed_systems:
            path = output_root / domain / f"{system}.csv"
            try:
                with path.open(newline="", encoding="utf-8-sig") as handle:
                    outputs_complete &= csv.DictReader(handle).fieldnames == COLUMNS
            except OSError:
                outputs_complete = False
    if args.reset or not processed.issubset(current) or not outputs_complete:
        for domain in ("bacteria", "archaea"):
            shutil.rmtree(output_root / domain, ignore_errors=True)
        processed = set()
    pending = current - processed
    args.progress.set(total_queries=len(current), reused_queries=len(processed), new_queries=len(pending))
    args.progress.message(f"total_queries={len(current)} reused_queries={len(processed)} new_queries={len(pending)}")
    if pending:
        with tempfile.TemporaryDirectory(prefix="rep_gtdb_queries_") as temp_dir:
            query_faa = Path(temp_dir) / "queries.faa"
            query_lookup = write_queries(records, pending, query_faa)
            for domain in ("bacteria", "archaea"):
                search_domain(
                    domain, query_faa, query_lookup, Path(args.diamond_db_dir),
                    output_root, args,
                )
        processed.update(pending)
        write_csv_atomic(
            manifest, PROCESSED_COLUMNS,
            [
                {"defense_system": system, "afdb50_id": protein_id}
                for system, protein_id in sorted(processed)
            ],
        )
    args.progress.set(phase="complete")
    args.progress.finish()
    print(
        f"GTDB mapping: {len(pending)} new query protein(s); "
        f"identity >= {args.min_identity}%, qcov >= {args.min_qcov}%, "
        f"tcov >= {args.min_tcov}%, evalue <= {args.evalue}, "
        f"max_target_seqs={args.max_target_seqs}"
    )


if __name__ == "__main__":
    main()
