#!/usr/bin/env python3
"""Filter Replication instances by the intersection of component-hit contigs.

The recurrence count records contigs with a qualifying homologue of every
component; it does not verify target-gene distances, order or strand.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from step4_utils import *


import sys as _log_sys
from pathlib import Path as _LogPath

_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


EVALUE = "1e-5"
IDENTITY = "50"
QUERY_COVER = "70"
SUBJECT_COVER = "70"
SENSITIVITY_MODE = "sensitive"
MAX_TARGET_SEQS = "0"
CACHE_SCHEMA_VERSION = 3
DOMAINS = ("bacteria", "archaea")
EVIDENCE_COLUMNS = [
    "defense_system",
    "mining_domain",
    "mining_genome_id",
    "mining_contig_id",
    "mining_contig_system_index",
    "homologous_assembly_count",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step4_output", required=True)
    parser.add_argument("--qma", required=True)
    parser.add_argument("--bacteria_db", required=True)
    parser.add_argument("--archaea_db", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--work_dir")
    parser.add_argument("--min_count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--shard_size", type=int, default=128)
    add_progress_argument(parser)
    return parser.parse_args()


def contig(pid):
    bits = str(pid).rsplit("_", 1)
    return bits[0] if len(bits) == 2 and bits[1].isdigit() else str(pid)


def atomic_write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_json(path, value):
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_write_frame(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def diamond_version(diamond):
    result = subprocess.run(
        [diamond, "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    text = (result.stdout or result.stderr).strip()
    return text.splitlines()[0] if text else "unknown"


def database_signature(path):
    candidate = Path(path)
    if not candidate.is_file():
        dmnd = Path(str(candidate) + ".dmnd")
        if dmnd.is_file():
            candidate = dmnd
    if candidate.is_file():
        stat = candidate.stat()
        return {
            "path": str(candidate.resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    return {"path": str(candidate.resolve()), "size": None, "mtime_ns": None}


def cache_manifest(args, diamond):
    return {
        "schema_version": CACHE_SCHEMA_VERSION,
        "diamond_version": diamond_version(diamond),
        "evalue": EVALUE,
        "identity": IDENTITY,
        "query_cover": QUERY_COVER,
        "subject_cover": SUBJECT_COVER,
        "sensitivity_mode": SENSITIVITY_MODE,
        "max_target_seqs": MAX_TARGET_SEQS,
        "databases": {
            "bacteria": database_signature(args.bacteria_db),
            "archaea": database_signature(args.archaea_db),
        },
    }


def prepare_cache(work_dir, manifest):
    work_dir = Path(work_dir)
    manifest_path = work_dir / "run_manifest.json"
    if manifest_path.is_file():
        try:
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            previous = None
        if previous != manifest:
            shutil.rmtree(work_dir / "cache", ignore_errors=True)
            shutil.rmtree(work_dir / "shards", ignore_errors=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_json(manifest_path, manifest)
    return work_dir / "cache"


def read_single_sequence(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one FASTA record, found {len(records)}")
    sequence = str(records[0].seq).strip()
    if not sequence:
        raise ValueError(f"{path}: empty protein sequence")
    return sequence


def query_id(sequence):
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def query_cache_path(cache_root, domain, qid):
    return Path(cache_root) / domain / f"{qid}.contigs"


def read_cached_contigs(path):
    return {
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def shard_id(qids):
    digest = hashlib.sha1("\n".join(qids).encode("utf-8")).hexdigest()
    return digest[:16]


def build_shards(pending, shard_size):
    shards = []
    for domain in DOMAINS:
        qids = pending.get(domain, [])
        for start in range(0, len(qids), shard_size):
            batch = qids[start : start + shard_size]
            shards.append(
                {
                    "domain": domain,
                    "qids": batch,
                    "id": shard_id(batch),
                }
            )
    return shards


def write_shard_fasta(path, qids, query_sequences, domain):
    text = "".join(
        f">{qid}\n{query_sequences[(domain, qid)]}\n"
        for qid in qids
    )
    atomic_write_text(path, text)


def run_shard(shard, diamond, args, query_sequences, cache_root, work_dir, threads):
    domain = shard["domain"]
    qids = shard["qids"]
    shard_root = Path(work_dir) / "shards" / domain
    fasta_path = shard_root / f"{shard['id']}.faa"
    tsv_path = shard_root / f"{shard['id']}.tsv"
    database = args.bacteria_db if domain == "bacteria" else args.archaea_db

    if not tsv_path.is_file():
        write_shard_fasta(fasta_path, qids, query_sequences, domain)
        temporary = tsv_path.with_name(tsv_path.name + ".tmp")
        temporary.parent.mkdir(parents=True, exist_ok=True)
        if temporary.exists():
            temporary.unlink()
        subprocess.run(
            [
                diamond,
                "blastp",
                "-d",
                str(database),
                "-q",
                str(fasta_path),
                "-o",
                str(temporary),
                "-e",
                EVALUE,
                f"--{SENSITIVITY_MODE}",
                "--id",
                IDENTITY,
                "--query-cover",
                QUERY_COVER,
                "--subject-cover",
                SUBJECT_COVER,
                "--max-target-seqs",
                MAX_TARGET_SEQS,
                "--threads",
                str(threads),
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "--quiet",
            ],
            check=True,
        )
        os.replace(temporary, tsv_path)

    hits = {qid: set() for qid in qids}
    with tsv_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{tsv_path}: malformed DIAMOND output line: {line.rstrip()}")
            qid, subject = fields[0], fields[1]
            if qid in hits:
                hits[qid].add(contig(subject))

    for qid in qids:
        text = "".join(f"{value}\n" for value in sorted(hits[qid]))
        atomic_write_text(query_cache_path(cache_root, domain, qid), text)

    fasta_path.unlink(missing_ok=True)
    tsv_path.unlink(missing_ok=True)
    return len(qids)


def parse_instances(args, assembly):
    qma = Path(args.qma)
    query_sequences = {}
    instances = {}
    if "mining_domain" not in assembly.columns:
        raise ValueError(
            "Replication frequency input is missing required column: mining_domain"
        )
    domains = assembly["mining_domain"].astype(str).str.strip()
    invalid = ~domains.isin(DOMAINS)
    if invalid.any():
        counts = domains[invalid].value_counts(dropna=False).to_dict()
        details = ", ".join(
            f"{value or '<empty>'}={count}" for value, count in sorted(counts.items())
        )
        raise ValueError(
            f"Replication frequency input has {int(invalid.sum())} instance(s) "
            f"with invalid mining_domain; expected bacteria or archaea: {details}"
        )
    for index, row in assembly.iterrows():
        domain = domains.at[index]
        key = tuple(str(value) for value in instance_key(row))
        leaf = qma / key[0] / f"{key[1]}__{key[2]}" / key[3] / "mining_assembly"
        qids = []
        for fasta_path in sorted(leaf.glob("*.faa")):
            sequence = read_single_sequence(fasta_path)
            qid = query_id(sequence)
            qids.append(qid)
            query_sequences[(domain, qid)] = sequence
        instances[index] = {"key": key, "domain": domain, "qids": qids}
    return instances, query_sequences


def main():
    args = parse_args()
    if args.min_count < 1:
        raise ValueError("--min_count must be positive")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if args.shard_size < 1:
        raise ValueError("--shard_size must be positive")

    diamond = shutil.which("diamond")
    if not diamond:
        raise FileNotFoundError("DIAMOND executable not found")

    step4_output = Path(args.step4_output)
    output_dir = Path(args.output)
    work_dir = Path(args.work_dir) if args.work_dir else output_dir / "work"
    assembly = pd.concat(
        [
            read_csv(step4_output / f"assembly_{status}.csv")
            for status in ("qualified", "disqualified")
        ],
        ignore_index=True,
    )

    cache_root = prepare_cache(work_dir, cache_manifest(args, diamond))
    instances, query_sequences = parse_instances(args, assembly)

    cached_qids = {
        query_key
        for query_key in query_sequences
        if query_cache_path(cache_root, query_key[0], query_key[1]).is_file()
    }
    pending = defaultdict(list)
    for domain, qid in sorted(query_sequences):
        if (domain, qid) not in cached_qids:
            pending[domain].append(qid)

    shards = build_shards(pending, args.shard_size)
    workers = max(1, min(args.workers, len(shards) or 1))
    threads_per_search = args.threads
    total_threads = threads_per_search * workers
    progress = Progress(
        args.log_label or "Homologous-frequency",
        total=len(query_sequences),
        processed=len(cached_qids),
        cache_reused=len(cached_qids),
        completed_shards=0,
        shards=len(shards),
        retained=0,
        workers=workers,
        threads_per_search=threads_per_search,
        total_threads=total_threads,
    )
    print(
        "[S09] INFO "
        f"unique_queries={len(query_sequences)} "
        f"cache_hit_queries={len(cached_qids)} "
        f"pending_queries={len(query_sequences) - len(cached_qids)} "
        f"shards={len(shards)} workers={workers} "
        f"threads_per_search={threads_per_search} "
        f"total_threads={total_threads}"
    )

    completed_queries = 0
    completed_shards = 0
    errors = []
    if shards:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    run_shard,
                    shard,
                    diamond,
                    args,
                    query_sequences,
                    cache_root,
                    work_dir,
                    threads_per_search,
                ): shard
                for shard in shards
            }
            for future in as_completed(futures):
                shard = futures[future]
                try:
                    completed_queries += future.result()
                except Exception as exc:
                    errors.append(f"{shard['domain']}/{shard['id']}: {exc}")
                    continue
                completed_shards += 1
                progress.set(
                    processed=len(cached_qids) + completed_queries,
                    completed_shards=completed_shards,
                )

    if errors:
        progress.close()
        raise RuntimeError("DIAMOND shard failure(s): " + "; ".join(errors))

    missing = [
        query_key
        for query_key in query_sequences
        if not query_cache_path(cache_root, query_key[0], query_key[1]).is_file()
    ]
    if missing:
        progress.close()
        raise RuntimeError(f"{len(missing)} query result(s) remain unavailable after DIAMOND")

    query_hits = {
        query_key: read_cached_contigs(
            query_cache_path(cache_root, query_key[0], query_key[1])
        )
        for query_key in query_sequences
    }

    evidence = {}
    passed = set()
    for index, spec in instances.items():
        hit_sets = [
            query_hits[(spec["domain"], qid)]
            for qid in spec["qids"]
        ]
        count = len(set.intersection(*hit_sets)) if hit_sets and all(hit_sets) else 0
        key = spec["key"]
        evidence[index] = {
            "defense_system": key[0],
            "mining_domain": spec["domain"],
            "mining_genome_id": key[1],
            "mining_contig_id": key[2],
            "mining_contig_system_index": key[3],
            "homologous_assembly_count": count,
        }
        if count >= args.min_count:
            passed.add(key)

    rows = [
        evidence[index]
        for index in assembly.index
        if index in evidence
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_frame(
        pd.DataFrame(rows, columns=EVIDENCE_COLUMNS),
        output_dir / "homologous_frequency_all.csv",
    )
    for status in ("qualified", "disqualified"):
        for kind in ("mapping", "assembly", "comparison"):
            path = step4_output / f"{kind}_{status}.csv"
            atomic_write_frame(subset_keys(read_csv(path), passed), path)

    shutil.rmtree(work_dir / "shards", ignore_errors=True)
    progress.set(retained=len(passed))
    progress.finish()
    print(
        f"replication frequency complete: {len(passed)} instance(s) retained "
        f"from {len(rows)}; searched_queries={len(query_sequences) - len(cached_qids)}"
    )


if __name__ == "__main__":
    main()
