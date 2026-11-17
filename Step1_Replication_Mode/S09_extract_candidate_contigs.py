#!/usr/bin/env python3
"""Extract GTDB contigs containing retained Replication mapping candidates."""
import argparse
import gzip
import zlib
from collections import defaultdict
from pathlib import Path

import pandas as pd


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--gtdb_genomes", required=True)
    parser.add_argument("--output", required=True)
    add_progress_argument(parser)
    return parser.parse_args()


def load_genome_paths(genome_root):
    path = Path(genome_root) / "genome_paths.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"GTDB genome path table not found: {path}")
    lookup = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split(maxsplit=1)
            if len(fields) < 2:
                continue
            filename = fields[0]
            accession = filename.removesuffix("_genomic.fna.gz")
            if accession.startswith("GCA_"):
                genome_id = f"GB_{accession}"
            elif accession.startswith("GCF_"):
                genome_id = f"RS_{accession}"
            else:
                continue
            relative = Path(fields[1].strip("/"))
            archive = Path(genome_root) / relative
            if archive.name != filename:
                archive = archive / filename
            previous = lookup.get(genome_id)
            if previous is not None and previous != archive:
                raise ValueError(
                    f"Duplicate GTDB genome path entries for {genome_id}: "
                    f"{previous} and {archive}"
                )
            lookup[genome_id] = archive
    return lookup


def locate_genome(genome_id, genome_paths):
    archive = genome_paths.get(str(genome_id))
    if archive is None:
        raise FileNotFoundError(
            f"GTDB genome ID absent from genome_paths.tsv: {genome_id}"
        )
    if not archive.is_file():
        raise FileNotFoundError(
            f"GTDB genome archive path does not exist for {genome_id}: {archive}"
        )
    return archive


def extract_contigs(archive, needed):
    found = {}
    current_id = None
    sequence = []
    with gzip.open(archive, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current_id in needed:
                    found[current_id] = "".join(sequence)
                current_id = line[1:].split()[0]
                sequence = []
            elif current_id is not None:
                sequence.append(line.strip())
        if current_id in needed:
            found[current_id] = "".join(sequence)
    return found


def main():
    args = parse_args()
    targets = defaultdict(set)
    for csv_path in Path(args.mapping).glob("*/*.csv"):
        frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        for _, row in frame.iterrows():
            genome = str(row.get("genome_id", "")).strip()
            contig = str(row.get("contig_id", "")).strip()
            if genome and contig:
                targets[genome].add(contig)

    genome_paths = load_genome_paths(args.gtdb_genomes)
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    progress = Progress(args.log_label or "Contig-extraction", total=sum(map(len, targets.values())),
                        processed=0, written=0, reused=0, unavailable_targets=0)
    extracted = new = reused = 0
    processed_genomes = 0
    processed_targets = unavailable_targets = 0
    failures = []
    for genome_id, contigs in sorted(targets.items()):
        progress.set(genome=genome_id)
        archive = None
        try:
            archive = locate_genome(genome_id, genome_paths)
            found = extract_contigs(archive, contigs)
        except (FileNotFoundError, OSError, EOFError, zlib.error) as exc:
            failures.append((
                genome_id, "", str(archive) if archive is not None else "", str(exc)
            ))
            processed_genomes += 1
            processed_targets += len(contigs)
            unavailable_targets += len(contigs)
            progress.set(processed=processed_targets, genomes_processed=processed_genomes,
                         unavailable_targets=unavailable_targets)
            continue
        except Exception as exc:
            failures.append((
                genome_id,
                "",
                str(archive) if archive is not None else "",
                f"unexpected {type(exc).__name__}: {exc}",
            ))
            processed_genomes += 1
            processed_targets += len(contigs)
            unavailable_targets += len(contigs)
            progress.set(processed=processed_targets, genomes_processed=processed_genomes,
                         unavailable_targets=unavailable_targets)
            continue
        missing = sorted(contigs - set(found))
        unavailable_targets += len(missing)
        for contig_id in missing:
            failures.append((
                genome_id, contig_id, str(archive),
                "contig absent from genome archive",
            ))
        genome_output = output_root / genome_id
        genome_output.mkdir(parents=True, exist_ok=True)
        for contig_id, sequence in found.items():
            target = genome_output / f"{contig_id}.fasta"
            if target.is_file() and target.stat().st_size > 0:
                extracted += 1
                reused += 1
                continue
            target.write_text(f">{contig_id}\n{sequence}\n", encoding="utf-8")
            extracted += 1
            new += 1
        processed_genomes += 1
        processed_targets += len(contigs)
        progress.set(processed=processed_targets, genomes_processed=processed_genomes,
                     written=new, reused=reused, successful=extracted,
                     unavailable_targets=unavailable_targets)
    progress.finish()
    failure_path = output_root / "candidate_contig_extraction_failures.csv"
    if failures:
        pd.DataFrame(
            failures, columns=["genome_id", "contig_id", "archive_path", "reason"]
        ).to_csv(failure_path, index=False)
        print(
            f"WARNING: skipped {len(failures)} unavailable candidate target(s); "
            f"see {failure_path}"
        )
    else:
        failure_path.unlink(missing_ok=True)
    print(f"Candidate contig extraction complete: {extracted} contigs")


if __name__ == "__main__":
    main()
