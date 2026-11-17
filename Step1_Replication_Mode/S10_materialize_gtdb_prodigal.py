#!/usr/bin/env python3
"""Materialize GTDB gene calls as per-contig FAA, FFN and GFF files."""
import argparse
import re
import shutil
from collections import defaultdict
from pathlib import Path

import pandas as pd


HEADER_RE = re.compile(
    r"^(?P<protein_id>\S+)\s+#\s+(?P<start>\d+)\s+#\s+"
    r"(?P<end>\d+)\s+#\s+(?P<strand>-?1)\s+#\s*(?P<attrs>.*)$"
)


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtdb_mapping", required=True)
    parser.add_argument("--protein_faa_reps", required=True)
    parser.add_argument("--candidate_genomes", required=True)
    parser.add_argument("--output", required=True)
    add_progress_argument(parser)
    return parser.parse_args()


def fasta_records(path):
    description = None
    sequence = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if description is not None:
                    yield description, "".join(sequence)
                description = line[1:]
                sequence = []
            elif description is not None:
                sequence.append(line)
        if description is not None:
            yield description, "".join(sequence)


def read_single_fasta(path):
    records = list(fasta_records(path))
    if len(records) != 1:
        raise ValueError(f"Expected one contig in {path}, found {len(records)}")
    return records[0][0].split()[0], records[0][1].upper()


def contig_from_protein(protein_id):
    parts = protein_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    raise ValueError(f"Cannot derive contig ID from GTDB protein ID: {protein_id}")


def parse_description(description):
    match = HEADER_RE.match(description)
    if match is None:
        raise ValueError(f"Unparseable GTDB FAA header: {description}")
    record = match.groupdict()
    record["start"] = int(record["start"])
    record["end"] = int(record["end"])
    return record


DNA_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def reverse_complement(sequence):
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def gff_attributes(protein_id, raw_attrs):
    retained = []
    for attribute in raw_attrs.strip().strip(";").split(";"):
        attribute = attribute.strip()
        if attribute and not attribute.startswith("ID="):
            retained.append(attribute)
    return ";".join([f"ID={protein_id}", *retained]) + ";"


def load_targets(mapping_root):
    targets = set()
    target_columns = ["domain", "genome_id", "contig_id"]
    for path in sorted(Path(mapping_root).glob("*/*.csv")):
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if frame.empty:
            continue
        missing = set(target_columns) - set(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        for row in frame[target_columns].itertuples(index=False, name=None):
            values = tuple(str(value).strip() for value in row)
            if all(values):
                targets.add(values)
    if not targets:
        raise ValueError(f"No GTDB candidate contigs found under {mapping_root}")
    return targets


def load_documented_upstream_losses(candidate_root):
    """Return genome-wide and contig-specific losses recorded by S09."""
    failure_path = Path(candidate_root) / "candidate_contig_extraction_failures.csv"
    if not failure_path.is_file():
        return set(), set()
    frame = pd.read_csv(failure_path, dtype=str, keep_default_na=False)
    required = {"genome_id", "contig_id"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"{failure_path} is missing columns: {sorted(missing)}"
        )
    failed_genomes = set()
    failed_contigs = set()
    for genome_id, contig_id in frame[["genome_id", "contig_id"]].itertuples(
            index=False, name=None):
        genome_id = str(genome_id).strip()
        contig_id = str(contig_id).strip()
        if not genome_id:
            continue
        if contig_id:
            failed_contigs.add((genome_id, contig_id))
        else:
            failed_genomes.add(genome_id)
    return failed_genomes, failed_contigs


def load_source_records(source_faa, needed_contigs):
    records = defaultdict(list)
    for description, sequence in fasta_records(source_faa):
        protein_id = description.split()[0]
        contig_id = contig_from_protein(protein_id)
        if contig_id not in needed_contigs:
            continue
        parsed = parse_description(description)
        if parsed["protein_id"] != protein_id:
            raise ValueError(f"FAA header ID mismatch in {source_faa}: {description}")
        parsed["description"] = description
        parsed["sequence"] = sequence.replace("*", "")
        records[contig_id].append(parsed)
    return records


def write_contig_output(genome_id, contig_id, dna, records, output_root):
    if not records:
        raise ValueError(f"No GTDB proteins found for {genome_id}/{contig_id}")
    target = output_root / genome_id / contig_id
    temporary = target.with_name(f".{target.name}.tmp")
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True)
    faa_path = temporary / f"{contig_id}.faa"
    ffn_path = temporary / f"{contig_id}.ffn"
    gff_path = temporary / f"{contig_id}.gff"

    seen_ids = set()
    with faa_path.open("w", encoding="utf-8") as faa, \
            ffn_path.open("w", encoding="utf-8") as ffn, \
            gff_path.open("w", encoding="utf-8") as gff:
        gff.write("##gff-version 3\n")
        gff.write(f"##sequence-region {contig_id} 1 {len(dna)}\n")
        gff.write("##source=GTDB original Prodigal protein predictions\n")
        for record in records:
            protein_id = record["protein_id"]
            if protein_id in seen_ids:
                raise ValueError(f"Duplicate GTDB protein ID: {protein_id}")
            seen_ids.add(protein_id)
            start, end = record["start"], record["end"]
            if start < 1 or end < start or end > len(dna):
                raise ValueError(
                    f"Coordinate outside contig for {protein_id}: {start}-{end}/{len(dna)}"
                )
            cds = dna[start - 1:end]
            if record["strand"] == "-1":
                cds = reverse_complement(cds)
                strand = "-"
            elif record["strand"] == "1":
                strand = "+"
            else:
                raise ValueError(f"Invalid strand for {protein_id}: {record['strand']}")

            faa.write(f">{record['description']}\n{record['sequence']}\n")
            ffn.write(f">{record['description']}\n{cds}\n")
            attributes = gff_attributes(protein_id, record["attrs"])
            gff.write("\t".join([
                contig_id, "GTDB_Prodigal", "CDS", str(start), str(end), ".",
                strand, "0", attributes,
            ]) + "\n")

    if target.exists():
        shutil.rmtree(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary.replace(target)
    return len(records)


def main():
    args = parse_args()
    targets = load_targets(args.gtdb_mapping)
    failed_genomes, failed_contigs = load_documented_upstream_losses(
        args.candidate_genomes
    )
    by_genome = defaultdict(set)
    for domain, genome_id, contig_id in targets:
        by_genome[(domain, genome_id)].add(contig_id)

    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    progress = Progress(args.log_label or "GTDB-materialization", total=len(targets),
                        processed=0, written_contigs=0, reused_contigs=0, upstream_losses=0, failures=0)
    processed_genomes = reused_contigs = new_contigs = 0
    processed_contigs = 0
    failures = []
    contig_count = protein_count = upstream_loss_count = 0
    for (domain, genome_id), contigs in sorted(by_genome.items()):
        progress.set(domain=domain, genome=genome_id)
        source_faa = Path(args.protein_faa_reps) / domain / f"{genome_id}_protein.faa"
        if not source_faa.is_file():
            failures.append((domain, genome_id, "", f"source FAA missing: {source_faa}"))
            processed_genomes += 1
            processed_contigs += len(contigs)
            progress.set(processed=processed_contigs, genomes_processed=processed_genomes, failures=len(failures))
            continue
        try:
            source_records = load_source_records(source_faa, contigs)
        except Exception as exc:
            failures.append((domain, genome_id, "", str(exc)))
            processed_genomes += 1
            processed_contigs += len(contigs)
            progress.set(processed=processed_contigs, genomes_processed=processed_genomes, failures=len(failures))
            continue
        for contig_id in sorted(contigs):
            candidate = Path(args.candidate_genomes) / genome_id / f"{contig_id}.fasta"
            if not candidate.is_file():
                if genome_id in failed_genomes or (genome_id, contig_id) in failed_contigs:
                    upstream_loss_count += 1
                    continue
                failures.append((
                    domain, genome_id, contig_id,
                    f"candidate contig missing without an S09 loss record: {candidate}",
                ))
                continue
            base = output_root / genome_id / contig_id / contig_id
            if all(
                Path(str(base) + suffix).is_file()
                and Path(str(base) + suffix).stat().st_size > 0
                for suffix in (".faa", ".ffn", ".gff")
            ):
                contig_count += 1
                protein_count += len(source_records.get(contig_id, []))
                reused_contigs += 1
                continue
            try:
                fasta_id, dna = read_single_fasta(candidate)
                if fasta_id != contig_id:
                    raise ValueError(
                        f"Candidate FASTA ID differs from expected contig: {fasta_id}"
                    )
                written = write_contig_output(
                    genome_id, contig_id, dna, source_records.get(contig_id, []),
                    output_root,
                )
                contig_count += 1
                protein_count += written
                new_contigs += 1
            except Exception as exc:
                failures.append((domain, genome_id, contig_id, str(exc)))

        processed_genomes += 1
        processed_contigs += len(contigs)
        progress.set(processed=processed_contigs, genomes_processed=processed_genomes,
                     written_contigs=new_contigs, reused_contigs=reused_contigs,
                     successful_contigs=contig_count, proteins=protein_count,
                     upstream_losses=upstream_loss_count, failures=len(failures))
    progress.finish()
    failure_path = output_root / "gtdb_prodigal_materialization_failures.csv"
    if failures:
        pd.DataFrame(
            failures, columns=["domain", "genome_id", "contig_id", "reason"]
        ).to_csv(failure_path, index=False)
        raise ValueError(
            f"{len(failures)} GTDB gene-call materialization failure(s); see {failure_path}"
        )
    failure_path.unlink(missing_ok=True)
    print(
        f"GTDB Prodigal materialization complete: {contig_count} contig(s), "
        f"{protein_count} protein(s), {upstream_loss_count} upstream loss(es) skipped"
    )


if __name__ == "__main__":
    main()
