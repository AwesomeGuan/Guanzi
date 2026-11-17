#!/usr/bin/env python3
"""Validate Step 1 outputs before a stage completion marker is reused or written."""
import argparse
import csv
import math
from pathlib import Path


FOLDSEEK_COLUMNS = [
    "query_id", "subject_id", "alntmscore", "fident",
    "qcov", "tcov", "evalue", "lddt",
]
CUSTOM_MAPPING_COLUMNS = [
    "defense_system", "query_protein_id", "mining_protein_id",
    "foldseek_alntmscore", "foldseek_sequence_identity",
    "foldseek_query_coverage", "foldseek_target_coverage",
    "foldseek_evalue", "foldseek_lddt",
]
GTDB_COLUMNS = [
    "defense_system", "domain", "genome_id", "contig_id", "afdb50_id",
    "gtdb_protein_id", "gtdb_pident", "gtdb_qcov", "gtdb_tcov",
    "gtdb_evalue", "gtdb_bitscore",
]
FINAL_COLUMNS = [
    "defense_system", "query_protein_id", "afdb50_id", "mining_domain",
    "mining_genome_id", "mining_contig_id", "mining_protein_id",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def nonempty(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def csv_rows(path, columns):
    path = Path(path)
    require(nonempty(path), f"missing or empty CSV: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == columns, f"unexpected columns in {path}: {reader.fieldnames}")
        return list(reader)


def fasta_ids(path):
    ids = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                token = line[1:].split()[0]
                if token:
                    ids.append(token)
    return ids


def validate_prodigal(args):
    root = Path(args.output)
    faa_files = list(root.rglob("*.faa")) if root.is_dir() else []
    require(faa_files, f"no Prodigal FAA files under {root}")
    for faa in faa_files:
        require(nonempty(faa), f"empty Prodigal FAA: {faa}")
        for suffix in (".ffn", ".gff"):
            require(nonempty(faa.with_suffix(suffix)), f"missing companion for {faa}: {suffix}")


def validate_esmfold(args):
    rows = csv_rows(args.status, [
        "protein_id", "faa_path", "pdb_path", "sequence_length", "status", "message"
    ])
    require(rows, f"ESMFold status has no records: {args.status}")
    failures = [row for row in rows if row["status"] in {"out_of_memory", "runtime_error", "error"}]
    require(not failures, f"ESMFold status contains {len(failures)} failed prediction(s)")
    if args.input:
        expected = []
        for faa in Path(args.input).rglob("*.faa"):
            if args.include_parent and faa.parent.name != args.include_parent:
                continue
            expected.extend((str(faa.resolve()), protein_id) for protein_id in fasta_ids(faa))
        actual = [
            (str(Path(row["faa_path"]).resolve()), row["protein_id"])
            for row in rows
        ]
        require(sorted(actual) == sorted(expected), "ESMFold status does not match the current FAA inputs")
    expected = sum(row["status"] in {"success", "skipped"} for row in rows)
    pdb_files = [path for path in Path(args.output).rglob("*.pdb") if nonempty(path)]
    require(len(pdb_files) == expected, f"ESMFold PDB count is {len(pdb_files)}, expected {expected}")


def validate_foldseek_db(args):
    prefix = Path(args.prefix)
    require(nonempty(str(prefix) + ".dbtype"), f"missing FoldSeek database: {prefix}.dbtype")
    require(nonempty(str(prefix) + "_pad.dbtype"), f"missing padded FoldSeek database: {prefix}_pad.dbtype")


def raw_foldseek_rows(root):
    by_system = {}
    skipped_nonfinite = 0
    for system_dir in sorted(path for path in Path(root).iterdir() if path.is_dir()):
        rows = []
        for hit_file in sorted(path for path in system_dir.iterdir() if path.is_file()):
            if not nonempty(hit_file):
                continue
            with hit_file.open(encoding="utf-8", errors="replace") as handle:
                for line_no, line in enumerate(handle, 1):
                    fields = line.rstrip("\n").split("\t")
                    require(len(fields) == 8, f"{hit_file}:{line_no} has {len(fields)} fields, expected 8")
                    try:
                        numeric_values = [float(value) for value in fields[2:]]
                    except ValueError:
                        require(False, f"non-numeric FoldSeek value in {hit_file}:{line_no}")
                    if not all(math.isfinite(value) for value in numeric_values):
                        skipped_nonfinite += 1
                        continue
                    rows.append(dict(zip(FOLDSEEK_COLUMNS, fields)))
        by_system[system_dir.name] = rows
    if skipped_nonfinite:
        print(f"WARNING: skipped {skipped_nonfinite} non-finite FoldSeek row(s)")
    return by_system


def validate_foldseek_search(args):
    queries = list(Path(args.query).rglob("*.pdb"))
    require(queries, f"no Query PDB files under {args.query}")
    raw_foldseek_rows(args.output)
    for query in queries:
        system = query.parent.name
        require((Path(args.output) / system / query.stem).is_file(), f"missing FoldSeek result for {system}/{query.name}")


def validate_coverage(args):
    source = raw_foldseek_rows(args.input)
    output = Path(args.output)
    expected_files = set()
    for system, rows in source.items():
        if not rows:
            continue
        expected_files.add(f"{system}.csv")
        actual = csv_rows(output / f"{system}.csv", FOLDSEEK_COLUMNS)
        expected_count = sum(
            float(row["qcov"]) >= args.min_qcov and float(row["tcov"]) >= args.min_tcov
            for row in rows
        )
        require(len(actual) == expected_count, f"coverage row count mismatch for {system}")
        for row in actual:
            require(float(row["qcov"]) >= args.min_qcov, f"qcov below threshold in {system}")
            require(float(row["tcov"]) >= args.min_tcov, f"tcov below threshold in {system}")
    actual_files = {path.name for path in output.glob("*.csv")} if output.is_dir() else set()
    require(actual_files == expected_files, "coverage output file set does not match FoldSeek systems")


def validate_fident(args):
    source = Path(args.input)
    output = Path(args.output)
    expected_files = set()
    for path in source.glob("*.csv"):
        rows = csv_rows(path, FOLDSEEK_COLUMNS)
        if not rows:
            continue
        expected_files.add(path.name)
        actual = csv_rows(output / path.name, FOLDSEEK_COLUMNS)
        blacklisted = {row["subject_id"] for row in rows if float(row["fident"]) > args.threshold}
        allowed = {row["subject_id"] for row in rows} - blacklisted
        require({row["subject_id"] for row in actual} == allowed, f"fident target set mismatch in {path.name}")
        require(len(actual) == len(allowed), f"duplicate target in {path.name}")
        for row in actual:
            require(float(row["fident"]) <= args.threshold, f"blacklisted target retained in {path.name}")
            best = max(float(item["alntmscore"]) for item in rows if item["subject_id"] == row["subject_id"])
            require(float(row["alntmscore"]) == best, f"non-best TM-score retained in {path.name}")
    actual_files = {path.name for path in output.glob("*.csv")} if output.is_dir() else set()
    require(actual_files == expected_files, "fident output file set does not match nonempty coverage files")


def validate_custom_mapping(args):
    expected = {}
    for path in Path(args.input).glob("*.csv"):
        rows = csv_rows(path, FOLDSEEK_COLUMNS)
        if rows:
            expected[path.name] = len({(row["query_id"], row["subject_id"]) for row in rows})
    actual_files = {path.name for path in Path(args.output).glob("*.csv")} if Path(args.output).is_dir() else set()
    require(actual_files == set(expected), "custom mapping file set mismatch")
    for name, count in expected.items():
        require(len(csv_rows(Path(args.output) / name, CUSTOM_MAPPING_COLUMNS)) == count, f"custom mapping row count mismatch: {name}")


def validate_ids(args):
    expected = {}
    for path in Path(args.input).glob("*.csv"):
        rows = csv_rows(path, FOLDSEEK_COLUMNS)
        expected[path.stem] = {row["subject_id"] for row in rows}
    actual_files = {path.stem for path in Path(args.output).glob("*.txt")} if Path(args.output).is_dir() else set()
    require(actual_files == set(expected), "AFDB50 ID file set mismatch")
    for system, ids in expected.items():
        actual = {line.strip() for line in (Path(args.output) / f"{system}.txt").read_text().splitlines() if line.strip()}
        require(actual == ids, f"AFDB50 ID set mismatch for {system}")


def validate_pdb(args):
    expected = set()
    for txt in Path(args.input).glob("*.txt"):
        for raw_id in (line.strip() for line in txt.read_text().splitlines()):
            if raw_id:
                expected.add((txt.stem, raw_id))
    require(expected, f"no AFDB50 IDs under {args.input}")
    rows = csv_rows(args.status, [
        "defense_system", "protein_id", "status", "message"
    ])
    accounted = set()
    successful = set()
    for row in rows:
        key = (row["defense_system"], row["protein_id"])
        require(key not in accounted, f"duplicate PDB status record: {key}")
        require(key in expected, f"unexpected PDB status record: {key}")
        require(row["status"] in {"success", "failed"}, f"invalid PDB status: {row['status']}")
        if row["status"] == "success":
            require(
                nonempty(Path(args.output) / key[0] / f"{key[1]}.pdb"),
                f"successful PDB status has no file: {key[0]}/{key[1]}",
            )
            successful.add(key)
        accounted.add(key)
    require(accounted == expected, "PDB success/failure accounting does not cover the input ID set")
    actual = {
        (path.parent.name, path.stem)
        for path in Path(args.output).glob("*/*.pdb")
        if nonempty(path)
    }
    require(
        actual == successful,
        f"existing PDB set ({len(actual)}) does not match status success set "
        f"({len(successful)})",
    )
    print(
        f"PDB accounting valid: existing={len(actual)}, success={len(successful)}, "
        f"processed={len(accounted)}, expected={len(expected)}"
    )


def validate_faa(args):
    expected = {
        (path.parent.name, path.stem)
        for path in Path(args.input).glob("*/*.pdb")
        if nonempty(path)
    }
    actual = set()
    for faa in Path(args.output).glob("*/*.faa"):
        require(nonempty(faa) and fasta_ids(faa), f"invalid AFDB50 FAA: {faa}")
        actual.add((faa.parent.name, faa.stem))
    require(actual == expected, "AFDB50 FAA set does not match the current successful PDB set")
    print(f"FAA accounting valid: existing={len(actual)}, expected={len(expected)}")


def validate_query_exclusion(args):
    for faa in Path(args.input).glob("*/*.faa"):
        marker = Path(args.markers) / faa.parent.name / f"{faa.stem}.ok"
        require(nonempty(marker), f"missing query-exclusion task marker: {marker}")
        kept = Path(args.output) / faa.parent.name / faa.name
        if kept.exists():
            require(nonempty(kept) and fasta_ids(kept), f"invalid retained FAA: {kept}")


def validate_gtdb(args):
    root = Path(args.output)
    expected_queries = set()
    for path in Path(args.input).glob("*/*.faa"):
        expected_queries.update((path.parent.name, record_id) for record_id in fasta_ids(path))
    systems = {system for system, _ in expected_queries}
    require(systems, f"no GTDB mapping query FAA files under {args.input}")
    processed = csv_rows(root / "processed_queries.csv", [
        "defense_system", "afdb50_id"
    ])
    actual_queries = {
        (row["defense_system"], row["afdb50_id"]) for row in processed
    }
    require(len(actual_queries) == len(processed), "duplicate GTDB processed-query record")
    require(actual_queries == expected_queries, "GTDB processed-query set is stale")
    for domain in ("bacteria", "archaea"):
        directory = root / domain
        require(directory.is_dir(), f"missing GTDB mapping domain directory: {directory}")
        files = {path.stem: path for path in directory.glob("*.csv")}
        require(set(files) == systems, f"GTDB {domain} output systems do not match query systems")
        for path in files.values():
            for row in csv_rows(path, GTDB_COLUMNS):
                require(float(row["gtdb_pident"]) >= args.min_identity, f"identity below threshold in {path}")
                require(float(row["gtdb_qcov"]) >= args.min_qcov, f"query coverage below threshold in {path}")
                require(float(row["gtdb_tcov"]) >= args.min_tcov, f"target coverage below threshold in {path}")
                require(float(row["gtdb_evalue"]) <= args.evalue, f"E-value above threshold in {path}")


def mapping_targets(root):
    targets = set()
    for path in Path(root).glob("*/*.csv"):
        for row in csv_rows(path, GTDB_COLUMNS):
            if row["genome_id"] and row["contig_id"]:
                targets.add((row["genome_id"], row["contig_id"]))
    return targets


def documented_contig_losses(root):
    path = Path(root) / "candidate_contig_extraction_failures.csv"
    genomes, contigs = set(), set()
    if not path.is_file():
        return genomes, contigs
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames and {"genome_id", "contig_id"} <= set(reader.fieldnames), f"invalid failure report: {path}")
        for row in reader:
            if row["contig_id"]:
                contigs.add((row["genome_id"], row["contig_id"]))
            else:
                genomes.add(row["genome_id"])
    return genomes, contigs


def validate_contigs(args):
    failed_genomes, failed_contigs = documented_contig_losses(args.output)
    for genome, contig in mapping_targets(args.mapping):
        target = Path(args.output) / genome / f"{contig}.fasta"
        require(nonempty(target) or genome in failed_genomes or (genome, contig) in failed_contigs, f"unexplained missing candidate contig: {genome}/{contig}")


def validate_materialized(args):
    failed_genomes, failed_contigs = documented_contig_losses(args.contigs)
    for genome, contig in mapping_targets(args.mapping):
        if genome in failed_genomes or (genome, contig) in failed_contigs:
            continue
        base = Path(args.output) / genome / contig / contig
        for suffix in (".faa", ".ffn", ".gff"):
            require(nonempty(str(base) + suffix), f"missing materialized GTDB file: {base}{suffix}")


def validate_final(args):
    root = Path(args.output)
    require(root.is_dir(), f"missing replication mapping directory: {root}")
    protein_index = set()
    for faa in Path(args.prodigal).rglob("*.faa"):
        genome, contig = faa.parent.parent.name, faa.parent.name
        protein_index.update((genome, contig, protein) for protein in fasta_ids(faa))
    links = {}
    for path in Path(args.foldseek).glob("*.csv"):
        by_target = {}
        for row in csv_rows(path, FOLDSEEK_COLUMNS):
            by_target.setdefault(row["subject_id"], set()).add(row["query_id"])
        links[path.stem] = by_target
    expected = {}
    for path in Path(args.gtdb).glob("*/*.csv"):
        for row in csv_rows(path, GTDB_COLUMNS):
            key = (row["genome_id"], row["contig_id"], row["gtdb_protein_id"])
            if key not in protein_index:
                continue
            for query_id in links.get(row["defense_system"], {}).get(row["afdb50_id"], set()):
                record = (
                    row["defense_system"], query_id, row["afdb50_id"], row["domain"],
                    row["genome_id"], row["contig_id"], row["gtdb_protein_id"],
                )
                expected.setdefault(row["defense_system"], set()).add(record)
    actual_files = {path.stem: path for path in root.glob("*.csv")}
    require(set(actual_files) == set(expected), "replication mapping file set mismatch")
    for system, records in expected.items():
        rows = csv_rows(actual_files[system], FINAL_COLUMNS)
        actual = {
            tuple(row[column] for column in FINAL_COLUMNS)
            for row in rows
        }
        require(len(rows) == len(actual), f"duplicate replication mapping row for {system}")
        require(actual == records, f"replication mapping rows mismatch for {system}")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="stage", required=True)
    def command(name, function, *arguments):
        sub = subs.add_parser(name)
        for flags, kwargs in arguments:
            sub.add_argument(*flags, **kwargs)
        sub.set_defaults(function=function)
    path = lambda name: ((f"--{name}",), {"required": True})
    command("prodigal", validate_prodigal, path("output"))
    command(
        "esmfold", validate_esmfold, path("output"), path("status"),
        (("--input",), {}),
        (("--include-parent",), {}),
    )
    command("foldseek-db", validate_foldseek_db, path("prefix"))
    command("foldseek-search", validate_foldseek_search, path("query"), path("output"))
    command("coverage", validate_coverage, path("input"), path("output"), (("--min-qcov",), {"type": float, "required": True}), (("--min-tcov",), {"type": float, "required": True}))
    command("fident", validate_fident, path("input"), path("output"), (("--threshold",), {"type": float, "required": True}))
    command("custom-mapping", validate_custom_mapping, path("input"), path("output"))
    command("ids", validate_ids, path("input"), path("output"))
    command("pdb", validate_pdb, path("input"), path("output"), path("status"))
    command("faa", validate_faa, path("input"), path("output"))
    command("query-exclusion", validate_query_exclusion, path("input"), path("output"), path("markers"))
    command("gtdb", validate_gtdb, path("input"), path("output"), (("--min-identity",), {"type": float, "required": True}), (("--min-qcov",), {"type": float, "required": True}), (("--min-tcov",), {"type": float, "required": True}), (("--evalue",), {"type": float, "required": True}))
    command("contigs", validate_contigs, path("mapping"), path("output"))
    command("materialized", validate_materialized, path("mapping"), path("contigs"), path("output"))
    command("final", validate_final, path("gtdb"), path("foldseek"), path("prodigal"), path("output"))
    return parser


def main():
    args = build_parser().parse_args()
    try:
        args.function(args)
    except (OSError, ValueError) as exc:
        print(f"OUTPUT VALIDATION FAILED: {exc}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
