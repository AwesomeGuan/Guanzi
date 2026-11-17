#!/usr/bin/env python3
"""Assemble remaining candidates with homologues of their exact Query instance."""

import argparse
import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.query_location import get_instance, read_query_location
from step3_assembly_utils import (
    assembly_records_to_dataframe,
    check_completeness,
    load_rules_from_dir,
    make_mapping_rows,
    make_assembly_record,
    read_csv_or_empty,
    step3_mapping_columns,
    write_csv_with_columns,
)


import sys as _log_sys
from pathlib import Path as _LogPath
_log_sys.path.insert(0, str(_LogPath(__file__).resolve().parents[1]))
from utils.progress import Progress, add_progress_argument


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--query_location", required=True)
    parser.add_argument("--prodigal", required=True)
    parser.add_argument("--rules_dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--diamond_bin", required=True)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--min_identity", type=float, default=50.0)
    parser.add_argument("--min_qcov", type=float, default=70.0)
    parser.add_argument("--min_tcov", type=float, default=70.0)
    parser.add_argument("--evalue", type=float, default=1e-5)
    add_progress_argument(parser)
    return parser.parse_args()


def load_query_tables(root):
    tables = {}
    for path in sorted(Path(root).glob("*.csv")):
        tables[path.stem] = read_query_location(path, path.stem)
    return tables


def contig_faa(root, genome_id, contig_id):
    return Path(root) / genome_id / contig_id / f"{contig_id}.faa"


def database_prefix(cache_dir, faa):
    digest = hashlib.sha1(str(faa.resolve()).encode("utf-8")).hexdigest()[:12]
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", faa.stem)
    return cache_dir / f"{safe}_{digest}"


def ensure_database(diamond, faa, cache_dir):
    prefix = database_prefix(cache_dir, faa)
    database = Path(f"{prefix}.dmnd")
    if not database.is_file() or database.stat().st_mtime < faa.stat().st_mtime:
        cache_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [diamond, "makedb", "--in", str(faa), "--db", str(prefix), "--quiet"],
            check=True,
        )
    return database


def write_partner_queries(instance, own_component, target):
    query_map = {}
    with target.open("w", encoding="utf-8") as handle:
        serial = 0
        for _, member in instance.iterrows():
            component = str(member["component"]).strip()
            if component == own_component:
                continue
            protein_id = str(member["protein_id"]).strip()
            sequence = str(member["protein_seq"]).replace("*", "").strip()
            if not protein_id or not sequence:
                continue
            serial += 1
            synthetic_id = f"Q{serial:06d}"
            query_map[synthetic_id] = (component, protein_id)
            handle.write(f">{synthetic_id}\n{sequence}\n")
    return query_map


def run_search(args, database, query_faa, result_tsv):
    subprocess.run(
        [
            args.diamond_bin, "blastp", "--db", str(database),
            "--query", str(query_faa), "--out", str(result_tsv),
            "--sensitive",
            "--threads", str(args.threads), "--id", str(args.min_identity),
            "--query-cover", str(args.min_qcov),
            "--subject-cover", str(args.min_tcov), "--evalue", str(args.evalue),
            "--max-target-seqs", "0", "--outfmt", "6", "qseqid", "sseqid",
            "pident", "qcovhsp", "scovhsp", "evalue", "bitscore", "--quiet",
        ],
        check=True,
    )


def read_hits(path, query_map):
    columns = ["synthetic_id", "target_gene_id", "pident", "qcov", "tcov", "evalue", "bitscore"]
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns + ["subunit_name", "partner_query_id"])
    frame = pd.read_csv(path, sep="\t", names=columns, dtype=str, keep_default_na=False)
    frame["subunit_name"] = frame["synthetic_id"].map(lambda value: query_map[value][0])
    frame["partner_query_id"] = frame["synthetic_id"].map(lambda value: query_map[value][1])
    return frame


def best_edges(hits, excluded_target):
    """Return component -> ordered distinct target IDs, excluding the candidate itself."""
    usable = hits[hits["target_gene_id"].astype(str) != str(excluded_target)].copy()
    if usable.empty:
        return {}, usable
    usable["_bitscore"] = pd.to_numeric(usable["bitscore"], errors="coerce").fillna(-1)
    usable["_evalue"] = pd.to_numeric(usable["evalue"], errors="coerce").fillna(float("inf"))
    usable = usable.sort_values(
        ["_bitscore", "_evalue", "pident", "qcov", "tcov"],
        ascending=[False, True, False, False, False],
        kind="stable",
    )
    edges = {}
    for component, group in usable.groupby("subunit_name", sort=False):
        edges[str(component)] = list(dict.fromkeys(group["target_gene_id"].astype(str)))
    return edges, usable.drop(columns=["_bitscore", "_evalue"])


def one_to_one_assignment(expected_components, edges):
    """Find a target-disjoint assignment using augmenting paths."""
    target_to_component = {}

    def augment(component, visited):
        for target in edges.get(component, []):
            if target in visited:
                continue
            visited.add(target)
            owner = target_to_component.get(target)
            if owner is None or augment(owner, visited):
                target_to_component[target] = component
                return True
        return False

    for component in sorted(expected_components, key=lambda value: (len(edges.get(value, [])), value)):
        if not augment(component, set()):
            return None
    return {component: target for target, component in target_to_component.items()}


def selected_targets(hits, assignment):
    """Keep mandatory assignment targets and assign every other target to its best component."""
    selected = dict(assignment)
    assigned_targets = set(selected.values())
    for target, group in hits.groupby("target_gene_id", sort=False):
        target = str(target)
        if target in assigned_targets:
            continue
        best = group.iloc[0]
        selected[f"{str(best['subunit_name'])}\0{target}"] = target
    target_components = {}
    for key, target in selected.items():
        component = key.split("\0", 1)[0]
        target_components[str(target)] = component
    return target_components


def failure_reason(expected, edges, assignment):
    missing = sorted(component for component in expected if not edges.get(component))
    if missing:
        return "missing_component_homologue:" + "|".join(missing)
    if assignment is None:
        return "component_target_one_to_one_conflict"
    return "query_instance_not_rule_complete"


def main():
    args = parse_args()
    diamond = Path(args.diamond_bin)
    if not diamond.is_file():
        raise FileNotFoundError(f"DIAMOND executable not found: {diamond}")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    database_cache = output_dir / "diamond_contig_db"
    frame = read_csv_or_empty(args.input)
    mapping_columns = list(frame.columns)
    success_columns = step3_mapping_columns(frame)
    failure_columns = mapping_columns
    tables = load_query_tables(args.query_location)
    rules = load_rules_from_dir(args.rules_dir)

    required = {
        "defense_system", "query_protein_id", "mining_protein_id", "mining_component",
        "mining_genome_id", "mining_contig_id", "query_genome_id", "query_contig_id",
        "query_contig_system_index",
    }
    missing_columns = required - set(frame.columns)
    if missing_columns:
        raise ValueError(f"missing mapping columns: {sorted(missing_columns)}")

    successes = []
    failures = []
    assembly_records = []
    search_cache = {}
    progress = Progress(args.log_label or "Query-guided-assembly", processed_groups=0, searches=0, cache_reused=0, assembled=0, unmatched=0)
    search_count = cache_reused = 0
    progress.message(f"mode=sensitive identity>={args.min_identity}% qcov>={args.min_qcov}% tcov>={args.min_tcov}% evalue<={args.evalue}")

    with tempfile.TemporaryDirectory(prefix="step3_query_guided_") as temp_name:
        temp = Path(temp_name)
        group_columns = [
            "defense_system", "mining_genome_id", "mining_contig_id",
            "mining_protein_id", "mining_component",
        ]
        if "mining_domain" in frame.columns:
            group_columns.insert(1, "mining_domain")
        groups = frame.groupby(group_columns, sort=False)
        progress.set(total_groups=groups.ngroups)
        for row_number, (_, group) in enumerate(groups, start=1):
            progress.set(current_group=row_number)
            ranked = group.copy()
            ranked["_component_count"] = pd.to_numeric(
                ranked.get("query_component_count", ""), errors="coerce"
            ).fillna(0)
            ranked = ranked.sort_values(
                ["_component_count", "query_genome_id", "query_contig_id",
                 "query_contig_system_index", "query_protein_id"],
                ascending=[False, True, True, True, True], kind="stable",
            )
            row = ranked.iloc[0]
            values = {column: str(row.get(column, "")).strip() for column in required}
            mining_domain = str(row.get("mining_domain", "")).strip()
            system = values["defense_system"]
            table = tables.get(system)
            key = (
                values["query_genome_id"], values["query_contig_id"],
                values["query_contig_system_index"],
            )
            reason = ""
            instance = pd.DataFrame()
            if table is None:
                reason = "query_system_resource_missing"
            else:
                instance = get_instance(table, *key)
                if instance.empty:
                    reason = "query_instance_missing"

            query_rows = instance[
                instance.get("protein_id", pd.Series(dtype=str)).astype(str).eq(values["query_protein_id"])
            ] if not instance.empty else instance
            if not reason and len(query_rows) != 1:
                reason = f"query_id_count_in_instance:{len(query_rows)}"
            if not reason and str(query_rows.iloc[0]["component"]).strip() != values["mining_component"]:
                reason = "query_component_context_mismatch"

            faa = contig_faa(args.prodigal, values["mining_genome_id"], values["mining_contig_id"])
            if not reason and (not faa.is_file() or faa.stat().st_size == 0):
                reason = "candidate_contig_faa_missing"

            own_component = values["mining_component"]
            expected = set(instance["component"].astype(str)) - {own_component} if not instance.empty else set()
            cache_key = (str(faa), system, *key, own_component)
            if not reason:
                if cache_key not in search_cache:
                    query_faa = temp / f"query_{len(search_cache):06d}.faa"
                    result_tsv = temp / f"hits_{len(search_cache):06d}.tsv"
                    query_map = write_partner_queries(instance, own_component, query_faa)
                    represented = {component for component, _ in query_map.values()}
                    absent_sequences = expected - represented
                    if absent_sequences:
                        search_cache[cache_key] = (None, "query_component_sequence_missing:" + "|".join(sorted(absent_sequences)))
                    elif expected:
                        database = ensure_database(args.diamond_bin, faa, database_cache)
                        search_count += 1
                        progress.set(searches=search_count)
                        run_search(args, database, query_faa, result_tsv)
                        search_cache[cache_key] = (read_hits(result_tsv, query_map), "")
                    else:
                        search_cache[cache_key] = (pd.DataFrame(), "")
                else:
                    cache_reused += 1
                hits, cached_reason = search_cache[cache_key]
                reason = cached_reason

            assignment = None
            usable_hits = pd.DataFrame()
            edges = {}
            if not reason:
                edges, usable_hits = best_edges(hits, values["mining_protein_id"])
                assignment = one_to_one_assignment(expected, edges)
                components_for_rule = set(expected) | {own_component}
                if assignment is None or not check_completeness(components_for_rule, system, rules):
                    reason = failure_reason(expected, edges, assignment)

            if reason:
                failures.append(group)
                progress.set(processed_groups=row_number, cache_reused=cache_reused,
                             assembled=len(assembly_records), unmatched=len(failures),
                             remaining_groups=groups.ngroups - row_number)
                continue

            target_components = selected_targets(usable_hits, assignment)
            instance_index = 300000 + len(assembly_records) + 1
            successes.append(make_mapping_rows(group, instance_index, "query_guided"))
            assembly_records.append(
                make_assembly_record(
                    system, mining_domain, values["mining_genome_id"], values["mining_contig_id"],
                    instance_index, values["mining_protein_id"], own_component, target_components,
                    "query_guided",
                )
            )

            progress.set(processed_groups=row_number, cache_reused=cache_reused,
                         assembled=len(assembly_records), unmatched=len(failures),
                         remaining_groups=groups.ngroups - row_number)
    progress.finish()
    success_frame = pd.concat(successes, ignore_index=True) if successes else pd.DataFrame(columns=success_columns)
    failure_frame = pd.concat(failures, ignore_index=True) if failures else pd.DataFrame(columns=failure_columns)
    write_csv_with_columns(
        success_frame, output_dir / "mapping.csv", success_columns
    )
    write_csv_with_columns(
        failure_frame, output_dir / "remain_unmatched.csv", failure_columns
    )
    assembly_frame = assembly_records_to_dataframe(assembly_records)
    if "mining_domain" in frame.columns and "mining_domain" not in assembly_frame.columns:
        assembly_frame.insert(1, "mining_domain", "")
    assembly_frame.to_csv(output_dir / "assembly.csv", index=False)
    print(
        f"query-guided assembly complete: {len(success_frame)} assembled, "
        f"{len(failure_frame)} remain unmatched"
    )


if __name__ == "__main__":
    main()
