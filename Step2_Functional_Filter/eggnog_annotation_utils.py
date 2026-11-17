"""Strict, dependency-free reading of native and normalized eggNOG annotations."""

import csv
import hashlib
import re
from pathlib import Path

REQUIRED_COLUMNS = {"Description", "PFAMs", "COG_category"}
COG_LETTERS = set("ABCDEFGHIJKLMNOPQRSTUVWYZ")


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_identifiers(path):
    identifiers = set()
    current = None
    has_sequence = False
    with Path(path).open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None and not has_sequence:
                    raise ValueError(f"Empty FASTA sequence: {current}")
                fields = line[1:].split()
                if not fields or fields[0] in identifiers:
                    raise ValueError(f"Empty or duplicate FASTA identifier at {path}:{number}")
                current = fields[0]
                identifiers.add(current)
                has_sequence = False
            else:
                if current is None:
                    raise ValueError(f"FASTA sequence precedes its header: {path}:{number}")
                has_sequence = True
    if current is not None and not has_sequence:
        raise ValueError(f"Empty FASTA sequence: {current}")
    return identifiers


def read_annotation_rows(path, allowed_ids=None, require_complete=False):
    """Keep all annotation fields, ignoring metadata but never malformed records."""
    header = None
    rows = []
    seen = set()
    complete = False
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            first = line.split("\t", 1)[0]
            if first in ("#query", "query"):
                complete = False
                candidate = next(csv.reader([line], delimiter="\t", quoting=csv.QUOTE_NONE))
                if header is not None:
                    if candidate != header:
                        raise ValueError(f"Conflicting eggNOG headers: {path}:{number}")
                    continue  # Native --resume can repeat its metadata/header.
                header = candidate
                if len(set(header)) != len(header):
                    raise ValueError(f"Duplicate eggNOG columns: {path}")
                missing = REQUIRED_COLUMNS - set(header)
                if missing:
                    raise ValueError(f"eggNOG table lacks required columns {sorted(missing)}: {path}")
                continue
            if line.startswith("#"):
                if line.startswith("## emapper-"):
                    complete = False
                if line.startswith("## Total time (seconds):"):
                    complete = True
                continue
            if header is None:
                raise ValueError(f"eggNOG record precedes #query/query header: {path}:{number}")
            complete = False  # A previous resume footer does not complete new records.
            row = next(csv.reader([line], delimiter="\t", quoting=csv.QUOTE_NONE))
            if len(row) != len(header):
                raise ValueError(f"Malformed eggNOG row: {path}:{number}; expected {len(header)} fields")
            query = row[0].strip()
            if not query or query in seen:
                raise ValueError(f"Empty or duplicate eggNOG query: {path}:{number}: {query}")
            if allowed_ids is not None and query not in allowed_ids:
                raise ValueError(f"eggNOG query is absent from the input FASTA: {query}")
            category = row[header.index("COG_category")].strip()
            if category not in ("", "-"):
                letters = re.sub(r"[,;\s]", "", category)
                if not letters or set(letters) - COG_LETTERS:
                    raise ValueError(f"Invalid COG functional category {category!r}: {path}:{number}")
            seen.add(query)
            rows.append(row)
    if header is None:
        raise ValueError(f"eggNOG table is missing #query/query header: {path}")
    if require_complete and not complete:
        raise ValueError(f"Native annotation has no final completion footer: {path}; wait for emapper to finish")
    return header, rows


def publish_annotation(path, header, rows):
    """Publish atomically, preserving the timestamp of an identical existing TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
            writer.writerow(header)
            writer.writerows(rows)
        if path.is_file() and file_digest(path) == file_digest(temporary):
            temporary.unlink()
            return False
        temporary.replace(path)
        return True
    finally:
        temporary.unlink(missing_ok=True)
