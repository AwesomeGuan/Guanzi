#!/usr/bin/env python3
"""Predict one PDB for every protein in recursive Prodigal FAA outputs."""

import argparse
import csv
import multiprocessing as mp
import os
import re
import sys
from pathlib import Path

MAX_SEQUENCE_LENGTH = 3500
TOKENIZER = None
MODEL = None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein_faa", required=True)
    parser.add_argument("--db_pdb", required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--gpus", required=True)
    parser.add_argument("--max_seq_len", type=int, default=3500)
    parser.add_argument("--thread_per_gpu", type=int, default=1)
    parser.add_argument("--status")
    parser.add_argument(
        "--include_parent",
        help="Process only FAA files whose immediate parent has this name",
    )
    return parser.parse_args()


def init_worker(gpu_id, model_path, max_sequence_length):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    print(f"[ESMFold/GPU{gpu_id}] INFO loading model", flush=True)

    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    global MAX_SEQUENCE_LENGTH, TOKENIZER, MODEL
    MAX_SEQUENCE_LENGTH = max_sequence_length
    TOKENIZER = AutoTokenizer.from_pretrained(model_path)
    MODEL = EsmForProteinFolding.from_pretrained(model_path, low_cpu_mem_usage=True)
    if not torch.cuda.is_available():
        raise RuntimeError(f"GPU {gpu_id} is unavailable")
    MODEL = MODEL.cuda()
    MODEL.esm = MODEL.esm.half()
    MODEL.trunk.set_chunk_size(32)
    torch.backends.cuda.matmul.allow_tf32 = True
    print(f"[ESMFold/GPU{gpu_id}] INFO model ready", flush=True)


def outputs_to_pdb(outputs):
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
    from transformers.models.esm.openfold_utils.protein import Protein, to_pdb

    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {key: value.to("cpu").numpy() for key, value in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for index in range(outputs["aatype"].shape[0]):
        protein = Protein(
            aatype=outputs["aatype"][index],
            atom_positions=final_atom_positions[index],
            atom_mask=final_atom_mask[index],
            residue_index=outputs["residue_index"][index] + 1,
            b_factors=outputs["plddt"][index],
            chain_index=(outputs["chain_index"][index] if "chain_index" in outputs else None),
        )
        pdbs.append(to_pdb(protein))
    return pdbs


def choose_chunk_size(sequence_length):
    if sequence_length > 2000:
        return 8
    if sequence_length > 1000:
        return 16
    return 32


def process_protein(task):
    faa_path, protein_id, sequence, pdb_path = task
    sequence_length = len(sequence)
    result = {
        "protein_id": protein_id,
        "faa_path": faa_path,
        "pdb_path": pdb_path,
        "sequence_length": sequence_length,
        "message": "",
    }
    if os.path.isfile(pdb_path) and os.path.getsize(pdb_path) > 0:
        return {**result, "status": "skipped"}
    if sequence_length > MAX_SEQUENCE_LENGTH:
        return {
            **result,
            "status": "skipped_too_long",
            "message": f"length exceeds {MAX_SEQUENCE_LENGTH}",
        }

    sequence = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYX]", "X", sequence.upper())
    if not sequence:
        return {**result, "status": "invalid_sequence", "message": "empty sequence"}

    import torch

    try:
        MODEL.trunk.set_chunk_size(choose_chunk_size(sequence_length))
        tokenized = TOKENIZER(
            [sequence], return_tensors="pt", add_special_tokens=False, padding=True
        )
        input_ids = tokenized["input_ids"].cuda()
        with torch.no_grad():
            outputs = MODEL(input_ids)
        pdb_text = "".join(outputs_to_pdb(outputs))
        Path(pdb_path).parent.mkdir(parents=True, exist_ok=True)
        temporary = f"{pdb_path}.partial"
        with open(temporary, "w", encoding="utf-8") as handle:
            handle.write(pdb_text)
        os.replace(temporary, pdb_path)
        del input_ids, outputs, pdb_text
        torch.cuda.empty_cache()
        return {**result, "status": "success"}
    except RuntimeError as exc:
        torch.cuda.empty_cache()
        message = str(exc)
        status = "out_of_memory" if "out of memory" in message.lower() else "runtime_error"
        return {**result, "status": status, "message": message[-1000:]}
    except Exception as exc:
        return {**result, "status": "error", "message": str(exc)[-1000:]}


def read_fasta(path):
    protein_id = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if protein_id is not None:
                    yield protein_id, "".join(sequence)
                protein_id = line[1:].split()[0]
                sequence = []
            elif protein_id is not None:
                sequence.append(line)
    if protein_id is not None:
        yield protein_id, "".join(sequence)


def collect_tasks(input_root, output_root, include_parent=None):
    tasks = []
    for faa_path in sorted(input_root.rglob("*.faa")):
        if include_parent and faa_path.parent.name != include_parent:
            continue
        output_dir = output_root / faa_path.parent.relative_to(input_root)
        for protein_id, sequence in read_fasta(faa_path):
            tasks.append(
                (
                    str(faa_path),
                    protein_id,
                    sequence,
                    str(output_dir / f"{protein_id}.pdb"),
                )
            )
    return tasks


def write_status(path, rows):
    columns = [
        "protein_id", "faa_path", "pdb_path", "sequence_length", "status", "message"
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    input_root = Path(args.protein_faa).resolve()
    output_root = Path(args.db_pdb).resolve()
    if not input_root.is_dir():
        raise SystemExit(f"Protein FAA directory not found: {input_root}")
    gpu_ids = [item.strip() for item in args.gpus.split(",") if item.strip()]
    if not gpu_ids:
        raise SystemExit("No GPU IDs were provided")
    if args.thread_per_gpu < 1:
        raise SystemExit("--thread_per_gpu must be at least 1")

    output_root.mkdir(parents=True, exist_ok=True)
    tasks = collect_tasks(input_root, output_root, args.include_parent)
    if not tasks:
        raise SystemExit(f"No protein records found under: {input_root}")
    expected_pdbs = {Path(task[3]).resolve() for task in tasks}
    for pdb_path in output_root.rglob("*.pdb"):
        if pdb_path.resolve() not in expected_pdbs:
            pdb_path.unlink()

    results = []
    pending = []
    for faa_path, protein_id, sequence, pdb_path in tasks:
        base = {
            "protein_id": protein_id,
            "faa_path": faa_path,
            "pdb_path": pdb_path,
            "sequence_length": len(sequence),
            "message": "",
        }
        if os.path.isfile(pdb_path) and os.path.getsize(pdb_path) > 0:
            results.append({**base, "status": "skipped"})
        elif len(sequence) > args.max_seq_len:
            results.append({
                **base,
                "status": "skipped_too_long",
                "message": f"length exceeds {args.max_seq_len}",
            })
        elif not sequence:
            results.append({**base, "status": "invalid_sequence", "message": "empty sequence"})
        else:
            pending.append((faa_path, protein_id, sequence, pdb_path))

    print(
        f"Proteins: {len(tasks)}, existing/terminal: {len(results)}, "
        f"pending: {len(pending)}, GPUs: {','.join(gpu_ids)}"
    )
    print(
        "[ESMFold] INFO "
        f"reused={sum(row['status'] == 'skipped' for row in results)} "
        f"too_long={sum(row['status'] == 'skipped_too_long' for row in results)} "
        f"invalid_sequence={sum(row['status'] == 'invalid_sequence' for row in results)} "
        f"processes_per_gpu={args.thread_per_gpu}",
        flush=True,
    )
    task_groups = [pending[index::len(gpu_ids)] for index in range(len(gpu_ids))]
    context = mp.get_context("spawn")
    pools = []
    asynchronous = []
    for gpu_id, gpu_tasks in zip(gpu_ids, task_groups):
        if not gpu_tasks:
            continue
        pool = context.Pool(
            processes=args.thread_per_gpu,
            initializer=init_worker,
            initargs=(gpu_id, args.model_path, args.max_seq_len),
        )
        pools.append(pool)
        asynchronous.append(pool.map_async(process_protein, gpu_tasks))

    for pool in pools:
        pool.close()
    try:
        for result in asynchronous:
            results.extend(result.get())
    finally:
        for pool in pools:
            pool.join()

    results.sort(key=lambda row: (row["faa_path"], row["protein_id"]))
    status_path = Path(args.status) if args.status else (
        Path(__file__).resolve().parent.parent
        / "output" / "step1_raw" / "esmfold_status.csv"
    )
    write_status(status_path, results)
    counts = {}
    for row in results:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print("ESMFold status: " + ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    print(f"[ESMFold] INFO status={status_path}", flush=True)

    retryable = sum(counts.get(key, 0) for key in ("out_of_memory", "runtime_error", "error"))
    if retryable:
        print(f"ESMFold has {retryable} retryable failure(s); rerun to retry them", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
