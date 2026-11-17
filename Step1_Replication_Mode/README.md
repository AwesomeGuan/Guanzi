# Step 1: Replication Mode

Replication Mode searches the GPU-padded AFDB50 database, applies the same
coverage and target-level fident rules as Custom Mode, excludes Query homologs,
and maps the retained AFDB50 proteins to GTDB r220 representatives.

## Complete run

```bash
bash Step1_Replication_Mode/run_prepare_gtdb_diamond_db.sh
bash Step1_Replication_Mode/run_replication_mode.sh --gpus 0,1
```

`--gpus` is required for a fresh run. The paths passed between operations
inside the chapter are fixed and are not exposed as arguments.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_prepare_gtdb_diamond_db.sh` | `--threads N` | 12 | DIAMOND threads for the GTDB reference build |
| | `--force` | off | Rebuild the GTDB DIAMOND database even when it exists |
| `run_replication_mode.sh` | `--gpus IDS` | required | Comma-separated GPU IDs (use `0` for GPU 0) |
| | `--query_pdb PATH` | `data/query_data/query_pdb` | Query PDB directory |
| | `--num_threads N` | 12 | Parallel FoldSeek processes |
| | `--download_threads N` | 8 | Concurrent AlphaFold downloads |
| | `--max_rounds N` | 5 | Maximum FoldSeek search rounds; `0` means unlimited |

## Stage commands

```bash
bash Step1_Replication_Mode/run_foldseek_search.sh
bash Step1_Replication_Mode/run_afdb50_filter.sh
bash Step1_Replication_Mode/run_afdb50_pdb_download.sh
bash Step1_Replication_Mode/run_query_exclusion.sh
bash Step1_Replication_Mode/run_gtdb_mapping.sh
bash Step1_Replication_Mode/run_gtdb_prodigal_materialization.sh
bash Step1_Replication_Mode/run_mapping_extraction.sh
```

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_foldseek_search.sh` | `--gpus IDS` | required | Comma-separated GPU IDs |
| | `--query_pdb PATH` | `data/query_data/query_pdb` | Query PDB directory |
| | `--num_threads N` | 12 | Parallel FoldSeek processes |
| | `--max_rounds N` | 5 | Maximum search rounds; `0` means unlimited |
| `run_afdb50_filter.sh` | `--min_qcov FLOAT` | 0.70 | Minimum FoldSeek query coverage |
| | `--min_tcov FLOAT` | 0.70 | Minimum FoldSeek target coverage |
| | `--fident_blacklist_threshold FLOAT` | 0.30 | Remove a target when any coverage-passing row has `fident` above the value; `fident == 0.30` is retained |
| `run_afdb50_pdb_download.sh` | `--threads N` | 8 | Concurrent AlphaFold downloads |
| `run_query_exclusion.sh` | `--query_location PATH` | `data/query_data/query_location` | Query-location reference for homolog exclusion |
| `run_gtdb_mapping.sh` | `--threads N` | 12 | DIAMOND threads |
| | `--min_identity FLOAT` | 90 | Minimum GTDB DIAMOND identity percent |
| | `--min_qcov FLOAT` | 90 | Minimum GTDB query coverage percent |
| | `--min_tcov FLOAT` | 90 | Minimum GTDB target coverage percent |
| | `--evalue FLOAT` | 1e-5 | Maximum GTDB DIAMOND E-value |
| | `--max_target_seqs N` | 0 | Maximum reporting targets; `0` means unlimited |
| `run_gtdb_prodigal_materialization.sh` | – | – | No optional arguments |
| `run_mapping_extraction.sh` | – | – | No optional arguments |

`run_gtdb_mapping.sh` builds the GTDB DIAMOND database through
`run_prepare_gtdb_diamond_db.sh --threads "$THREADS"` before mapping, so the
database build shares the same thread count and is not listed separately.

`--gpus` selects the GPUs used by the FoldSeek search stage:

- **FoldSeek** `--num_threads` controls concurrent FoldSeek search subprocesses
  that share the selected GPUs. It is not a per-GPU worker count and not the
  CPU-thread count of one FoldSeek alignment.
- **AlphaFold downloads** `--download_threads` controls concurrent AlphaFold
  database downloads and has no GPU meaning.

The PDB download pass is incremental and fully accounted for. Its status table
lists every input pair:

```text
defense_system,protein_id,status,message
```

Accounting rules:

- Every pair appears as `success` or `failed`, and the nonempty PDB set must
  equal the `success` set. If the two counts differ, the S05 and S06 states are
  reconciled.
- Valid expected PDB and FAA files are reused; only missing files are added, and
  stale or invalid files are removed.
- Removing only `05_download_isozyme_pdb.done` retries the missing files; newly
  successful files are processed incrementally by the later stages.
- GTDB mapping records handled Query proteins in `processed_queries.csv` and
  runs DIAMOND only for newly retained proteins.

## Key thresholds

| Parameter | Default | Rule |
|---|---|---|
| `--foldseek_min_qcov` | 0.70 | Retain FoldSeek rows with query coverage at or above the value |
| `--foldseek_min_tcov` | 0.70 | Retain FoldSeek rows with target coverage at or above the value |
| `--foldseek_fident_blacklist_threshold` | 0.30 | Remove a target when any coverage-passing row has `fident` above the value; `fident == 0.30` is retained |
| Query-exclusion DIAMOND `--id` | 70 | Exclude close Query homologs at or above 70% identity |
| Query-exclusion DIAMOND `--query-cover` | 70 | Require at least 70% candidate-sequence coverage |
| Query-exclusion DIAMOND `--subject-cover` | 70 | Require at least 70% Query-reference coverage |
| `--gtdb_min_identity` | 90 | Minimum GTDB DIAMOND identity percent |
| `--gtdb_min_qcov` | 90 | Minimum GTDB query coverage percent |
| `--gtdb_min_tcov` | 90 | Minimum GTDB target coverage percent |
| `--gtdb_evalue` | 1e-5 | Maximum GTDB DIAMOND E-value |
| `--gtdb_max_target_seqs` | 0 | Maximum reporting targets; `0` means unlimited |

## Outputs

```text
output/rep_step1_raw/
├── foldseek_search_results/<system>/<query_id>
├── coverage_filter/<system>.csv
├── fident_filter/<system>.csv
├── afdb50_isozyme_ids/<system>.txt
├── afdb50_isozyme_pdb/<system>/<afdb50_id>.pdb
├── afdb50_pdb_download_status.csv
├── afdb50_isozyme_faa/<system>/<afdb50_id>.faa
├── diamond_filtered_faa/<system>/<afdb50_id>.faa
├── gtdb_mapping/{bacteria,archaea}/<system>.csv
└── candidate_genomes/<genome>/<contig>.fasta
output/rep_prodigal_raw/<genome>/<contig>/<contig>.{faa,ffn,gff}
output/rep_step1_mapping/<system>.csv
```

GTDB mapping columns:

```text
defense_system,domain,genome_id,contig_id,afdb50_id,gtdb_protein_id,
gtdb_pident,gtdb_qcov,gtdb_tcov,gtdb_evalue,gtdb_bitscore
```

Final Mapping columns:

```text
defense_system,query_protein_id,afdb50_id,mining_domain,
mining_genome_id,mining_contig_id,mining_protein_id
```

No GTDB mapping class is generated. Empty markers are stored in
`output/.done/replication/step1/`.

## Empty candidate sets

The chapter ends normally with exit code 0 when there are no valid FoldSeek
hits, candidate IDs, successful PDBs, retained Query-excluded proteins, GTDB
mapping rows, or materialized proteins. The PDB/FAA wrapper also stops before
FAA extraction if the validated status table records zero successful downloads.
Failures remain accounted for; later operations are not marked complete. GTDB
mapping counts combine bacteria and archaea and exclude CSV headers.
