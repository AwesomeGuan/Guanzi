# Step 1: Custom Mode

Custom Mode annotates user genomes with Prodigal, predicts one structure per
protein with ESMFold, builds a custom FoldSeek database and maps preset Query
structures to retained Mining proteins.

## Complete run

```bash
bash Step1_Custom_Mode/run_custom_mode.sh --genome_dir input/genomes --model_path esmfold_v1 --gpus 0,1
```

`--genome_dir`, `--model_path` and `--gpus` are required for a fresh run. The
paths passed between operations inside the chapter are fixed and are not
exposed as arguments.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_custom_mode.sh` | `--genome_dir PATH` | required | Directory containing MAG/genome FASTA files |
| | `--model_path PATH` | required | ESMFold model path |
| | `--gpus IDS` | required | Comma-separated GPU IDs (use `0` for GPU 0) |
| | `--query_pdb PATH` | `data/query_data/query_pdb` | Query PDB directory |
| | `--db_prefix NAME` | `custom` | FoldSeek database prefix |
| | `--num_threads N` | 12 | Parallel FoldSeek searches |
| | `--max_seq_len N` | 3500 | Maximum ESMFold sequence length |
| | `--thread_per_gpu N` | 1 | ESMFold model processes per GPU |
| | `--max_rounds N` | 5 | FoldSeek search rounds; `0` means unlimited |

## Stage commands

```bash
bash Step1_Custom_Mode/run_prodigal.sh
bash Step1_Custom_Mode/run_esmfold_prediction.sh
bash Step1_Custom_Mode/run_foldseek_db.sh --db_prefix custom
bash Step1_Custom_Mode/run_foldseek_search.sh --db_prefix custom
bash Step1_Custom_Mode/run_step1_custom.sh --input output/step1_raw/foldseek_search_results
```

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_prodigal.sh` | `--genome_dir PATH` | required | Directory containing `.fna`/`.fasta` genomes |
| `run_esmfold_prediction.sh` | `--model_path PATH` | required | ESMFold model path |
| | `--gpus IDS` | required | Comma-separated GPU IDs |
| | `--max_seq_len N` | 3500 | Maximum ESMFold sequence length |
| | `--thread_per_gpu N` | 1 | ESMFold model processes per GPU |
| | `--mode MODE` | `custom` | Log mode label |
| `run_foldseek_db.sh` | `--db_prefix NAME` | `custom` | Database name; writes `foldseek/{NAME}` and `foldseek/{NAME}_pad` |
| | `--threads N` | 16 | CPU threads |
| `run_foldseek_search.sh` | `--db_prefix NAME` | `custom` | FoldSeek database prefix |
| | `--gpus IDS` | required | Comma-separated GPU IDs |
| | `--query_pdb PATH` | `data/query_data/query_pdb` | Query PDB directory |
| | `--num_threads N` | 12 | Parallel FoldSeek processes |
| | `--max_rounds N` | 5 | Maximum search rounds; `0` means unlimited |
| `run_step1_custom.sh` | `--input PATH` | required | FoldSeek raw output from S01 (`output/step1_raw/foldseek_search_results`) |
| | `--output DIR` | `output/step1_mapping` | Final mapping output directory |
| | `--min_qcov FLOAT` | 0.70 | Minimum FoldSeek query coverage |
| | `--min_tcov FLOAT` | 0.70 | Minimum FoldSeek target coverage |
| | `--fident_blacklist_threshold FLOAT` | 0.30 | Target blacklist boundary |

Prodigal processes contigs sequentially and skips an existing `.faa`; it does
not expose a misleading thread option.

`--gpus` selects the GPUs used by both GPU stages:

- **ESMFold** distributes protein tasks round-robin across the selected GPUs.
  `--thread_per_gpu` is the number of ESMFold model processes loaded per GPU;
  the default of one avoids unnecessary model-memory duplication.
- **FoldSeek** `--num_threads` controls concurrent FoldSeek search subprocesses
  that share the selected GPUs. It is not a per-GPU worker count and not the
  CPU-thread count of one FoldSeek alignment.

## Key thresholds

| Parameter | Default | Rule |
|---|---|---|
| FoldSeek TM-score (built in) | 0.70 | FoldSeek keeps hits at or above the value |
| FoldSeek lDDT (built in) | 0.70 | FoldSeek keeps hits at or above the value |
| `--min_qcov` | 0.70 | Retain rows with FoldSeek query coverage at or above the value |
| `--min_tcov` | 0.70 | Retain rows with FoldSeek target coverage at or above the value |
| `--fident_blacklist_threshold` | 0.30 | Remove a target when any coverage-passing row has `fident` above the value; `fident == 0.30` is retained |

FoldSeek applies its built-in TM-score, lDDT and E-value filters first.
Post-processing then keeps the highest-TM-score row for each remaining target.

## Outputs

```text
output/prodigal_raw/<genome>/<contig>/<contig>.{faa,ffn,gff}
output/custom_pdb/<genome>/<contig>/<protein_id>.pdb
output/step1_raw/
├── esmfold_status.csv
├── foldseek_search_results/<system>/<query_id>
├── coverage_filter/<system>.csv
└── fident_filter/<system>.csv
output/step1_mapping/<system>.csv
```

FoldSeek intermediate CSV columns:

```text
query_id,subject_id,alntmscore,fident,qcov,tcov,evalue,lddt
```

Final Mapping columns:

```text
defense_system,query_protein_id,mining_protein_id,
foldseek_alntmscore,foldseek_sequence_identity,foldseek_query_coverage,
foldseek_target_coverage,foldseek_evalue,foldseek_lddt
```

Empty markers are stored in `output/.done/custom/step1/`. Existing outputs are
validated before a marker is reused. Delete the corresponding marker before an
intentional input, model, threshold or logic change.

## Empty candidate sets

The chapter stops normally when successful Prodigal processing produces no
protein records, ESMFold has no usable PDBs, FoldSeek has no valid hits, or the
coverage/fident filters retain no candidates. Counts exclude CSV headers and
preserve the existing thresholds. Later operations are not marked complete.
