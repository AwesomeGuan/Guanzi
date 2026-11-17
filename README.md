# Guanzi

Guanzi mines structurally similar proteins that may replace
components of known prokaryotic defense systems. The Linux pipeline combines
Prodigal, ESMFold, FoldSeek, eggNOG, DefenseFinder, DIAMOND, CD-HIT and TM-align.

## Installation

Create the six conda environments, then deploy the resources required by the
selected mode:

```bash
bash setup_environments.sh
```

See [setup_environments.md](setup_environments.md) for software details and
[setup_resources.md](setup_resources.md) for models and reference data. Pipeline
configuration is passed through command-line arguments; no project-wide YAML or
environment-variable configuration is required.

This repository contains the Step 1–5 workflow and the fixed input tables.
Models, reference databases, Query structures and precomputed GTDB
DefenseFinder results are installed separately. See the resource setup guide
for which resources each mode requires.

## Pipeline modes

Custom Mode starts from user genomes, predicts their proteins and structures,
and searches the resulting custom FoldSeek database:

```bash
bash Step1_Custom_Mode/run_custom_mode.sh \
  --genome_dir input/genomes --model_path esmfold_v1 --gpus 0,1
```

Replication Mode searches AFDB50 and anchors retained proteins to bacterial and
archaeal GTDB r220 representatives:

```bash
bash Step1_Replication_Mode/run_prepare_gtdb_diamond_db.sh
bash Step1_Replication_Mode/run_replication_mode.sh --gpus 0,1
```

Step 2 performs local eggNOG annotation and selects the DefenseFinder source by
mode. Follow [Step2 Functional Filter](Step2_Functional_Filter/README.md) for
individual operations. The shared chapter commands are:

```bash
bash Step2_Functional_Filter/run_step2.sh --mode custom
bash Step3_Assembly_Analysis/run_step3.sh --mode custom
bash Step4_Quality_Control/run_step4.sh --mode custom
bash Step5_Information_Integration/run_step5.sh --mode custom
```

For Replication Mode, replace `custom` with `replication`. Step 4 predicts
structures for the retained GTDB Mining components before TM-align:

```bash
bash Step4_Quality_Control/run_step4.sh --mode replication \
  --model_path esmfold_v1 --gpus 0,1
```

There is no whole-project runner. Database deployment and the split
CPU/network/GPU server workflow remain explicit boundaries; Step2 runs locally
on a CPU host without a required web-service handoff.

## Chapters

1. [Step1 Custom Mode](Step1_Custom_Mode/README.md) or
   [Step1 Replication Mode](Step1_Replication_Mode/README.md): structural search
   and Query--Mining protein mapping.
2. [Functional Filter](Step2_Functional_Filter/README.md): eggNOG and
   DefenseFinder novelty filtering and system-style routing.
3. [Assembly Analysis](Step3_Assembly_Analysis/README.md): independent,
   complemented, expanded and Query-guided candidate assemblies.
4. [Quality Control](Step4_Quality_Control/README.md): distance, CD-HIT and
   TM-align filters, plus Replication homologous-frequency filtering.
5. [Information Integration](Step5_Information_Integration/README.md): component
   coordinates, transcription order and syntax for single co-transcribed instances,
   plus the final tables.

## Resume state

Empty stage markers are stored under `output/.done/<mode>/stepN/`. A marker is
reused only when the concrete stage output passes validation. Delete the
relevant marker when changing input data, thresholds, models or processing
logic. `output/.active_mode` prevents shared Step 2-5 outputs from being reused
across modes.

The AFDB50 PDB download is intentionally incremental. Its status CSV accounts
for every `(defense_system, protein_id)` as `success` or `failed`; existing PDBs
are skipped on retry. Newly downloaded structures propagate through missing FAA,
Query-exclusion and GTDB tasks without clearing every downstream marker.

After updating an existing server copy that used the former Step 1 marker
locations, migrate those markers once:

```bash
bash scripts/migrate_done_markers.sh
```

## GPU concurrency

GPU IDs are explicit for GPU stages: use `--gpus 0` on a single-GPU host or a
comma-separated list such as `--gpus 0,1`. The distinction between ESMFold
workers and FoldSeek concurrent searches is documented in
[Step1 Custom Mode](Step1_Custom_Mode/README.md) and
[Step1 Replication Mode](Step1_Replication_Mode/README.md).

## Main outputs

Step 2 writes protein-level novelty summaries and routes system mappings under
`output/step2_output/`. Step 4 and Step 5 use the canonical Mapping, Assembly and
Comparison schemas. The final deliverables are:

```text
output/step5_output/
├── mapping_qualified.csv
├── assembly_qualified.csv
├── comparison_qualified.csv
├── mapping_disqualified.csv
├── assembly_disqualified.csv
└── comparison_disqualified.csv
```

Mapping rows connect Query and Mining proteins. Assembly rows describe Mining
system instances. Comparison rows hold component-level Query--Mining structural
evidence.

## Runtime messages

Each operation reports `START` and `END` with its mode, CPU or selected GPU
IDs, and elapsed time measured by a monotonic clock. Shell-owned messages do
not include wall-clock timestamps. Original tool output remains visible.
CPU thread options are optional and retain their existing defaults; for
example, GTDB mapping uses 12 threads unless `--threads N` is supplied.

Long operations report an in-memory `PROGRESS` snapshot every 300 seconds (5 minutes).
Operations without reliable item counts report `RUNNING` instead; ESMFold
retains its existing GPU batches and does not claim per-protein live counts.
Fast operations need no intermediate message. `SUMMARY` preserves
operation-specific counts, including reuse, new results, exclusions and
failures. `SKIP` means validated output reuse with no processing; its elapsed
value is zero processing time, not the preceding validation duration.
Nonzero exits retain their original code and report `ERROR`. Local eggNOG
annotation reports unannotated proteins separately from execution failures.

Individual commands are documented in each chapter README. In the shared
Step 2-5 chapters, `--mode` selects the mode-specific default inputs as well as
the log label. Step 1 modes have fixed mode labels.
Logging does not modify completion markers or processing decisions.

## Empty candidate sets

Chapter runners count operation-specific outputs after successful execution or
validated reuse. A valid zero candidate set reports `SUMMARY` and `STOP` and
ends the chapter with exit code 0. It does not mark unexecuted operations as
complete or delete downstream outputs. Missing or malformed required files and
execution failures remain errors. Do not continue to the next chapter after a
zero-input stop unless it has another valid input source.

CSV counts exclude headers; FASTA counts use sequence records; AlphaFold success
counts use the current status table and existing successful PDBs. Step 3 treats
remaining candidates separately from already assembled instances. Step 4 and
Step 5 consider both qualified and disqualified instances.
