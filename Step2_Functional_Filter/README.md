# Step 2: Functional Filter

Step2 prepares unique candidate proteins, annotates them locally with eggNOG,
excludes defense-related annotations, and classifies the retained mappings after
DefenseFinder filtering. Run all commands from the project root.

## Complete run

Custom Mode:

```bash
bash Step2_Functional_Filter/run_step2.sh --mode custom
```

Replication Mode:

```bash
bash Step2_Functional_Filter/run_step2.sh --mode replication
```

Wrappers automatically use `defense_system` for tables, `eggnog` for annotation,
and `defensefinder` for Custom annotation. Install environments with
`setup_environments.sh` and deploy databases through `setup_resources.md`.
Annotation requires no web-service pause or runtime database download.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_step2.sh` | `--mode MODE` | required | Custom or Replication |
| | `--threads N` | 12 | eggNOG/DIAMOND CPUs, not concurrent jobs |

## Stage commands

Use the same mode throughout. Operations are numbered S01-S08 in execution order.

### S01: Prepare candidate proteins

Custom Mode:

```bash
bash Step2_Functional_Filter/run_prepare_eggnog.sh --mode custom
```

Replication Mode:

```bash
bash Step2_Functional_Filter/run_prepare_eggnog.sh --mode replication
```

This operation prepares `output/step2_raw/eggnog_upload.faa` for local annotation.

### S02: Local eggNOG annotation

```bash
bash Step2_Functional_Filter/run_eggnog_annotation.sh
```

This operation is shared by both modes and uses the active mode selected by S01.
It runs eggNOG-mapper 2.1.13 with database 5.0.2 and DIAMOND 2.2.8, using:

```text
--itype proteins -m diamond
--sensmode sensitive --dmnd_iterate yes --dmnd_ignore_warnings
--matrix BLOSUM62 --gapopen 11 --gapextend 1
--tax_scope auto --evalue 0.001 --target_orthologs all
```

These are fixed annotation settings, not optional pipeline arguments. No extra
identity or coverage cutoff is imposed; PFAM terms are transferred from orthologs.
The native annotations are preserved and automatically converted to
`output/step2_raw/eggnog_results.tsv` for downstream filtering.

### S03: eggNOG functional exclusion

Custom Mode:

```bash
bash Step2_Functional_Filter/run_eggnog_filter.sh --mode custom
```

Replication Mode:

```bash
bash Step2_Functional_Filter/run_eggnog_filter.sh --mode replication
```

Existing rules are unchanged: system-specific Description/PFAM keywords and
COG category V. `Description`, `PFAMs` and `COG_category` must all exist.
An individual protein can be unannotated; missing columns or COG identifiers
instead of category letters are errors.

### S04: DefenseFinder annotation or import

Custom Mode:

```bash
bash Step2_Functional_Filter/run_defensefinder.sh --mode custom
```

Replication Mode:

```bash
bash Step2_Functional_Filter/run_replication_defensefinder_import.sh --mode replication
```

Custom Mode annotates retained contigs; Replication Mode imports the
release-matched GTDB tables.

### S05: DefenseFinder novelty filtering

Custom Mode:

```bash
bash Step2_Functional_Filter/run_defensefinder_filter.sh --mode custom
```

Replication Mode:

```bash
bash Step2_Functional_Filter/run_defensefinder_filter.sh --mode replication
```

### S06-S08: Classification and location context

Custom Mode:

```bash
bash Step2_Functional_Filter/run_classify_mining_systems.sh --mode custom
```

Replication Mode:

```bash
bash Step2_Functional_Filter/run_classify_mining_systems.sh --mode replication
```

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| Mode-specific stage wrappers | `--mode MODE` | required | Custom or Replication |
| `run_prepare_eggnog.sh` | `--input DIR` | mode-selected mapping | Step1 input |
| | `--protein_faa DIR` | mode-selected Prodigal | Protein input |
| | `--output DIR` | `output/step2_raw` | Preparation output |
| `run_eggnog_annotation.sh` | `--threads N` | 12 | eggNOG/DIAMOND CPUs |
| `run_eggnog_filter.sh` | `--eggnog_results FILE` | `output/step2_raw/eggnog_results.tsv` | Annotation input |
| | `--description_keywords FILE` | `data/eggNOG_Description_out_function.csv` | Description rules |
| | `--pfam_keywords FILE` | `data/eggNOG_PFAMs_out_function.csv` | PFAM rules |
| `run_defensefinder.sh` | `--models_dir DIR` | bundled | DefenseFinder models |
| `run_replication_defensefinder_import.sh` | `--source DIR` | `data/gtdb/defensefinder_results` | Reference results |
| `run_classify_mining_systems.sh` | `--style_file FILE` | `data/defense_system_style.csv` | System classification |
| | `--rules DIR` | `data/DefenseFinder_rules` | Rules |
| | `--query_location DIR` | `data/query_data/query_location` | Query metadata |

Other retained path overrides are documented by each wrapper's `--help`.
The chapter runner connects operations through their standard paths.

## Outputs

| Input | Custom | Replication |
|---|---|---|
| Step1 mapping | `output/step1_mapping/` | `output/rep_step1_mapping/` |
| Prodigal tree | `output/prodigal_raw/` | `output/rep_prodigal_raw/` |
| DefenseFinder | Run on retained contigs | Import GTDB r220 results |

```text
output/step2_raw/
├── eggnog_upload.faa
├── eggnog_extraction_failures.csv               # only when needed
├── eggnog_results/
│   ├── query.emapper.annotations                # preserved native table
│   ├── query.emapper.hits
│   └── query.emapper.seed_orthologs
├── eggnog_results.tsv                          # normalized downstream input
├── eggnog_filtered/<system>.csv
├── defensefinder_results/<genome>/<contig>/...
└── defensefinder_filtered/<system>.csv
output/
├── isozyme_candidates_eggnog_annotated.csv
├── isozyme_candidates_defensefinder_annotated.csv
├── isozyme_candidates_novel.csv
└── step2_output/
    ├── step2_mapping_single/<system>.csv
    ├── step2_mapping_multi/<system>.csv
    ├── step2_mapping_complex/<system>.csv
    ├── step2_mapping_cas/<system>.csv
    ├── step2_summary_single.csv
    ├── step2_summary_multi.csv
    ├── step2_summary_complex.csv
    └── step2_summary_cas.csv
```

## Resume and empty candidate sets

Completion markers reside under `output/.done/<mode>/step2/`. Completed operations
are skipped after checking their outputs. S02 reuses complete native annotations
and generates the standard TSV automatically; missing or invalid outputs are
not treated as success. Long-running operations report a heartbeat every five
minutes and completion messages include elapsed time.

Empty candidates stop normally with exit code 0. Only `step2_mapping_multi`
enters Step3; other classifications remain available as Step2 outputs.
