# Guanzi

Guanzi uses protein structure searches to find potential replacements for
components of known prokaryotic defense systems. It combines structural matches,
gene neighborhoods, and system rules to reconstruct candidate multi-gene
defense systems and report the evidence for each component.

## Installation

On a Linux system with Conda, run this command from the repository root:

```bash
bash setup_environments.sh
```

See the [environment guide](setup_environments.md) for software requirements and
the [resource guide](setup_resources.md) for the models and reference data each
mode needs. Deploy the resources for your chosen mode before running Step 1.
This repository includes the workflow and input tables. Large resources must be
installed separately; the resource guide lists their sources and archive
availability.

## Quick start

Choose one Step 1 entry and use its mode for Steps 2–5. Run the commands in
order from the repository root. For GPU stages, use `--gpus 0` for one GPU or
`--gpus 0,1` for two GPUs, adjusting the IDs to match your host.

### Custom Mode: analyze your own genomes

Place your genome FASTA files in `input/genomes/` or pass another directory with
`--genome_dir`. Guanzi predicts protein structures from those genomes and
searches a Foldseek database built from them.

```bash
bash Step1_Custom_Mode/run_custom_mode.sh \
  --genome_dir input/genomes --model_path esmfold_v1 --gpus 0,1
bash Step2_Functional_Filter/run_step2.sh --mode custom
bash Step3_Assembly_Analysis/run_step3.sh --mode custom
bash Step4_Quality_Control/run_step4.sh --mode custom
bash Step5_Information_Integration/run_step5.sh --mode custom
```

### Replication Mode: run the GTDB r220 reference workflow

Guanzi searches AFDB50 and maps retained hits to GTDB r220 representative
genomes. Step 4 also predicts structures for retained candidate components, so
a fresh run needs `--model_path` and `--gpus` at that step.

```bash
bash Step1_Replication_Mode/run_prepare_gtdb_diamond_db.sh
bash Step1_Replication_Mode/run_replication_mode.sh --gpus 0,1
bash Step2_Functional_Filter/run_step2.sh --mode replication
bash Step3_Assembly_Analysis/run_step3.sh --mode replication
bash Step4_Quality_Control/run_step4.sh --mode replication \
  --model_path esmfold_v1 --gpus 0,1
bash Step5_Information_Integration/run_step5.sh --mode replication
```

If a chapter reports `STOP` because no candidates remain, end that run there.
See the chapter guides below for individual commands and how to resume a run.

## Workflow layout

```text
Guanzi/
├── Step1_Custom_Mode/              # Search protein structures from user genomes
├── Step1_Replication_Mode/         # Search AFDB50 and map hits to GTDB r220
├── Step2_Functional_Filter/        # Annotate and filter protein candidates
├── Step3_Assembly_Analysis/        # Reconstruct candidate defense systems
├── Step4_Quality_Control/          # Check gene proximity and structural similarity
├── Step5_Information_Integration/  # Resolve component coordinates and order
├── data/                           # Query metadata, system rules, and filter terms
├── input/genomes/                  # Default location for user genome FASTA files
├── setup_environments.md           # Software and environment setup
└── setup_resources.md              # Model and reference data setup
```

Chapter guides: [Step 1 Custom](Step1_Custom_Mode/README.md) ·
[Step 1 Replication](Step1_Replication_Mode/README.md) ·
[Step 2](Step2_Functional_Filter/README.md) ·
[Step 3](Step3_Assembly_Analysis/README.md) ·
[Step 4](Step4_Quality_Control/README.md) ·
[Step 5](Step5_Information_Integration/README.md).

## Outputs

Step 5 reports candidate systems with resolved component coordinates on the same
coding strand. The final tables retain their structural qualification labels
from Step 4.

```text
output/step5_output/
├── mapping_qualified.csv
├── assembly_qualified.csv
├── comparison_qualified.csv
├── mapping_disqualified.csv
├── assembly_disqualified.csv
└── comparison_disqualified.csv
```

**Mapping** links reference (Query) and candidate (Mining) proteins.
**Assembly** describes candidate system instances. **Comparison** reports
component-level structural evidence and the resolved component order. See the
[Step 5 output guide](Step5_Information_Integration/README.md#outputs) for the
table fields.
