# Guanzi resource setup

Run all commands from the repository root on Linux after creating the conda
environments. This document installs models and data only. Server-specific
storage layouts are outside the reproducible project procedure; every resource
must appear at the canonical path shown here when the pipeline runs.

## 1. Resources distributed with the repository

The repository includes the normalized Query instance tables and the small,
fixed rule tables:

```text
data/query_data/query_location/*.csv
data/DefenseFinder_rules/DefenseFinder_rules.tsv
data/DefenseFinder_rules/DefenseFinder_rules_special.tsv
data/defense_system_style.csv
data/eggNOG_Description_out_function.csv
data/eggNOG_PFAMs_out_function.csv
```

`query_location` is the runtime dataset. No pre-normalization source directory
is required by the published pipeline.

The curated Query PDB collection must be installed as one directory per defense
system:

```text
data/query_data/query_pdb/<defense_system>/<query_protein_id>.pdb
```

Copy the Query PDB archive supplied with the project release into this path.
The release URL should be recorded here when the public data archive is issued;
the pipeline does not manufacture missing Query structures.

Download the DefenseFinder model definitions after installing its environment:

```bash
bash scripts/run_in_env.sh defensefinder defense-finder update
```

## 2. ESMFold model and Custom Mode input

Place input genome or MAG assemblies under `input/genomes/`, or pass another
directory explicitly to the Step 1 wrapper. Accepted files are FASTA-family
sequence files documented in `Step1_Custom_Mode/README.md`.

Download ESMFold v1 through ModelScope and place the model under the repository
root. Custom Mode uses it in Step 1; Replication Mode uses it for retained GTDB
Mining components in Step 4:

```bash
bash scripts/run_in_env.sh esmfold modelscope download --model facebook/esmfold_v1
mv ./facebook/esmfold_v1 ./esmfold_v1
```

Expected model path:

```text
esmfold_v1/
```

Custom Mode creates its FoldSeek database under `foldseek/` during Step 1.

## 3. Replication Mode: FoldSeek AFDB50

Download and extract the fixed AFDB50 archive, then create the index and GPU
padded database used by the structural search:

```bash
mkdir -p foldseek
wget -c https://foldseek-cdn.mmseqs.org/afdb50.tar.gz \
  -O foldseek/afdb50.tar.gz
tar -xzf foldseek/afdb50.tar.gz -C foldseek
bash scripts/run_in_env.sh foldseek foldseek createindex \
  foldseek/afdb50 foldseek/tmp
bash scripts/run_in_env.sh foldseek foldseek makepaddedseqdb \
  foldseek/afdb50 foldseek/afdb50_pad
```

Replication Step 1 reads the prefix `foldseek/afdb50_pad`.

## 4. Replication Mode: GTDB r220 genomes

```bash
mkdir -p data/gtdb/downloads
wget -c \
  https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_genomes_reps_r220.tar.gz \
  -O data/gtdb/downloads/gtdb_genomes_reps_r220.tar.gz
tar -xzf data/gtdb/downloads/gtdb_genomes_reps_r220.tar.gz \
  -C data/gtdb
```

Confirm that the canonical root contains the genome archives and both lookup
tables:

```text
data/gtdb/gtdb_genomes_reps_r220/
├── database/
├── genome_paths.tsv
└── file_paths.tsv
```

## 5. Replication Mode: GTDB r220 representative proteins

```bash
mkdir -p data/gtdb/downloads
wget -c \
  https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_proteins_aa_reps_r220.tar.gz \
  -O data/gtdb/downloads/gtdb_proteins_aa_reps_r220.tar.gz
tar -xzf data/gtdb/downloads/gtdb_proteins_aa_reps_r220.tar.gz \
  -C data/gtdb

bash scripts/unpack_gtdb_proteins.sh
```

The tar archive contains individually compressed `*_protein.faa.gz` files.
The unpacking script decompresses both domains in place and removes each gzip
file only after successful publication. Existing FAA files are not overwritten;
when both forms exist, their decoded contents must match before the gzip file
is removed. Rerunning processes the remaining compressed files. No conda
environment or completion marker is required.

Confirm that the archive's top-level `protein_faa_reps/` directory now has this
layout without a separate move, preserving every `*_protein.faa` filename:

```text
data/gtdb/protein_faa_reps/
├── bacteria/*_protein.faa
└── archaea/*_protein.faa
```

Build the reusable DIAMOND databases and location indexes:

```bash
bash Step1_Replication_Mode/run_prepare_gtdb_diamond_db.sh --threads 12
```

This creates `data/gtdb/diamond_db/`. Replication Mode also creates reusable
per-system Query databases under `data/query_data/query_diamond_db/`.

## 6. Replication Mode: precomputed GTDB DefenseFinder results

Replication Step 2 imports DefenseFinder output produced for the same GTDB r220
representative proteins. Install the complete, version-matched result tree at:

```text
data/gtdb/defensefinder_results/
├── bacteria/<genome_id>_protein/*_defense_finder_systems.tsv
└── archaea/<genome_id>_protein/*_defense_finder_systems.tsv
```

This project does not currently publish or generate the full precomputed tree.
Do not substitute output from another GTDB release. The public data archive
location must be added here when it is issued.

## 7. Replication Mode: Mining structure prediction

No external GTDB Mining PDB collection is required. After assembly filtering
and CD-HIT clustering, Replication Step 4 extracts the retained GTDB Mining
component sequences and predicts their structures with the ESMFold model
installed in Section 2. Structures are written under
`output/step4_raw/mining_pdb/`.

## 8. eggNOG 5.0.2 local annotation databases

After creating the `eggnog` environment, deploy the DIAMOND-route resources:

```bash
bash scripts/download_eggnog_data.sh
```

The script downloads from the official eggNOG 5 host at
<http://eggnog5.embl.de/download/emapperdb-5.0.2/>. It resumes incomplete
downloads, checks compressed archives, validates extracted resources and reuses
existing valid databases without overwriting them. Invalid existing resources
cause a clear error and require manual repair. No optional HMMER, MMseqs2,
novel-family or de novo Pfam databases are downloaded.

```text
data/eggnog_5.0.2/
├── eggnog.db                 # annotation database
├── eggnog.taxa.db            # taxonomy database
├── eggnog_proteins.dmnd      # DIAMOND protein database
└── ...                      # companion taxonomy files, if supplied
```

The compressed downloads total approximately 11.2 GiB. Extracted resources need
roughly 50 GB; allow additional space for archives, temporary extraction and
annotation results. Archives are retained and may be removed after validation.
Existing manually deployed databases at the canonical path can be reused.

Step2 uses eggNOG-mapper 2.1.13 with database 5.0.2. Native annotation output is
preserved and normalized to `output/step2_raw/eggnog_results.tsv` only after
successful validation. The required annotation fields are `Description`,
`PFAMs` and `COG_category`; v3 output missing these fields is not a substitute.
See `Step2_Functional_Filter/README.md` for execution.

## 9. Resource checks

Before starting a mode, inspect its chapter README and run the corresponding
wrapper. Each wrapper checks the concrete resources it consumes. Large download
archives may be removed after extraction and validation.
