# Step 4: Quality Control

Step 4 applies genomic distance filtering, CD-HIT instance clustering and
component-level TM-align structural quality control. Replication Mode adds a
GTDB co-contig recurrence filter.

The inputs are initial candidate multi-gene defense systems (CMDSs) constructed
in Step 3. Retained outputs remain CMDSs, with qualified/disqualified labels
recorded separately. CMDS identifies the candidate object, not a quality-control
pass or experimentally established function. Assembly continues to name the
construction process and the unchanged table/file class.

## Complete run

Custom Mode:

```bash
bash Step4_Quality_Control/run_step4.sh --mode custom
```

Replication Mode:

```bash
bash Step4_Quality_Control/run_step4.sh --mode replication \
  --model_path esmfold_v1 --gpus 0,1
```

`--mode` selects the mode's completion markers and mode-specific inputs.
Replication Mode predicts the retained Mining structures with ESMFold, so
`--model_path` and `--gpus` are required for a fresh Replication run. Empty
markers are stored in `output/.done/<mode>/step4/`.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_step4.sh` | `--mode MODE` | required | Select the Custom or Replication input layout |
| | `--prodigal_dir DIR` | derived from `--mode` | Override the mode-specific Prodigal directory |
| | `--model_path PATH` | – | ESMFold model; required for fresh Replication structures |
| | `--gpus IDS` | – | Comma-separated GPU IDs; required for fresh Replication structures |
| | `--max_seq_len N` | 3500 | Maximum ESMFold sequence length |
| | `--thread_per_gpu N` | 1 | ESMFold model processes per GPU |
| | `--threads N` | 12 | CPU threads/workers for CD-HIT, TM-align and each DIAMOND worker |
| | `--workers N` | 4 | Concurrent DIAMOND workers for the Replication frequency filter |
| | `--shard_size N` | 128 | Unique queries per DIAMOND shard |
| | `--max_gene_gap N` | 10 | Maximum gene-index difference for a connectivity edge |

## Stage commands

The lower-level entry points support step-by-step or recovery runs. Distance
filtering, clustering and TM-align run on a CPU host; the Replication Mining
ESMFold prediction needs a GPU host.

### Distance filter (S01)

```bash
bash Step4_Quality_Control/run_distance_filter.sh
```

S01 treats component genes on the candidate contig as nodes. An edge joins
nodes whose absolute gene-index difference is at most `--max_gene_gap`
(default 10), and the operation retains the connected component containing
the Mining anchor. The numeric index is parsed from the protein-ID suffix;
it is not a nucleotide distance. Chained connections are allowed: indices
1, 10 and 19 can remain connected even though the endpoints differ by 18.

Within that connected component, the anchor is retained and redundant
rule-defined component roles are represented by the copy nearest the anchor,
with deterministic index/ID tie-breaks. The resulting set must satisfy the
mandatory, accessory and forbidden DefenseFinder rules. An unresolved anchor
index or a rule-incomplete set removes the instance.

### Clustering (S02–S04)

Custom Mode:

```bash
bash Step4_Quality_Control/run_cdhit_cluster.sh --mode custom
```

Replication Mode:

```bash
bash Step4_Quality_Control/run_cdhit_cluster.sh --mode replication
bash Step4_Quality_Control/run_mining_esmfold.sh \
  --model_path esmfold_v1 --gpus 0,1
```

The Replication branch then predicts structures for the Mining components
retained by clustering with ESMFold on a GPU host.

Stage 2 owns its working tree: before S02 runs it clears `qma_faa/`,
`qma_faa_clustered/`, `clustering_temp/` and the two `clustered_*.csv` tables,
so the result always reflects the current Step 3 output and can never reuse
instances or components that Step 3 no longer contains. The stage is skipped
when its completion marker and both `clustered_*.csv` tables are present, so the
cleanup is never applied to a completed run.

### TM-align (S05–S08)

Custom Mode:

```bash
bash Step4_Quality_Control/run_tmalign.sh --mode custom
```

Replication Mode:

```bash
bash Step4_Quality_Control/run_tmalign.sh --mode replication
```

Custom Mode compares against the Mining structures predicted across all proteins
in Step 1 (`output/custom_pdb`); Replication Mode compares against the structures
predicted only for the components retained by clustering
(`output/step4_raw/mining_pdb`). `--mode` selects between the two libraries.

Components are paired by exact name first, then one-to-one within the same
DefenseFinder rule-defined role, including wildcard and interchangeable-profile
groups. Remaining same-role alternatives follow the historical name-distance
ordering with deterministic name tie-breaks; names alone never permit
cross-role pairing. S05 gives the paired Query PDB the Mining filename expected
by S06. The source FAA names are unchanged, and S07 reports the actual Query
and Mining component names separately. Missing reference components or PDBs
are reported separately from TM-align execution failures.

### Co-contig recurrence (S09, Replication only)

```bash
bash Step4_Quality_Control/run_replication_homologous_frequency.sh
```

For each component, S09 collects the distinct GTDB contigs containing a
qualifying DIAMOND hit. The co-contig recurrence count is the size of the
intersection of these sets across all components. An instance passes when
at least five contigs contain a qualifying hit for every component.
This operation checks co-contig detection; it does not re-check distances,
strand directions or component order on those target contigs. The existing
script, directory and `homologous_assembly_count` column names are retained
for compatibility; the column records this contig-intersection count.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_distance_filter.sh` | `--input DIR` | `output/step3_output` | Step 3 output directory |
| | `--rules DIR` | `data/DefenseFinder_rules` | DefenseFinder rules directory |
| | `--output DIR` | `output/step4_raw/distance_filter` | Distance-filter output directory |
| | `--max_gene_gap N` | 10 | Maximum gene-index difference for a connectivity edge |
| | `--mode MODE` | `unspecified` | Log mode label |
| `run_cdhit_cluster.sh` | `--mode MODE` | – | Select the mode-specific Prodigal directory; required unless `--prodigal_dir` is given |
| | `--prodigal_dir DIR` | derived from `--mode` | Override the Prodigal protein tree |
| | `--threads N` | 12 | CD-HIT threads per process |
| `run_mining_esmfold.sh` | `--model_path PATH` | required | ESMFold model path |
| | `--gpus IDS` | required | Comma-separated GPU IDs |
| | `--max_seq_len N` | 3500 | Maximum sequence length |
| | `--thread_per_gpu N` | 1 | ESMFold model processes per GPU |
| `run_tmalign.sh` | `--mode MODE` | – | Select the Mining PDB directory; required unless `--mining_pdb` is given |
| | `--mining_pdb DIR` | derived from `--mode` | Override the Mining PDB directory |
| | `--output DIR` | `output/step4_output/tmalign` | TM-align working directory |
| | `--threads N` | 12 | Concurrent TM-align workers |
| `run_replication_homologous_frequency.sh` | `--step4_output DIR` | `output/step4_output` | Step 4 formal output directory |
| | `--diamond_db_dir DIR` | `data/gtdb/diamond_db` | GTDB DIAMOND database directory |
| | `--threads N` | 12 | Threads per DIAMOND worker |
| | `--workers N` | 4 | Concurrent DIAMOND workers |
| | `--shard_size N` | 128 | Unique queries per DIAMOND shard |
| | `--work_dir DIR` | `output/step4_raw/homologous_frequency` | Internal shard and cache directory |

## Key thresholds

| Parameter | Default | Rule |
|---|---|---|
| `--max_gene_gap` (S01) | 10 | Maximum gene-index difference for an edge in the anchor-connected component |
| CD-HIT `-c` (S03) | 0.5 | Sequence-identity threshold for instance clustering |
| CD-HIT `-n` (S03) | 3 | CD-HIT word length paired with `-c 0.5` |
| CD-HIT `-d` (S03) | 0 | Report the full sequence name in the `.clstr` file |
| `--tm_score` (S08) | 0.7 | Minimum Mining-anchor CDP-length-normalized TM-score for structural qualification |
| `--max_identity` (S08) | 0.3 | Maximum Mining-anchor aligned-region identity for structural qualification |
| S09 DIAMOND `-e` | 1e-5 | E-value for the homologous-frequency search |
| S09 DIAMOND `--id` | 50 | Percent identity for the homologous-frequency search |
| S09 DIAMOND `--query-cover` | 70 | Minimum query coverage for the homologous-frequency search |
| S09 DIAMOND `--subject-cover` | 70 | Minimum target coverage for the homologous-frequency search |
| `--min_count` (S09) | 5 | Minimum distinct contigs in the intersection of all component hit sets |

CD-HIT runs with 12 CPU threads (`-T 12`) and a 16000 MB memory cap (`-M 16000`)
by default; both are throughput and resource settings rather than scientific
thresholds.
`-d 0` is functional rather than a threshold: the clustering step reads the full
sequence names back from the `.clstr` file, and the 20-character CD-HIT default
would truncate long identifiers.

`--tm_score` and `--max_identity` (S08) and `--min_count` (S09) are Python
defaults; the Shell wrappers do not expose or override them.
S08 assigns structural qualification using the original Mining-anchor pair
only: both the TM-score and identity criteria must pass. Qualified instances
retain every assessed component row, including partners that fail the same
component-level criterion. Disqualified instances have an assessed anchor
that fails one or both criteria; missing structural assessments are separate.
The four anchor/partner conformity combinations can therefore occur in the
final collection.

S06 passes Mining first and Query second to TM-align. S07 reads the first
reported TM-score, normalized by the Mining (CDP) length. Sequence identity
is measured over aligned residues. Query-length coverage is a separate
quantity and does not redefine the stored TM-score or the qualification rule.

S09 runs DIAMOND in fixed `--sensitive` mode; its identity and bilateral
coverage thresholds are fixed in the Python operation rather than exposed as
chapter-level arguments.

S03 passes `--threads` directly to CD-HIT. S06 treats `--threads` as the number
of concurrent single-threaded TM-align workers. S09 uses sequence-level
deduplication and batches each domain into shards; `--threads` is the number of
threads per DIAMOND worker and `--workers` is the number of concurrent DIAMOND
workers, so the maximum requested CPU parallelism is `threads * workers`.

## Outputs

```text
output/step4_raw/
├── distance_filter/
│   ├── step4_distance_filtered_mapping.csv
│   ├── step4_distance_filtered_assembly.csv
│   └── step4_cleaning_log.txt
├── cdhit_cluster/
│   ├── clustered_mapping.csv
│   ├── clustered_assembly.csv
│   ├── qma_faa/
│   ├── qma_faa_clustered/
│   └── clustering_temp/                   # CD-HIT working files
│       ├── all_clusters.json
│       ├── sequence_id_mappings.csv
│       └── <system>/
│           ├── <component>_input.faa
│           ├── <component>_clusters
│           └── <component>_clusters.clstr
├── mining_pdb/                    # Replication ESMFold output
├── mining_pdb_status.csv          # Replication prediction status
└── homologous_frequency/          # Replication S09 work area
    ├── run_manifest.json          # Cache parameter and database signature
    ├── cache/<domain>/<sha256>.contigs
    └── shards/                    # Temporary, normally removed on success
output/step4_output/
├── tmalign/
│   ├── mining_assembly_pdb/
│   ├── query_assembly_pdb/
│   ├── tmalign_raw/                       # overlays plus canonical metrics
│   │   └── <system>/<context>/<instance>/
│   │       ├── <component>.pdb            # current TM-align overlay output
│   │       ├── <component>.pml            # visualization script
│   │       └── <component>.tmalign.txt    # captured TM-align stdout/stderr
│   └── comparison.csv
├── mapping_qualified.csv
├── assembly_qualified.csv
├── comparison_qualified.csv
├── mapping_disqualified.csv
├── assembly_disqualified.csv
├── comparison_disqualified.csv
└── homologous_frequency/                  # Replication only
    └── homologous_frequency_all.csv
```

`clustering_temp/` is the CD-HIT working area: one merged `_input.faa` per
component, its `.clstr` output, and `all_clusters.json`, whose
`system -> component -> sequence` key chain S04 consumes. The two flat tables in
that directory are internal working manifests keyed by directory-derived names
(`mining_id` is the instance directory, `subunit_name` is the component FAA stem)
rather than by the canonical Mapping/Assembly columns.

Internal shard files and the sequence-level result cache are stored under
`output/step4_raw/homologous_frequency/`; shard files are removed after their
results are cached, while the cache remains reusable on the next run.
`homologous_frequency_all.csv` reports one row per searched instance:

```text
defense_system,mining_domain,mining_genome_id,mining_contig_id,
mining_contig_system_index,homologous_assembly_count
```

Mapping and Assembly keep the canonical Step 3 schemas. Comparison columns:

```text
defense_system,query_genome_id,query_contig_id,query_contig_system_index,
query_component_protein_id,query_component,
mining_genome_id,mining_contig_id,mining_contig_system_index,
mining_component_protein_id,mining_component,is_candidate_mining,
tm_score,sequence_identity,aligned_length,rmsd
```

In Replication Mode the schema adds a `mining_domain` column after
`query_component`. Structurally unassessed instances (missing PDB or no valid
TM-align) stay out of the disqualified tables.

S06 captures the textual output of every successful TM-align call in
`<component>.tmalign.txt`; S07 reads TM-score, sequence identity, aligned length
and RMSD from that canonical file. A result is reusable only when the metric
file is parseable and a non-empty structural overlay (`<component>.pdb` for the
current binary or legacy `<component>_all_atm`) also exists. Legacy
`*_all_atm` files that contain all four metrics are converted in place to the
canonical metric file without recomputing the alignment. TM-align return-code,
metric-parsing or overlay failures stop the stage and do not write its `.done`
marker.

## Empty candidate sets

Zero input, distance-filtered or clustered Assembly rows end the chapter
normally before the dependent operations. After TM-align, qualified and
disqualified Assembly counts are reported separately. `qualified=0` with
`disqualified>0` is a valid scientific outcome and still proceeds to the
Replication frequency operation. A combined zero count from a non-empty
clustered Assembly input is instead treated as an execution error, because it
means that no candidate Mining component received a valid structural
comparison. Structurally unassessed instances are not relabelled as
disqualified, and an all-unassessed run cannot be marked complete.
