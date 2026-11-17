# Step 5: Information Integration

Step 5 resolves Mining component coordinates from the mode-specific Prodigal GFF
tree, orders the components of each Mining instance by transcriptional order, and
appends the derived ordering fields and a paper-style syntax string to the
Comparison tables. Structural qualification is inherited from Step 4 without
an additional structural threshold. Coordinate-complete, single-strand
instances receive the final syntax; instances missing component coordinates
or spanning more than one strand are dropped whole, as described under S02.
The coordinate and strand checks establish a same-strand arrangement, not
experimental evidence of co-transcription.

Candidate multi-gene defense systems (CMDSs) are named from their initial
construction in Step 3. Step 5 reports the coordinate-resolved, single-strand
subset of the Step 4 CMDSs, retaining their structural quality labels. It does
not introduce CMDS as a new qualification category. Assembly remains the name
of the table class and existing file interface.

## Complete run

Custom Mode:

```bash
bash Step5_Information_Integration/run_step5.sh --mode custom
```

Replication Mode:

```bash
bash Step5_Information_Integration/run_step5.sh --mode replication
```

`--mode` selects the mode-specific Prodigal tree as well as the mode's
completion markers, which are stored as empty files under
`output/.done/<mode>/step5/` (`01_locations`, `02_finalize`). Step 4 must have
produced non-empty `comparison_qualified.csv` and `comparison_disqualified.csv`;
the chapter reports `ERROR` and exits if they are missing or empty.

A completion marker is never trusted on its own. The location operation is
re-run whenever the existing location table no longer covers the component keys
referenced by the current Step 4 Comparison tables, so a Step 4 rerun that
introduces new instances cannot leave the coordinate table behind. The
final-tables operation always rewrites all six tables.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_step5.sh` | `--mode MODE` | required | Select the Custom or Replication input layout |
| | `--prodigal_dir DIR` | derived from `--mode` | Override the mode-specific Prodigal directory |
| | `--tmalign_dir DIR` | `output/step4_output` | Step 4 formal output directory |
| | `--output DIR` | `output/step5_output` | Final output directory |

The internal location table always lives at
`output/step5_raw/gene_location/component_locations.csv`; the chapter runner does
not expose an override for it.

## Stage commands

The two stage wrappers support step-by-step or recovery runs. Each runs one
operation and does not read or write chapter completion markers, which belong to
the chapter runner.

### Component locations (S01)

```bash
bash Step5_Information_Integration/run_component_locations.sh --mode custom
```

`--mode` selects the mode-specific Prodigal GFF tree: `output/prodigal_raw` for
Custom Mode and `output/rep_prodigal_raw` for Replication Mode. Pass
`--prodigal_dir` only to point it at a tree outside the standard location.

`S01_build_component_locations.py` reads both Comparison tables, collects every
`(mining_genome_id, mining_contig_id, mining_component_protein_id)` key they
reference, and looks each key up in the Prodigal GFF tree
(`<prodigal_dir>/<genome>/<contig>/<contig>.gff`). It keeps the `CDS` rows whose
`ID`, `locus_tag` or `gene_id` attribute matches and writes the GFF `start`,
`end` and `strand` of each match. Components without a matching `CDS` row, and
genomes or contigs that are absent from the tree, are omitted and reported as
`[S01] WARN unlocated_component=<n>`, so an incomplete coordinate table is never
silent. Independent GFF files are parsed concurrently (`min(CPU count, 8)`
workers), which is an internal throughput setting and not a scientific threshold.

The underlying Python entry point also accepts `--check`, which performs the key
comparison only: it scans no GFF file, writes nothing, and exits `0` when the
existing table already covers every key. The chapter runner uses it to decide
whether the operation can be skipped. The check validates the table against the
Step 4 component keys, not against the GFF tree, so after re-running Step 1 with
a different Prodigal tree, delete `output/.done/<mode>/step5/` to force a
rebuild.

### Final tables and syntax (S02)

```bash
bash Step5_Information_Integration/run_finalize_syntax.sh
```

`S02_finalize_syntax.py` copies the Mapping and Assembly tables from Step 4 and
appends two columns to both Comparison tables:

| Column | Meaning |
|---|---|
| `transcription_order` | 1-based rank of the component within its instance in transcriptional order |
| `syntax_representation` | Per-instance component string in the paper nomenclature, repeated on every row of the instance |

Instances are the Mining instance key used throughout Steps 3-4
(`defense_system`, `mining_genome_id`, `mining_contig_id`,
`mining_contig_system_index`).

**Instance exclusion.** The current syntax describes single-strand component arrangements, so an
instance is dropped **whole** from all six tables when either

- any component has no row in `component_locations.csv`; or
- its located components span more than one strand (mixed strand).

The whole instance is removed, not the offending component, and the run reports
`[S02] WARN dropped_instances=<n> missing_coordinates=<a> mixed_strand=<b>`.
S01's `unlocated_component` warning is the earlier signal that coordinates are
missing upstream.

Whole-instance exclusion keeps Mapping, Assembly and Comparison consistent
and preserves the complete component set. A single `[Chain]` value describes
the coding strand of the retained same-strand arrangement. Mixed-strand
instances are outside this syntax's current reporting scope.

**Transcriptional order.** The components of an instance are ranked
forward-strand first, then reverse-strand: forward-strand components ascending by
GFF `start`, followed by reverse-strand components descending by `start`.
`transcription_order` is the 1-based rank in that order.

**Syntax.** Each component is lettered `A`, `B`, `C`… in transcriptional order.
The start point is the discovery subunit (`is_candidate_mining == "T"`); if an
instance has no such row, the first transcribed component is used. Subunits
meeting the isozyme standard (`tm_score >= 0.7` and `sequence_identity <= 0.3`)
are marked with `*`; the chain direction in `[+]`/`[-]` is taken from the **start
point**:

```text
(Pre-group)*-Start-(Post-group)*[Chain]
```

The start point receives `*` when its own pair satisfies the criterion and
remains unstarred otherwise. Components before and after it in coding-direction
order are wrapped in parenthesised groups, each followed by one `*` per
conforming component in that group. Group stars report the number of
conforming components, not necessarily that every member conforms. The
Comparison rows retain the individual component metrics. Empty groups are
omitted. The `tm_score` criterion uses the CDP-length-normalized value stored
by Step 4; aligned-region identity and Query-length coverage are distinct
metrics. `[+]` or `[-]` records coding direction, not measured transcription.

The operation owns all six output tables and removes them before writing the new
ones, so an interrupted run cannot leave a complete-looking stale set behind.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_component_locations.sh` | `--mode MODE` | `unspecified` | Select the default Prodigal GFF tree (`custom` or `replication`) |
| | `--step4 DIR` | `output/step4_output` | Step 4 formal output directory holding both Comparison tables |
| | `--prodigal_dir DIR` | derived from `--mode` | Override the Prodigal GFF tree |
| | `--output FILE` | `output/step5_raw/gene_location/component_locations.csv` | Component-location table |
| `run_finalize_syntax.sh` | `--mode MODE` | `unspecified` | Log mode label |
| | `--step4 DIR` | `output/step4_output` | Step 4 formal output directory |
| | `--locations FILE` | `output/step5_raw/gene_location/component_locations.csv` | Component-location table written by S01 |
| | `--output DIR` | `output/step5_output` | Final output directory |

Each wrapper accepts `--mode custom` or `--mode replication`; without it the
standard path under `output/` is used. Running a stage wrapper directly neither
reads nor writes chapter completion markers.

## Key thresholds

Step 5 introduces no additional structural threshold: the `tm_score`/`sequence_identity`
criteria and the homologous-frequency filter are applied in Step 4. The
`tm_score >= 0.7` / `sequence_identity <= 0.3` values used for the syntax
conformance marker are a re-reading of the same Step 4 criteria, not a new
filter. The only instances Step 5 removes are those it cannot describe: an
instance missing a coordinate row, or one spanning more than one strand.

## Outputs

```text
output/step5_raw/
└── gene_location/
    └── component_locations.csv
output/step5_output/
├── mapping_qualified.csv
├── assembly_qualified.csv
├── comparison_qualified.csv
├── mapping_disqualified.csv
├── assembly_disqualified.csv
└── comparison_disqualified.csv
```

`output/step5_output/` holds the final pipeline deliverables. No further chapter
of the public workflow consumes them.

### component_locations.csv (internal)

One row per located Mining component, in GFF traversal order (`rglob` order, then
line order within each GFF); the table is not sorted by coordinate.

```text
mining_genome_id,mining_contig_id,mining_component_protein_id,start,end,strand
```

`start` and `end` are the 1-based inclusive GFF coordinates kept as text;
`strand` is `+` or `-`. The table is a lookup keyed by
`(mining_genome_id, mining_contig_id, mining_component_protein_id)`.

### mapping_{qualified,disqualified}.csv

The Step 4 Mapping rows for the surviving instance set (instances removed by the
exclusion rules above are absent), with the canonical Mapping schema:
`defense_system`, the five Query columns, then the six Mining columns
(12 columns in Custom Mode). Replication Mode inserts `mining_domain` after the
Query columns (13 columns):

```text
defense_system,query_genome_id,query_contig_id,query_contig_system_index,
query_protein_id,query_component_count,[mining_domain,]
mining_genome_id,mining_contig_id,mining_contig_system_index,
mining_protein_id,mining_component,assembly_type
```

### assembly_{qualified,disqualified}.csv

The Step 4 Assembly rows for the surviving instance set. The Assembly schema is
variable-width: a fixed header block followed by a
`mining_component_protein_id_<i>` / `mining_component_<i>` pair per component
slot, so the column count grows with the largest populated component index.

```text
defense_system,[mining_domain,]mining_genome_id,mining_contig_id,
mining_contig_system_index,mining_component_count,assembly_type,
mining_component_protein_id_1,mining_component_1,
mining_component_protein_id_2,mining_component_2,...
```

### comparison_{qualified,disqualified}.csv

The Step 4 Comparison schema with `transcription_order` and
`syntax_representation` appended at the end (18 columns in Custom Mode, 19 in
Replication Mode):

```text
defense_system,query_genome_id,query_contig_id,query_contig_system_index,
query_component_protein_id,query_component,[mining_domain,]
mining_genome_id,mining_contig_id,mining_contig_system_index,
mining_component_protein_id,mining_component,is_candidate_mining,
tm_score,sequence_identity,aligned_length,rmsd,
transcription_order,syntax_representation
```

An instance contributes one Comparison row per component, and every row of an
instance carries the same `syntax_representation`. Instance membership follows
Step 4: an instance qualified by its candidate component keeps all of its
component rows in `comparison_qualified.csv`. An instance removed by the
exclusion rules is absent from all three tables of its status, so the six tables
never disagree on which instances they describe.

## Empty candidate sets

The chapter counts qualified plus disqualified input Assembly rows. A combined
zero reports `STOP` and exits normally before coordinate extraction. Qualified
being empty alone does not stop information integration, and structurally
unassessed instances are never relabelled as disqualified. Missing or empty
Step 4 Comparison tables remain errors. The row count of the location table is
not a candidate count: an instance that is incomplete or mixed-strand is dropped
from the final tables and reported as
`[S02] WARN dropped_instances=<n> missing_coordinates=<a> mixed_strand=<b>`,
not silently kept.
