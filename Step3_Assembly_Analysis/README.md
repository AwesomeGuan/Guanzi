# Step 3: Assembly Analysis

Step 3 builds dictionaries from the candidate-contig DefenseFinder output and
assembles candidate multi-gene defense systems (CMDSs) through four ordered strategies:
independent novel-hit assembly, partial complementation, complete-system
expansion and Query-guided homologous assembly.

CMDS names the candidate system instance from its construction in Step 3;
subsequent quality control assigns retention and structural-qualification
states. The term does not require Step 4 qualification or demonstrated defense
function. Each instance contains multiple component genes. Instance counts are
separate from defense-system type counts and source-genome counts.

Use CMDSs for the candidate objects and their counts, and assembly for the
construction operation, strategies and data class. Existing directory names,
`assembly.csv`, Assembly tables and `assembly_type` retain their interfaces.

## Complete run

Custom Mode:

```bash
bash Step3_Assembly_Analysis/run_step3.sh --mode custom
```

Replication Mode:

```bash
bash Step3_Assembly_Analysis/run_step3.sh --mode replication
```

`--mode` selects the mode's completion markers and default protein tree, and
stops the chapter normally when the Step 2 multi-component input is empty. Empty
markers are stored in `output/.done/<mode>/step3/`.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_step3.sh` | `--mode MODE` | required | Select the Custom or Replication input layout |
| | `--threads N` | 12 | DIAMOND threads for Query-guided assembly (S07) |

## Stage commands

The lower-level entry points support step-by-step or recovery runs. The chapter
entry runs them in order and additionally manages completion markers.

### Dictionaries (S01–S03)

```bash
bash Step3_Assembly_Analysis/run_step3_dict.sh
```

The dictionaries are built from the shared Step 2 DefenseFinder results, so the
dictionary step is identical for both modes and takes no mode argument.

### Assembly (S00, S04–S07)

Custom Mode:

```bash
bash Step3_Assembly_Analysis/run_step3_assembly.sh --mode custom
```

Replication Mode:

```bash
bash Step3_Assembly_Analysis/run_step3_assembly.sh --mode replication
```

`--mode` selects the default protein tree used by Query-guided assembly. Pass
`--prodigal` only to point it at a protein tree outside the standard location.

### Optional arguments

| Script | Argument | Default | Description |
|---|---|---|---|
| `run_step3_dict.sh` | `--input DIR` | `output/step2_raw/defensefinder_results` | DefenseFinder results tree |
| | `--output DIR` | `output/step3_raw` | Dictionary output directory |
| `run_step3_assembly.sh` | `--mode MODE` | `unspecified` | Select the default protein tree (`custom` or `replication`) |
| | `--prodigal DIR` | derived from `--mode` | Override the protein tree used by Query-guided assembly |
| | `--input DIR` | `output/step2_output/step2_mapping_multi` | Step 2 multi-component mapping directory |
| | `--dict_dir DIR` | `output/step3_raw` | Dictionary directory written by S01–S03 |
| | `--rules DIR` | `data/DefenseFinder_rules` | DefenseFinder rule directory |
| | `--query_location DIR` | `data/query_data/query_location` | Normalized Query-location directory |
| | `--diamond_bin PATH` | auto-detected | Optional DIAMOND executable override |
| | `--threads N` | 12 | DIAMOND threads for S07 |
| | `--output DIR` | `output/step3_output` | Assembly output directory |

The `--prodigal` default is `output/prodigal_raw` for `--mode custom` and
`output/rep_prodigal_raw` for `--mode replication`. Only Query-guided assembly
(S07) reads the protein tree. Without `--mode` or `--prodigal`, the wrapper exits
with an error.

## Key thresholds

Query-guided assembly (S07) searches Query proteins against the Mining proteome
with DIAMOND `--sensitive` and retains:

| Parameter | Default | Rule |
|---|---|---|
| `--min_identity` | 50.0 | Retain hits with percent identity at or above the value |
| `--min_qcov` | 70.0 | Retain hits with Query coverage at or above the value |
| `--min_tcov` | 70.0 | Retain hits with target coverage at or above the value |
| `--evalue` | 1e-5 | Retain hits with E-value at or below the value |

These are S07 Python defaults; the Shell wrappers do not override them.

## Outputs

```text
output/step3_raw/
├── defensefinder_genes_dict.pkl
├── defensefinder_systems_dict.pkl
├── defensefinder_incomplete_dict.pkl
└── step2_merged_mapping.csv
output/step3_output/
├── independent_novel_hit_assembly/{mapping.csv,assembly.csv,remain_for_partial_complementation.csv}
├── partial_complementation/{mapping.csv,assembly.csv,remain_for_complete_expansion.csv}
├── complete_expansion/{mapping.csv,assembly.csv,remain_for_query_guided_assembly.csv}
└── query_guided_assembly/{mapping.csv,assembly.csv,remain_unmatched.csv}
```

Mapping columns:

```text
defense_system,query_genome_id,query_contig_id,query_contig_system_index,
query_protein_id,query_component_count,
mining_genome_id,mining_contig_id,mining_contig_system_index,
mining_protein_id,mining_component,assembly_type
```

Assembly columns:

```text
defense_system,mining_genome_id,mining_contig_id,mining_contig_system_index,
mining_component_count,assembly_type,
mining_component_protein_id_1,mining_component_1,
mining_component_protein_id_2,mining_component_2,...
```

Assembly holds one row per Mining CMDS instance, and component columns expand
in pairs `mining_component_protein_id_N` / `mining_component_N`. In Replication
Mode the schema adds a `mining_domain` column: after `query_component_count` in
Mapping and after `defense_system` in Assembly.

## Empty candidate sets

Zero multi-component input rows end the chapter normally before dictionary
building. Inside the assembly wrapper, the remaining CSV data rows decide
whether the next strategy runs; a header-only remain file counts as zero, the
later strategies are skipped, and earlier assemblies stay available. The skipped
strategy directories are still created for the chapter interface. The final
summary counts Assembly rows across all four strategies and does not treat
remaining candidates as successful assemblies.
