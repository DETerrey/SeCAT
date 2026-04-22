# SeCAT — Sequence Consensus Amplicon Trimming

SeCAT harmonises multi-primer 16S rRNA amplicon datasets for cross-study
meta-analysis. It aligns study ASVs to a common SILVA reference space,
detects the consensus trimming region, validates trim decisions against
null model simulations, and produces a unified feature table ready for
downstream community ecology analyses.

Developed for *Zostera marina* rhizosphere microbiome research at the
University of Derby. Applicable to any multi-primer 16S amplicon
meta-analysis.

---

## Requirements

- [Nextflow](https://nextflow.io/) ≥ 23.04 (`curl -s https://get.nextflow.io | bash`)
- [Singularity](https://sylabs.io/singularity/) ≥ 3.8 **or** Docker ≥ 20
- SILVA 138 full-length aligned reference FASTA
  (`SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta`)

No other software installation is required. All R packages, VSEARCH, and
MAFFT are bundled in the container, which is pulled automatically from
GitHub Container Registry on first run.

---

## Quick start

**1. Prepare your manifest**

Copy `assets/manifest_template.tsv` and fill in one row per study:

```
study_name    otu_table    taxonomy    metadata    fasta
Study_A       ...          ...         ...         ...
```

**2. Edit `params.yaml`**

Set the two required fields:

```yaml
manifest:     "secat_manifest.tsv"
reference_db: "/path/to/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta"
```

**3. Run**

SLURM (generic):
```bash
nextflow run DerbyDT/SeCAT \
    -profile slurm,singularity \
    -params-file params.yaml
```

SGE:
```bash
nextflow run DerbyDT/SeCAT \
    -profile sge,singularity \
    -params-file params.yaml
```

Local (testing):
```bash
nextflow run DerbyDT/SeCAT \
    -profile local,singularity \
    -params-file params.yaml
```

Resume after failure:
```bash
nextflow run DerbyDT/SeCAT -profile slurm,singularity -params-file params.yaml -resume
```

---

## Site-specific HPC configuration

SeCAT uses generic resource labels so it runs on any SLURM or SGE cluster
without modification. Whether you need a site config depends on your cluster:

| Situation | What to do |
|-----------|-----------|
| Single default queue, no account flag | No config file needed |
| Multiple partitions (e.g. highmem) | Use `conf/custom.config` |
| Account/project billing required | Use `conf/custom.config` |
| JASMIN LOTUS2 | Use `conf/jasmin.config` |
| SGE cluster | Use `-profile sge,singularity` — no config usually needed |

**For most SLURM clusters**, copy `conf/custom.config`, set your queue name,
and pass it with `-c`:

```bash
nextflow run DerbyDT/SeCAT \
    -profile slurm,singularity \
    -c conf/custom.config \
    -params-file params.yaml
```

The config file contains a full list of the resource labels SeCAT uses
(memory, CPUs, time limits) so you know exactly what each process requests.

**For JASMIN LOTUS2**, a pre-configured file is provided:

```bash
nextflow run DerbyDT/SeCAT \
    -profile slurm,singularity \
    -c conf/jasmin.config \
    -params-file params.yaml
```

**If your data is on a non-standard filesystem** (not automounted inside
Singularity), add a bind mount to your config file:

```
singularity.runOptions = '--bind /path/to/your/data:/path/to/your/data'
```

---

## Pipeline stages

| Stage | Process | Key output |
|-------|---------|------------|
| 0 | Data cleaning & QC | Filtered ASV tables, taxonomies |
| 1 | Primer/study coordinate mapping | SILVA-aligned coordinates |
| 2 | Consensus region calculation | Global overlap coordinates |
| 3 | Simulation preparation | Task matrix, reference subset |
| 4 | Null model generation | Simulated degradation baselines |
| 5 | Real data trimming analysis | Study-level β-diversity curves |
| 6 | Aggregation & verdict table | `verdict_data_all_levels.csv` |
| 7 | Report generation | Per-study PDF reports |
| 8 | Study selection | `selected_studies_for_trim.txt` |
| 9 | Sequence trimming | `*_standardized.fasta` |
| 10 | Dataset merge | Combined feature table, taxonomy, FASTA |
| 11 | Validation | Multi-tier ecological QC report |

By default the pipeline pauses after stage 7 for manual verdict review.
Set `auto_trim: true` in `params.yaml` to run stages 8–11 automatically.

After manual review, continue with the STANDARDIZE sub-workflow:

```bash
# 1. Review output/aggregated_data/verdict_data_all_levels.csv
# 2. Select studies
Rscript R/11_select_studies.R

# 3. Resume trimming + merge + validation
nextflow run DerbyDT/SeCAT \
    -profile slurm,singularity \
    -params-file params.yaml \
    --entry STANDARDIZE \
    -resume
```

---

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `manifest` | — | Path to study manifest TSV (required) |
| `reference_db` | — | Path to SILVA aligned FASTA (required) |
| `num_simulations` | 100 | Simulations per study per taxonomic level |
| `changepoint_penalty_multiplier` | 1 | PELT penalty multiplier for degradation detection |
| `null_model_p_threshold` | 0.05 | Empirical p-value threshold for null model |
| `distance_cutoff_threshold` | 0.15 | Bray-Curtis Δ threshold for distance cutoff method |
| `auto_trim` | false | Skip manual verdict review |
| `run_validation` | true | Run multi-tier ecological validation after merge |

Full parameter reference: see `params.yaml` and `nextflow.config`.

---

## Output structure

```
output/
├── aggregated_data/
│   ├── verdict_data_all_levels.csv     # trimming verdict per study/level
│   └── master_verdict_table.csv
├── reports/
│   └── *.pdf                           # per-study diagnostic plots
├── standardized_datasets/
│   └── *_standardized.fasta
├── meta_analysis/
│   ├── combined_feature_table.tsv
│   ├── combined_taxonomy.tsv
│   ├── combined_sequences.fasta
│   └── combined_metadata.tsv
└── validation/
    └── outputs/                        # multi-tier validation results
```

---

## Citation

Terry D. et al. (in preparation). SeCAT: a bioinformatics pipeline for
cross-study harmonisation of 16S rRNA amplicon datasets.

---

## Licence

MIT
