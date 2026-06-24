# SeCAT — Sequence Consensus Amplicon Trimming

SeCAT is a Nextflow pipeline that **harmonises 16S rRNA amplicon datasets across studies that used different primer pairs**. It aligns every study's ASVs to a common SILVA coordinate space, computes the longest 16S region shared by all studies, validates that trimming to that region does not distort each study's beta-diversity, and emits a unified feature table ready for cross-study meta-analysis.

If you only ever read one paragraph: edit two paths in `params.yaml`, run one Nextflow command, and look in `output/final_outputs/` for your results.

> **A note on the defaults.** The statistical thresholds shipped with SeCAT (`null_model_p_threshold`, `changepoint_penalty_multiplier`, `distance_cutoff_threshold`, `num_simulations`, and the simulation-engine parameters) are **calibrated through a systematic parameter-sensitivity study**, not arbitrary picks. They were selected to give optimal Type I error control at α = 0.05 across the seagrass meta-analysis test set. If you change any of them, report the deviation explicitly in your methods section — see [Parameter calibration](#parameter-calibration) below.

---

## Table of contents

1. [What SeCAT does, in five lines](#what-secat-does-in-five-lines)
2. [Before you start](#before-you-start)
3. [If you have never used an HPC before](#if-you-have-never-used-an-hpc-before)
4. [Installing the dependencies](#installing-the-dependencies)
5. [Preparing your inputs](#preparing-your-inputs)
6. [Running SeCAT](#running-secat)
7. [Understanding the output folder](#understanding-the-output-folder)
8. [The manual-review pause and the STANDARDIZE step](#the-manual-review-pause-and-the-standardize-step)
9. [Site-specific HPC configuration](#site-specific-hpc-configuration)
10. [Resource and runtime expectations](#resource-and-runtime-expectations)
11. [Troubleshooting](#troubleshooting)
12. [Parameter calibration](#parameter-calibration)
13. [Parameter reference](#parameter-reference)
14. [Citation and licence](#citation-and-licence)

---

## What SeCAT does, in five lines

1. **Cleans** each study (removes chloroplasts/mitochondria/empty samples).
2. **Maps** every ASV onto the SILVA aligned reference and finds each study's amplicon coordinates.
3. **Finds the consensus 16S window** — the longest stretch of SILVA covered by *all* studies at once.
4. **Tests by simulation** whether trimming each study down to that window distorts its community structure beyond what null-model noise would produce.
5. **Trims and merges** the surviving studies into a single feature table, taxonomy file, FASTA, and metadata table that can be analysed as one cohort.

---

## Before you start

You will need:

- A computer (laptop, server, or HPC login node) with **Linux or macOS**. Windows works only through WSL.
- Around **50 GB of free disk** for a 10-study analysis (most of it is Nextflow's work-cache; SeCAT's own outputs are <1 GB).
- **At least one machine with 64 GB RAM** to run the SILVA alignment step (laptops will not do this). On an HPC this is automatic; locally you need a workstation or a cloud VM.
- The **SILVA 138.x aligned reference FASTA**: `SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta`. This is a ~24 GB file. Download once from the SILVA website (https://www.arb-silva.de/) and keep it somewhere persistent — every SeCAT run reads it.
- One **ASV table, taxonomy file, FASTA and metadata table per study**. These are the outputs of DADA2 / QIIME2 / Deblur. The exact format SeCAT expects is shown in `examples/study_A/` and described in `examples/README.md`.

You do **not** need to know R or Nextflow to run SeCAT. You do need a working shell and the ability to edit a YAML file.

---

## If you have never used an HPC before

Skip this section if you already use SLURM (or SGE/LSF).

An HPC (High-Performance Computer) is a cluster of large servers shared by many users. You log in over SSH to a **login node**, but you **never run computation there** — instead you submit jobs to a **scheduler** (most commonly SLURM, sometimes SGE or LSF), which queues them and runs them on **compute nodes** that have lots of CPU and RAM. Nextflow handles all the job submission for you: from your point of view, you type one `nextflow run …` command on the login node, and Nextflow submits ~100 jobs over the next 12–24 hours on your behalf, retrying any that fail, gathering results, and writing them back to a shared filesystem.

Vocabulary you will see:

| Term | Meaning |
|---|---|
| **Login node** | The machine you SSH into. Don't run heavy work here. |
| **Compute node** | Where jobs actually run. You get to it via the scheduler. |
| **SLURM** | The most widely used HPC scheduler. Reads job scripts and queues them. (SGE and LSF are alternatives; SeCAT supports all three.) |
| **Partition / queue** | A pool of compute nodes (e.g. `short`, `standard`, `highmem`). Each partition has different limits. |
| **QoS** | "Quality of Service" — a priority tier inside a partition (some sites only). Common tier names: `short`, `standard`, `high`. |
| **`sinfo`** | (SLURM) Shows partition status (free / down). SGE equivalent: `qhost`. |
| **`squeue -u $USER`** | (SLURM) Shows YOUR queued and running jobs. SGE equivalent: `qstat -u $USER`. |
| **`scancel <jobid>`** | (SLURM) Cancels a job. SGE equivalent: `qdel`. |
| **`scratch`** | Fast but often time-limited filesystem. Many sites auto-delete files older than 28–90 days. Use for working directories. |
| **`home`** | Slower but permanent. Use for your code and small reference files. |
| **Module system** | Many HPCs require `module load <name>` to access compilers, Java, R, etc. Common environment modules: `module load java`, `module load r`, `module load singularity`. |

What you need to find out about *your* cluster before launching SeCAT:

1. **Scheduler type** — almost always SLURM in academic HPCs, occasionally SGE or LSF.
2. **Default partition / queue name** — ask your site admin or check `sinfo`. Common names: `compute`, `batch`, `normal`, `short`, `all.q`.
3. **Whether you need an `--account=` flag** — most shared university clusters require this for billing/quota.
4. **How to make Java ≥ 11 available** — Nextflow needs Java to run. Either `module load java` (if your site provides it), or install OpenJDK into your home directory.
5. **Whether Singularity is available** — `singularity --version` from the login node. If not, ask your admin (it's standard on academic HPCs).
6. **Scratch path conventions** — where can you write large temporary files, and what's the retention policy?

With those six answers you can populate `conf/custom.config` (described in [Site-specific HPC configuration](#site-specific-hpc-configuration) below) and launch SeCAT on any SLURM cluster.

> **Worked example — JASMIN LOTUS2** (UK NERC HPC, the pipeline's primary test environment). Skip if you're on a different site; the same six questions apply to yours.
>
> - Log in: `ssh <username>@login1.jasmin.ac.uk` (or `login2`, `login3`).
> - Make Java available: `module load jaspy`.
> - Default partition: `standard`. QoS tiers: `short`/`standard`/`high`. The `--qos=high` flag is required for any multi-CPU job — `short`/`standard` cap at 1 CPU per task.
> - Scratch: `/work/scratch-nopw2/<username>/` — auto-deleted after 28 days untouched.
> - Long-lived workspace: `/work/xfc/vol7/user_cache/<username>/` (30-day retention, 10 TB quota).
> - Pre-configured site file: `conf/jasmin.config` (included in this repo) — pass with `-c conf/jasmin.config`.

---

## Installing the dependencies

You need three pieces of software on the machine from which you will launch the pipeline (typically the HPC login node):

### 1. Nextflow (>= 25.10; 26.04 recommended)

Nextflow needs Java 11 or later. How you make it available depends on your site:

```bash
# Option A — use your site's module system (most academic HPCs):
module load java          # or `module load java/11`, `module load openjdk`, etc.
# Examples for specific sites:
#   JASMIN:    module load jaspy
#   ARCHER2:   module load PrgEnv-gnu  java
#   Cirrus:    module load openjdk

# Option B — install OpenJDK into your home directory (no module system, or older sites):
#   Download OpenJDK 17 from https://adoptium.net/, extract, and add to PATH.

# Then install Nextflow itself (system-agnostic):
curl -s https://get.nextflow.io | bash

# Move the resulting `nextflow` binary somewhere on your PATH:
mkdir -p ~/bin && mv nextflow ~/bin/
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

nextflow -version       # confirm install
```

### 2. A container runtime (Singularity OR Docker)

- **On an HPC** you almost certainly need **Singularity** (it does not require root, unlike Docker). On JASMIN it is available system-wide — no install needed. Just run `singularity --version` to confirm.
- **On a laptop or cloud VM** use **Docker Desktop** (or `apt install docker.io` on Ubuntu).

You do **not** need to build any container yourself. SeCAT pulls a pre-built one from GitHub Container Registry (`ghcr.io/deterrey/secat:latest`) on first launch, ~10 minutes.

### 3. The SeCAT pipeline itself

You don't actually install SeCAT — Nextflow can pull it directly from GitHub:

```bash
nextflow pull DETerrey/SeCAT
```

This places the pipeline under `~/.nextflow/assets/DETerrey/SeCAT/`. If you prefer a manual checkout (for editing or for offline use), clone it directly:

```bash
cd ~
git clone https://github.com/DETerrey/SeCAT.git
cd SeCAT
```

The instructions below assume the manual-checkout layout (run commands from inside the cloned `SeCAT/` directory). If you used `nextflow pull`, replace `main.nf` with `DETerrey/SeCAT` in every `nextflow run` command.

---

## Preparing your inputs

You need:

1. **A manifest TSV** — one row per study. Template: `assets/manifest_template.tsv`. Worked example: `examples/manifest_example.tsv` (with column descriptions in `examples/README.md`).
2. **The SILVA aligned reference FASTA** — download once from https://www.arb-silva.de/. Use the NR99 (non-redundant 99%) full-length aligned and truncated FASTA, e.g. `SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta`.
3. **An edited `params.yaml`** — minimally, change the two paths near the top:
   ```yaml
   manifest:     "/absolute/path/to/your/secat_manifest.tsv"
   reference_db: "/absolute/path/to/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta"
   ```

Before launching anything heavy, validate your inputs with the bundled healthcheck:

```bash
./secat_healthcheck.sh
```

The script verifies that all paths listed in your manifest exist and are readable. It is fast (no compute) and catches the most common configuration mistakes — much cheaper than discovering them after a 24-hour SLURM job has run.

See `examples/README.md` for the format requirements of each input file (feature table, taxonomy, FASTA, metadata).

---

## Container bind paths (Nextflow 26+)

SeCAT runs every step inside a Singularity/Docker container, which can only see
host folders explicitly mounted into it. Tell SeCAT where your data lives with
the `bind_paths` parameter in `params.yaml`:

```yaml
bind_paths: "/path/to/your/data"     # folder(s) with study files, metadata, reference_db
```

Point it at the folder(s) containing the paths in your manifest plus your
`reference_db`. One folder for everything? name it. Scattered? comma-separate
the roots (`"/data/studies,/refs/silva"`).

**Why:** Nextflow 26 enabled the strict config parser by default, which can no
longer auto-derive these paths from the manifest and does not pass
`SINGULARITY_BIND` through to cluster jobs. `bind_paths` puts the mounts on the
container command line — works on Nextflow 25 and 26.

## Running SeCAT

### Standard launch (HPC, SLURM + Singularity)

The minimum command — works on any SLURM cluster with a single default partition that doesn't require billing/account flags:

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    -params-file params.yaml
```

For most real clusters you'll also need a **site config** that maps SeCAT's generic resource labels (`mem_4g`, `mem_16g`, `mem_64g_cpu4`, etc.) to your cluster's partition names, account flags, and time limits. Pass it with `-c`:

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    -c conf/custom.config \
    -params-file params.yaml
```

`conf/custom.config` is a fully commented template — copy it, fill in your cluster's queue name and any account flags, and use it as your launch config from then on. See [Site-specific HPC configuration](#site-specific-hpc-configuration) for what to put in it.

> **Worked example — JASMIN LOTUS2**: a pre-filled site config (`conf/jasmin.config`) is shipped with this repo. Use it directly:
> ```bash
> nextflow run main.nf -profile slurm,singularity -c conf/jasmin.config -params-file params.yaml
> ```

### Local laptop / workstation (Docker)

```bash
nextflow run main.nf \
    -profile local,docker \
    -params-file params.yaml
```

Only suitable for very small test datasets (<5 studies, <50 simulations). The SILVA alignment step needs 64 GB RAM, which most laptops do not have.

### Resume a failed run

Add `-resume` to any of the above. Nextflow will re-use cached results for every process that completed successfully and only re-run the failed step:

```bash
nextflow run main.nf -profile slurm,singularity -params-file params.yaml -resume
```

`-resume` works even after a power cut or a manually-killed pipeline, as long as the `work/` directory and `.nextflow/` cache still exist. Do **not** delete those between launches if you might want to resume.

### Submitting Nextflow itself as a background process

On a login node, the safest pattern is to launch Nextflow inside `tmux` (or `screen`) so the pipeline survives a disconnection:

```bash
tmux new -s secat
# inside tmux:
module load java                                     # or your site's equivalent
nextflow run main.nf \
    -profile slurm,singularity \
    -c conf/custom.config \
    -params-file params.yaml \
    -resume
# detach with Ctrl-b then d
# reattach later with:  tmux attach -t secat
```

Most HPCs allow login nodes to host long-running shell sessions like `tmux` indefinitely, but check your site's policy — some clusters reap idle sessions after a few days. If yours does, use the scheduler to submit Nextflow itself as a long-running job (a "head job" pattern); ask your admin for the local convention.

---

## Understanding the output folder

After a **successful full run** (`auto_trim: true`, or after STANDARDIZE), everything you need is under `output/final_outputs/`. Start with the generated `README.txt`, then `1_report/`:

```
output/final_outputs/
  README.txt                            # START HERE - index to everything below
  1_report/                             # read this first
    MESAP_Master_Summary_Report.pdf     # <- headline diagnostic - open this first
    per_study/<study>_report.pdf        # per-study trim diagnostic, one per study
    figures/                            # report pages as PNGs (master/ + per_study/<study>/)
    run_diagnostics/                    # nextflow_report.html / timeline.html / dag.html
  2_dataset/                            # merged, trimmed dataset for downstream analysis
    combined_feature_table.tsv          # <- merged ASV table - feed to phyloseq / qiime2 / vegan
    combined_taxonomy.tsv
    combined_metadata.tsv
    combined_sequences.fasta
    asv_mapping_final.tsv               # crosswalk: original ASV -> harmonised ASV
  3_verdicts/                           # per-study keep/trim decisions
    master_verdict_table.csv            # per study x level: KEEP / EXCLUDE
    verdict_data_all_levels.csv         # diagnostic table behind the verdicts
    final_trim_verdicts.csv             # per-study trim coordinates applied
    trim_summary.csv                    # bp trimmed from 5'/3' per study
    simulation_baseline_statistics.csv  # null-model summary stats
    simulation_retention_curves.csv     # null-model degradation curves
  4_verification/                       # did trimming preserve the biology? plots + stats per level
    validation_summary.csv
    asv/  genus/  family/               # each level: diagnostic PDFs + permanova/network/concordance together
    logs/validation.log
  5_comparison/                         # transparency: pre- vs post-consensus tables
    pre_consensus/
    post_consensus/
  supporting/                           # provenance + dig-deeper material
    provenance/                         # params_used.yaml, manifest copies, selection_roster.txt
    coordinates/                        # per-study SILVA windows + global consensus region
    per_study_data/                     # per-study .rds, aligned FASTAs, ASV coords, standardized_fastas/
    taxon_impact/                       # flattened per-taxon CSVs (+ README.txt for the fields)
```

Folders are numbered in the order you'll usually want them: the write-up (`1_report/`), the data for analysis (`2_dataset/`), the decisions (`3_verdicts/`), the evidence (`4_verification/`), then transparency and supporting material. Nextflow's HTML run diagnostics are written straight into `1_report/run_diagnostics/`.

The working `output/` tree also holds intermediate folders (`cleaned_data/`, `intermediate/`, `real_data_results/`, `simulation_results/`, `meta_analysis/`, `standardized_datasets/`, `validation/`, `reports/`) - these are **deleted at end of run** unless `keep_intermediates: true`. Everything a user needs is mirrored into `final_outputs/` first, and the Nextflow `work/` cache is always preserved so `-resume` still works.

## The manual-review pause and the STANDARDIZE step

By default (`auto_trim: false`), SeCAT runs the first half of the pipeline (clean → map → simulate → verdict → reports) and then **pauses**. This is deliberate: you should review the per-study verdicts before committing to a 12-hour merge step.

The review loop is:

1. Wait for the first run to finish (you'll see "Pipeline complete. Review verdicts then run:" in the log).
2. Open `output/final_outputs/1_report/MESAP_Master_Summary_Report.pdf` and the per-study PDFs.
3. Look at `output/final_outputs/3_verdicts/verdict_data_all_levels.csv` — each row is a study × taxonomic level with a KEEP/EXCLUDE verdict.
4. Edit `selection_roster.txt` to list the studies you want in the final merged dataset.
5. Re-launch with `--step standardize -resume`:
   ```bash
   nextflow run main.nf \
       -profile slurm,singularity \
       -c conf/jasmin.config \
       -params-file params.yaml \
       --step standardize -resume
   ```

`STANDARDIZE` reads the existing intermediates (which is why cleanup is deferred until after a successful STANDARDIZE run), runs trimming + merging + validation, and writes the final results.

If you trust the automatic KEEP-verdicts and want to skip the manual pause, set `auto_trim: true` and `selection_mode: "auto"` in `params.yaml`. Not recommended for first runs on new datasets.

---

## Site-specific HPC configuration

SeCAT's `nextflow.config` declares **generic resource labels** (`mem_4g`, `mem_16g`, `mem_64g_cpu4`, etc.) describing what each process needs. A site config file maps those labels to your cluster's actual partition names, billing accounts, and time limits. This is what makes the pipeline portable across HPCs without modifying its source.

### Building a config for your SLURM cluster

Copy `conf/custom.config` and edit:

```groovy
process {
    // 1. Default partition / queue name on your cluster
    queue = 'compute'                       // ← replace with yours

    // 2. Account flag if your cluster requires billing/quota tracking
    // clusterOptions = '--account=YOUR_PROJECT_CODE'

    // 3. (Optional) Per-label overrides — only needed if certain memory tiers
    //    require a different partition (e.g. 'highmem' for >32GB jobs).
    // withLabel: 'mem_64g' {
    //     queue          = 'highmem'
    //     clusterOptions = '--account=YOUR_PROJECT_CODE --mem=64G'
    // }
}
```

Then launch with `-c conf/custom.config`. The file is heavily commented; read it once before editing.

### Worked examples for specific sites

**JASMIN LOTUS2** (UK NERC HPC) — a complete `conf/jasmin.config` is shipped with the repo. Notable settings: `queue = 'standard'`, `--account=no-project` (for unattached users — change to your project code if you have one), and `--qos=high` for the heavy `STUDY_MAPPING` step.

**A generic university SLURM cluster** with default `compute` partition and project billing:

```groovy
process {
    queue          = 'compute'
    clusterOptions = '--account=myproject'
}
```

**A SLURM cluster with a separate high-memory partition**:

```groovy
process {
    queue = 'compute'
    withLabel: 'mem_64g'        { queue = 'highmem' }
    withLabel: 'mem_64g_cpu4'   { queue = 'highmem' }
    withLabel: 'mem_48g_cpu16'  { queue = 'parallel' }
}
```

**SGE / Oracle Grid Engine clusters**: use `-profile sge,singularity` instead of `slurm,singularity`. No custom config usually needed — the `sge` profile in `nextflow.config` already maps every label to appropriate `clusterOptions` flags (`-pe smp N`, `h_vmem`, etc.).

**LSF clusters** (IBM Platform LSF): use `-profile lsf,singularity`. You may need to add a custom config with `process.queue` and any project-specific bsub flags.

### Singularity bind mounts

If your input data lives on a filesystem that Singularity doesn't auto-mount inside containers (common on sites with dedicated storage volumes), add a bind to your custom config:

```groovy
singularity.runOptions = '--bind /scratch/projects/myproject:/scratch/projects/myproject'
```

### Scheduler rate limits

If your cluster throttles rapid job submission (some do, to protect the scheduler), reduce Nextflow's submit rate:

```groovy
executor {
    queueSize       = 100
    pollInterval    = '30 sec'
    submitRateLimit = '20/1min'
}
```

---

## Resource and runtime expectations

For a representative meta-analysis (10–15 studies, 100 simulations, full validation):

| Stage | Memory | CPUs | Wall time |
|---|---|---|---|
| Cleaning | 4 GB | 1 | <30 min |
| Study mapping (DECIPHER alignment to SILVA) | 16–128 GB per study | 1–4 | 2–8 h per study (parallel) |
| Simulation (per replicate) | 18 GB | 1 | 30 min – 2 h |
| Real-data analysis | 64 GB | 1 | 2–4 h per study (parallel) |
| Aggregation | 64 GB | 1 | 2–4 h |
| Trim & merge | 30 GB | 1 | 4–12 h |
| Validation | 30 GB | 1 | 1–4 h |

Total wall time on JASMIN LOTUS2: typically **12–24 hours** end-to-end. Most of the elapsed time is queue wait, not compute.

---

## Troubleshooting

**"ERROR: --manifest is required"** — the pipeline cannot see your `params.yaml`. Make sure you passed `-params-file params.yaml` and that the file is in the working directory.

**"ERROR: Missing study files referenced in manifest"** — at least one path in the manifest does not exist or is unreadable. Run `./secat_healthcheck.sh` to identify which.

**Container pull fails (`Failed to pull singularity image`)** — your compute nodes have no outbound HTTPS access (common on secure clusters). Pre-pull the container from the login node and point Nextflow at the local image:
```bash
singularity pull docker://ghcr.io/deterrey/secat:latest
# then in your site config:
process.container = '/absolute/path/to/secat_latest.sif'
```

**Jobs sit `PENDING` in the queue forever** — three things to check, in order:
1. Is your cluster busy? Run `sinfo` (SLURM) or `qstat` (SGE) and look at partition load.
2. Is your `--account=` flag (in your site config) valid and within quota? Ask your admin if unsure.
3. Are you exceeding a concurrent-job or CPU limit set by your QoS tier? Reduce `executor.queueSize` in your site config (default 200; try 50).

**STUDY_MAPPING fails with out-of-memory** — Increase `reference_subset_size` slowly, or switch to `reference_alignment_mode: "full"` on a high-memory node. As a workaround, set `use_all_asvs: false` and reduce `asv_sample_size`.

**Too many studies flagged EXCLUDE** — your dataset may genuinely be too heterogeneous for the chosen consensus region, or your thresholds may be too strict. See [Parameter calibration](#parameter-calibration) before deviating from defaults; a calibration-aware first try is `distance_cutoff_threshold: 0.20`.

**Too few studies flagged EXCLUDE (over-permissive)** — drop `null_model_p_threshold` to 0.01 and increase `null_model_min_consecutive` to 5. Read the calibration section before bigger changes — over-tuning these dismantles the FPR safety net.

**`STANDARDIZE` complains it cannot find `output/intermediate/`** — the cleanup step deleted it. Set `keep_intermediates: true` and re-launch from the beginning, or run STANDARDIZE in the same Nextflow invocation as the verdict step (`auto_trim: true`).

**Nextflow exits silently after the version banner, or `-resume` fails with a misleading error** (e.g. "--manifest is required" when your `params.yaml` clearly sets it). Symptom: you launch the pipeline, get the `Nextflow ~ version 25.x.x` line, then the shell prompt returns immediately with no SeCAT banner, no errors, no SLURM jobs submitted. Cause: a corrupted Nextflow cache, usually from an earlier `^C`'d or SIGHUP'd run that left the `.nextflow/history` and `.nextflow/cache/` in an inconsistent state. Recovery:

```bash
# Nuke the Nextflow runtime state (NOT your output directory):
rm -rf ~/SeCAT/.nextflow/

# Re-launch WITHOUT -resume so Nextflow rebuilds its history from scratch:
nextflow run main.nf \
    -profile slurm,singularity \
    -c conf/jasmin.config \
    -params-file params.yaml
# (add --step standardize if you're past the verdict step)
```

What's safe to delete and what isn't:

- **Safe to delete**: `~/SeCAT/.nextflow/` (Nextflow runtime cache — re-created on every launch).
- **DO NOT delete**: `output/` (your actual results) or `/work/scratch-nopw2/<user>/secat_work/` (the per-task work directories — Nextflow can still pick these up as cached even without `.nextflow/history`).

After `rm -rf .nextflow/`, the first run won't be able to `-resume` (no history yet), so omit the flag. From the *next* run onwards, `-resume` works normally because Nextflow rebuilds its history as it goes. You haven't lost any actual data — only Nextflow's bookkeeping.

> **JASMIN-specific troubleshooting**
>
> - **Jobs pending forever**: check the `--account=` flag in `conf/jasmin.config`. The default `--account=no-project` works for unattached users; for project allocations, change it to your project code.
> - **`QOSGroupMaxCpu` error**: your QoS tier caps total concurrent CPUs. Reduce `executor.queueSize` in `conf/jasmin.config` from 200 to 50.
> - **`STUDY_MAPPING` needs `--qos=high`**: it's already set in `conf/jasmin.config`. If you've overridden this, the default `standard` QoS caps CPUs at 1 and will fail multi-CPU jobs.

---

## Parameter calibration

SeCAT's defaults are **not arbitrary**. They are the outcome of a systematic FPR/sensitivity calibration study (false-positive rate measured under a no-signal null; sensitivity measured under controlled-magnitude spike-ins). Three detectors are used in combination, and each has a distinct sensitivity/specificity profile that the defaults were tuned against:

| Detector | Parameter (default) | Calibrated behaviour | Role in the ensemble |
|---|---|---|---|
| **Null model** | `null_model_p_threshold` = 0.05, `null_model_min_consecutive` = 3 | Near-zero FPR across the parameter grid (FPR = 0 at `p ≤ 0.01, min_consecutive ≥ 3`). Near-zero standalone sensitivity. | High-specificity confirmatory signal — almost never fires falsely. |
| **Distance cutoff** | `distance_cutoff_threshold` = 0.15 | FPR ≈ 0.006 at 0.15; remains well below 0.05 across the 0.10–0.30 range. Fires at lag = 0 (exactly at the true signal). Sensitivity is approximately linear with spike size (~80% at Δ = 0.3, ~5% at Δ = 0.1). | **Primary detector** — well-calibrated and the workhorse of the ensemble. |
| **Changepoint (PELT)** | `changepoint_penalty_multiplier` = 50 (params.yaml) / 1 (nextflow.config) | FPR ≈ 0.70 at low multipliers; drops below 0.05 only at multiplier ≈ 500. At multiplier = 500, sensitivity falls to ~15% for spike Δ = 0.20. | Sensitivity-driven confirmatory signal — penalised lightly on purpose, relies on ensemble agreement to keep FPR controlled. |

The conservative-voting ensemble (≥ 2 of 3 methods must agree before a study is flagged for trimming) is the **load-bearing mechanism** for FPR control: no individual detector is calibrated for FPR < 0.05 *standalone*, but the ensemble exploits their complementary error modes. Lowering the changepoint multiplier in isolation will not break the pipeline, but it will erode the ensemble's specificity floor.

### Why the changepoint multiplier looks "wrong"

The change-point penalty in `params.yaml` is 50 — orders of magnitude below the ~500 needed to control PELT's FPR by itself. This is **deliberate**. Calibration showed that at multiplier = 500 the changepoint method retains only ~15% sensitivity for moderate community shifts (Δ-Bray-Curtis = 0.20), which would defeat its purpose as an early-warning signal. The current value keeps it sensitive at the cost of standalone specificity, and the ≥ 2-methods-agree rule provides the FPR control that the penalty alone cannot.

### If you change a default

You must do two things:

1. **Document it.** Record the change in your methods section alongside the rationale, ideally with a sensitivity analysis (re-run with the published default for comparison).
2. **Keep the provenance file.** `output/final_outputs/provenance/params_used.yaml` records every effective parameter for the run — supply it as supplementary material so reviewers can reproduce your call set.

### When deviation is defensible

- **Under-powered datasets** (< 5 studies, < 30 samples per study). The null model's near-zero standalone sensitivity becomes a problem when the ensemble is degraded by missing votes — consider relaxing to `null_model_min_consecutive = 2` and raising `null_model_p_threshold` to 0.10. The DistCutoff remains the safety net.
- **Strongly heterogeneous primer sets** (e.g. V1–V3 + V4 + V6–V8 in one analysis). Reducing study count from the consensus calculation is often more defensible than relaxing the detectors — but if you need to retain studies, `distance_cutoff_threshold = 0.20` admits ~10–15 more studies in our test set without crossing FPR = 0.05.
- **Publication-quality runs.** Increase `num_simulations` from 100 to 500 to tighten the null distribution; this does *not* change the calibrated penalty multiplier (the multiplier scales against data variance, not simulation count).

For anything else, the defaults represent the most defensible starting point and were tuned specifically for the seagrass 16S meta-analysis benchmark set. See `Manuscripts/SeCAT/` (in-prep manuscript) for the full calibration methodology.

---

## Parameter reference

The full annotated reference for every parameter lives in `params.yaml` (in the project root). The `nextflow.config` file documents the same parameters in their default form. The `examples/params_example.yaml` shows the minimal subset you would typically set.

The three knobs that most affect verdicts (and thus which studies enter your merged dataset) are:

| Parameter | Default | Effect |
|---|---|---|
| `null_model_p_threshold` | 0.05 | Alpha for the null-model significance test. Lower = stricter. |
| `changepoint_penalty_multiplier` | 1 | PELT changepoint penalty scale. Higher = fewer changepoints = more KEEP verdicts. |
| `distance_cutoff_threshold` | 0.15 | Bray-Curtis cutoff for the secondary verdict test. Lower = stricter. |

**Document any deviation from these defaults in your methods section** — they materially affect which studies you ultimately include.

---

## Citation and licence

If you use SeCAT in published work, please cite:

> Terrey, D., et al. *SeCAT: Sequence Consensus Amplicon Trimming for cross-study harmonisation of 16S rRNA microbiome datasets.* (Manuscript in preparation.)

Licence: MIT. See `LICENSE` for details.

Issue tracker: https://github.com/DETerrey/SeCAT/issues
