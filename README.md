# SeCAT — Sequence Consensus Amplicon Trimming

**Version 5.0.6**

SeCAT is a Nextflow DSL2 pipeline that harmonises independently published 16S rRNA
amplicon datasets sequenced with *different primer pairs* into a single feature table of
exact sequence units (MetaASVs) that are directly comparable across studies. It
establishes, from the sequences alone, that studies read the same stretch of the gene:
every study's ASVs are mapped into the fixed SILVA 138.2 alignment frame, a consensus
window shared by the largest mutually overlapping study set is determined, simulation-
calibrated degradation verdicts quantify what trimming each study to that window would
cost, and only then are sequences trimmed, merged at exact identity and validated.

## How it works

The pipeline runs in two phases with a deliberate pause between them.

**Phase 1 — map, model, verdict** (CLEAN_DATA, STUDY_MAPPING, COLLECT_MAPS,
PREPARE_SIMS, RUN_SIMULATION, ANALYSE_REAL, AGGREGATE, GENERATE_VERDICTS,
GENERATE_REPORT): input QC, coordinate mapping onto the SILVA 138.2 alignment
(DECIPHER profile alignment), consensus-window determination via an overlap
compatibility graph, per-study simulation of a structure-matched null community with
PCR and sequencing artefact layers, three independent degradation detectors, a
consensus vote, and a per-study PDF report. The run then pauses for review.

**Phase 2 — standardise** (TRIM_SEQUENCES, MERGE_DATASETS, VALIDATE; rerun with
`--step standardize -resume`): the selected studies are trimmed to the consensus
window, collapsed at exact sequence identity into MetaASVs with a composite-
score taxonomy, and the harmonised dataset is validated against the pre-harmonisation
data across alpha diversity, ordination geometry, dispersion, rank abundance,
taxonomic composition and co-occurrence network structure. Validation runs by default
and can be disabled with run_validation: false.

## Requirements

- Nextflow >= 25.10 (26.04 recommended), DSL2
- Singularity/Apptainer (recommended on HPC) or Docker
- The published container: ghcr.io/deterrey/secat:5.0.6
  (R 4.3.2, DECIPHER 2.30.0, changepoint 2.3, vegan 2.7-5, microeco 2.3.0,
  SpiecEasi 1.1.2, VSEARCH 2.23.0, and dependencies; the image build pins a dated
  CRAN snapshot and Bioconductor 3.18, so these versions are reproducible)
- The SILVA 138.2 SSURef NR99 **full-length aligned** FASTA as the reference

R and all R packages ship inside the container — the host needs only Nextflow and a
container runtime.

## Quick start

    nextflow pull DETerrey/SeCAT -r v5.0.6
    # fetch the parameter template (nextflow pull stores the repo under $NXF_HOME/assets):
    cp ~/.nextflow/assets/DETerrey/SeCAT/params.yaml my_params.yaml
    # edit my_params.yaml: set manifest, reference_db, outdir, bind_paths
    nextflow run DETerrey/SeCAT -r v5.0.6 -profile slurm,singularity -params-file my_params.yaml
    # review 1_report/ and the verdict tables, choose studies (see Study selection), then:
    nextflow run DETerrey/SeCAT -r v5.0.6 -profile slurm,singularity -params-file my_params.yaml \
        --step standardize -resume

Executor profiles: slurm, sge, local. Container profiles: singularity, docker.
**The singularity profile requires `bind_paths`** (comma-separated host paths covering
your data, reference and output locations); if it is unset, every task fails at
submission with a Singularity bind error rather than a named pipeline message. Use an
**absolute `outdir`** — a relative one is resolved inconsistently between processes
(some resolve it against the pipeline directory, not your launch directory).

## Inputs

SeCAT uses a tab-separated manifest with one row per study and exactly these
10 columns (template in assets/manifest_template.tsv; full column reference in
examples/README.md):

    study_name  primer_name  initial_fwd_trim  initial_rev_trim  asv_fasta_path
    asv_counts_path  taxonomy_path  metadata_path  metadata_variable  metadata_value

For each study the manifest points to four files: the ASV representative-sequence
FASTA, a pre-processed feature (ASV x sample) count table, taxonomy assignments, and
sample metadata. The last two manifest columns are an optional per-study sample
filter — samples are kept where metadata_variable matches metadata_value, which may
be a single value or a comma-separated list. Leave them empty to keep all samples.

The reference alignment (SILVA 138.2 SSURef NR99, full-length aligned) defines the
common coordinate frame for ASV mapping across studies, so changing the reference
version invalidates any previously computed consensus — all studies in a synthesis
must be mapped against the same reference.

The first process, CLEAN_DATA, prepares inputs by removing chloroplast and
mitochondrial sequences, dropping samples below the read-depth threshold, and then
filtering ASVs by abundance and prevalence (computed over the surviving samples);
studies with fewer than three samples are excluded, checked on the metadata before
depth filtering. Sample identifiers are reconciled
between the metadata and the feature table locally (BioSample identifiers in the
metadata are matched to run accessions in the table headers), so either identifier
style works; a further online SRA-run-to-BioSample lookup via NCBI E-utilities runs
at the merge stage only if fewer than half the identifiers reconcile. The feature
table, metadata and FASTA are then intersected to a common ASV and sample set.

## Key parameters

The shipped params.yaml is a template carrying the defaults used in the SeCAT paper.
Detection and simulation parameters are exported by the modules and read by the R
stages end-to-end. (Two legacy knobs — asv_sample_size and max_primer_mismatch — are
currently not exported and fall back to internal defaults; they affect no published
result and will be wired through or removed in a future release.)

- **Input QC (CLEAN_DATA)**: min_sample_depth (10000), min_asv_prevalence (2),
  min_asv_reads (50); remove_organelles (true) drops chloroplast/mitochondrial
  ASVs and can be disabled when organellar sequences are the study target; the
  three-samples-per-study floor and a feature-table/metadata overlap check are
  always applied. The trim/merge (standardize) step builds the harmonised
  dataset from CLEAN_DATA's outputs, so these filters carry through to the
  final tables.
- **Coordinate mapping**: mapping_support_band (250 columns), mapping_support_warn
  (0.70), mapping_support_fail (0.40). With exclude_unstable_from_consensus: true
  (default), UNSTABLE studies are removed from the analysis set entirely: they take
  no part in consensus determination and receive no simulation, verdict or report.
- **Simulation null**: one simulated community per study and replicate, built as a
  *structure-matched* null (sim_structure_matched_null: true, cluster identity 0.97)
  so the null degrades realistically. Artefact layers: GC-dependent PCR bias
  (quadratic penalty, strength 0.65, floored at 10%, 25 cycles), position-dependent
  sequencing error (0.003 at 5', doubling by 3'; transitions favoured
  0.70:0.15:0.15), optional chimeras (sim_add_chimeras: false).
- **Degradation detectors**: PELT changepoint with the penalty selected per
  study x level by a bootstrap CV scan (changepoint_scan_min–max, 0.4–2.4, x log(n),
  changepoint_bootstrap_n 50, changepoint_seed 10010; changepoint_penalty_multiplier
  is legacy and affects only per-study diagnostics); step-wise Bray–Curtis cutoff
  distance_cutoff_threshold (0.20); null-model empirical p < null_model_p_threshold
  (0.05) for null_model_min_consecutive (2) consecutive steps.
- **Verdicts**: >= 2 detectors agreeing within verdict_consensus_tolerance_bp
  (500 columns) -> CONFIRMED; one -> WARNING_SINGLE; hard gates
  (bc_ceiling_threshold 0.20, retention_floor_threshold 0.70) override to
  HARD_FAIL. A verdict is not an inclusion decision: detected thresholds are
  compared with each study's required trim at study selection.

## Study selection (the pause)

The pause between phases is a safe stopping point: all Phase 1 outputs are written to
disk, and standardisation can be launched hours or days later without rerunning any
completed work. Phase 2 reads Phase 1's published outputs from outdir — the
intermediate/ tree (consensus and coordinate files, aligned FASTAs), aggregated_data/
and cleaned_data/ — so do not delete these between phases.

After Phase 1, review the per-study reports in 1_report/per_study/,
SeCAT_Master_Summary_Report.pdf and the verdict tables in 3_verdicts/. Then write the
studies to carry forward into a roster: a plain text file, one study_name per line
(names must match the manifest exactly; lines starting with # are ignored), and set
in your params file

    selection_mode: roster
    selection_file: /absolute/path/to/selection_roster.txt

before rerunning with --step standardize -resume. Phase 2 always requires
selection_file (it is checked at launch whatever the mode), and in roster mode the
path must be absolute and, under the singularity profile, covered by bind_paths,
because it is re-read inside the container. The roster keeps every inclusion decision
explicit and auditable, and a copy is archived with the run's provenance.

## Running on HPC

Each executor profile (slurm, sge, local) combines with a container profile
(singularity, docker). Per-process resource labels map to scheduler directives, so
queue/memory/time can be tuned per site in conf/ without touching the pipeline.
RUN_SIMULATION fans out one single-core task per study x replicate and is throttled
with maxForks (set in the process block; override per site with -c). Four processes
share the top (64 GB) memory tier — STUDY_MAPPING (the per-study DECIPHER
alignment), ANALYSE_REAL (which also fans out per study, so several such tasks can
run concurrently), AGGREGATE and GENERATE_REPORT — worth knowing when sizing queues. Host directories are exposed to the container through bind_paths.
Nextflow's -resume reuses cached work after interruptions, and the named STANDARDIZE
workflow is dispatched with the step parameter.

## Outputs

    outdir/final_outputs/
      1_report/        per_study/ PDF reports, SeCAT_Master_Summary_Report.pdf,
                       figures/, run_diagnostics/ (Nextflow report, timeline, DAG)
      2_dataset/       harmonised MetaASV feature table, taxonomy, metadata, FASTA,
                       ASV->MetaASV map (phyloseq-importable)
      3_verdicts/      master_verdict_table.csv, verdict_data_all_levels.csv,
                       trim_summary.csv, excluded_studies.csv, simulation
                       baseline/retention curves
      4_verification/  validation outputs per level (asv/, genus/, family/) +
                       validation_summary.csv
      5_comparison/    pre_consensus/ and post_consensus/ tables (cleaned untrimmed
                       vs harmonised) for before/after comparison
      supporting/      provenance (params_used.yaml, manifest copies, roster),
                       coordinates, per_study_data, taxon_impact

Intermediates are removed after a successful --step standardize run unless
keep_intermediates: true; a Phase-1-only run keeps them, because Phase 2 needs them.

## Practical guidance

- **Prefer the pause and the roster to automatic inclusion.** The window-support
  gate and the per-study verdicts make each inclusion decision auditable; automatic
  mode applies the verdicts without review and excludes any study that fails them.
- **Take the scope of legitimate cross-study comparison from the validation report
  generated for each synthesis rather than assuming it.** Restrict comparison to the
  components of community structure the report shows to be preserved, or state the
  expected drift explicitly using the taxon-impact outputs. Abundance-weighted and
  compositional analyses are typically the most robust; richness is systematically
  altered by harmonisation and should be compared only within the harmonised table,
  treating its values as their own baseline.
- **For co-occurrence networks, filter low-prevalence MetaASVs before inference** —
  correlations for features observed in a handful of samples cannot be estimated
  reliably, and because many MetaASVs are recovered from only a subset of
  contributing studies, the pooled table is enriched for features whose absences
  reflect study membership rather than ecology; network inference can misread this
  shared sparsity as co-occurrence structure. Networks are also better supported at
  aggregated ranks: in the paper's validation, hub structure was conserved at genus
  and family but not at MetaASV level. Note that a zero in the merged table means
  "not observed in that study's window", which is not the same as ecological absence.

## Reproducing the SeCAT paper

The complete inputs and outputs of the paper run — study manifests, params_used.yaml,
selection roster, the harmonised MetaASV dataset and all validation outputs — are
archived at Zenodo: https://doi.org/10.5281/zenodo.22113757 (produced with SeCAT
v5.0.1, https://doi.org/10.5281/zenodo.22113432; v5.0.6 differs only in report
labelling, container build pinning and metadata). The deposit's supporting/provenance
directory is the authoritative record of the run's manifest, parameters and roster;
a sanitised copy of that configuration (manifest, roster, effective parameters) is
kept in this repository under `examples/paper_run/`.

## Citation

If you use SeCAT, please cite the software (see CITATION.cff) and:

> Terrey, D., Huber, P., Pereira, E., Ferguson, R. F., Dumbrell, A. J., Sweet, M.,
> Bani, A. & Metz, S. SeCAT (Sequence Consensus Amplicon Trimming): a
> simulation-informed Nextflow pipeline for cross-study harmonisation of 16S rRNA
> metabarcoding datasets. (In review, Molecular Ecology Resources.)

Zenodo DOIs: this release https://doi.org/10.5281/zenodo.22133165; all versions
https://doi.org/10.5281/zenodo.22113431 (minted automatically from GitHub releases).

## Licence

MIT — see LICENSE.
