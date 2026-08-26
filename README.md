# SeCAT — Sequence Consensus Amplicon Trimming

**Version 5.0.1**

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
window, deduplicated at exact 100% identity (VSEARCH) into MetaASVs with a composite-
score taxonomy, and the harmonised dataset is validated against the pre-harmonisation
data across alpha diversity, ordination geometry, dispersion, rank abundance,
taxonomic composition and co-occurrence network structure.

## Requirements

- Nextflow >= 25.10 (26.04 recommended), DSL2
- Singularity/Apptainer (recommended on HPC) or Docker
- The published container: ghcr.io/deterrey/secat:v5.0.1
  (R 4.3.2, DECIPHER 2.30.0, changepoint 2.3, vegan 2.7-5, microeco 2.3.0,
  SpiecEasi 1.1.2, VSEARCH 2.23.0, and dependencies)
- The SILVA 138.2 SSURef NR99 **full-length aligned** FASTA as the reference

## Quick start

    nextflow pull DETerrey/SeCAT -r v5.0.1
    cp params.yaml my_params.yaml       # set manifest, reference_db, outdir, bind_paths
    nextflow run DETerrey/SeCAT -r v5.0.1 -profile slurm,singularity -params-file my_params.yaml
    # review 1_report/ and edit selection_roster.txt, then Phase 2:
    nextflow run DETerrey/SeCAT -r v5.0.1 -profile slurm,singularity -params-file my_params.yaml \
        --step standardize -resume

Executor profiles: slurm, sge, local. Container profiles: singularity, docker.
**The singularity profile requires `bind_paths`** (comma-separated host paths covering
your data, reference and output locations). `-resume` reuses cached work, and the run
can be stopped and resumed between the two phases.

## Input manifest

A tab-separated file with one row per study and exactly these 10 columns
(see assets/manifest_template.tsv and examples/): study_name, primer_name,
initial_fwd_trim, initial_rev_trim, asv_fasta_path, asv_counts_path,
taxonomy_path, metadata_path, metadata_variable, metadata_value. The last two are an
optional sample filter. Full column reference in examples/README.md.

## Key parameters

All parameters are honoured end-to-end (exported by the modules and read by the R
stages). The shipped params.yaml is a template carrying the defaults used in the
SeCAT paper.

- **Input QC (CLEAN_DATA)**: min_sample_depth (10000), min_asv_prevalence (2),
  min_asv_reads (50); chloroplast/mitochondria removal, a minimum of 3 samples per
  study and a feature-table/metadata overlap check are always applied.
- **Coordinate mapping**: mapping_support_band (250 columns), mapping_support_warn
  (0.70), mapping_support_fail (0.40); UNSTABLE studies are excluded from consensus
  determination (exclude_unstable_from_consensus: true).
- **Simulation null**: one synthetic community per study and replicate, built as a
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

After Phase 1, review the per-study reports in 1_report/ (including
SeCAT_Master_Summary_Report.pdf) and the verdict tables in 3_verdicts/. List the
studies to carry forward in selection_roster.txt and rerun with
`--step standardize -resume`. Automatic selection is available but conservative by
design; the roster gives you the final say.

## Outputs

    outdir/final_outputs/
      1_report/        per-study PDF reports + SeCAT_Master_Summary_Report.pdf
      2_dataset/       harmonised MetaASV feature table, taxonomy, metadata, FASTA,
                       ASV->MetaASV map (phyloseq-importable)
      3_verdicts/      master_verdict_table.csv, verdict_data_all_levels.csv, trim_summary.csv
      4_verification/  validation outputs per level (asv/, genus/, family/) + validation_summary.csv
      5_comparison/    pre-consensus (untrimmed) dataset for before/after comparison
      supporting/      provenance (params_used.yaml, manifest copies), coordinates,
                       per_study_data, taxon_impact

Intermediates are removed after a successful run unless keep_intermediates: true.

## Reproducing the SeCAT paper

The paper run's configuration is preserved under examples/paper_run/
(params_paper.yaml, manifest, roster). Run it against release v5.0.1.

## Citation

If you use SeCAT, please cite the software (see CITATION.cff) and:

> Terrey, D., Huber, P., Pereira, E., Ferguson, R. F., Dumbrell, A. J., Sweet, M.,
> Bani, A. & Metz, S. SeCAT (Sequence Consensus Amplicon Trimming): a
> simulation-informed Nextflow pipeline for cross-study harmonisation of 16S rRNA
> metabarcoding datasets. (In review.)

A Zenodo DOI for each release is minted automatically from GitHub releases.

## Licence

MIT — see LICENSE.
