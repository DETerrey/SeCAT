#!/usr/bin/env nextflow
// ============================================================================
// SeCAT: Sequence Consensus Amplicon Trimming
// Nextflow DSL2 Main Workflow (v4.1)
//
// Converted from SGE/bash orchestration (secat_pipeline.sh V3.4).
// v4.1 changes:
//   - Updated 08_gen_report.R (eff_inner/outer colouring, dg_outside_cons logic)
//   - Added VALIDATE process (validation_taxon_v2.R, multi-tier ecological QC)
//   - Added jasmin SLURM profile for LOTUS2
//
// USAGE:
//   nextflow run main.nf -profile jasmin,singularity -params-file params.yaml
//   nextflow run main.nf -profile slurm,singularity  --manifest secat_manifest.tsv
//   nextflow run main.nf -profile sge,singularity    --manifest secat_manifest.tsv
//   nextflow run main.nf -profile local              --manifest secat_manifest.tsv
//
// RESUME after failure:
//   nextflow run main.nf -profile jasmin,singularity -resume
//
// See nextflow.config and params.yaml for all configurable parameters.
// ============================================================================

nextflow.enable.dsl = 2

include { STUDY_MAPPING       } from './modules/local/study_mapping'
include { COLLECT_MAPS        } from './modules/local/collect_maps'
include { PREPARE_SIMS        } from './modules/local/prepare_sims'
include { RUN_SIMULATION      } from './modules/local/run_simulation'
include { ANALYSE_REAL        } from './modules/local/analyse_real'
include { AGGREGATE           } from './modules/local/aggregate'
include { GENERATE_VERDICTS   } from './modules/local/generate_verdicts'
include { GENERATE_REPORT     } from './modules/local/generate_report'
include { GENERATE_INDEX      } from './modules/local/generate_index'
include { TRIM_SEQUENCES      } from './modules/local/trim_sequences'
include { MERGE_DATASETS      } from './modules/local/merge_datasets'
include { VALIDATE            } from './modules/local/validate'
include { PRIMER_MAPPING      } from './modules/local/primer_mapping'
include { GENERATE_PRIMER_DBS } from './modules/local/primer_mapping'

workflow {

    if (!params.manifest) {
        error "ERROR: --manifest is required. E.g. --manifest secat_manifest.tsv"
    }
    if (!params.reference_db) {
        error "ERROR: --reference_db is required. Provide path to SILVA aligned FASTA."
    }

    log.info """
    ============================================================
      SeCAT v4.1 (Nextflow)
    ============================================================
      Manifest        : ${params.manifest}
      Reference DB    : ${params.reference_db}
      Analysis mode   : ${params.analysis_mode}
      Simulations/set : ${params.num_simulations}
      Trim step mode  : ${params.trim_step_mode}
      Run validation  : ${params.run_validation}
      Executor        : ${workflow.profile}
      Work dir        : ${workflow.workDir}
      Output dir      : ${params.outdir}
    ============================================================
    """.stripIndent()

    manifest_ch = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.study_name, row) }

    indexed_studies_ch = manifest_ch
        .collect()
        .flatMap { rows ->
            rows.withIndex().collect { item, idx ->
                tuple(idx + 1, item[0], item[1])
            }
        }

    if (params.analysis_mode == 'study') {

        STUDY_MAPPING(
            indexed_studies_ch,
            file(params.reference_db)
        )

        all_mapping_parts_ch = STUDY_MAPPING.out.mapping_part.collect()
        COLLECT_MAPS(all_mapping_parts_ch)

        collected_coords = COLLECT_MAPS.out.study_coords
        collected_maps   = COLLECT_MAPS.out.collected_maps

    } else if (params.analysis_mode == 'primer') {

        PRIMER_MAPPING(file(params.reference_db))

        primer_coords_ch = PRIMER_MAPPING.out.primer_coords
            .splitCsv(header: true)
            .map { row -> tuple(row.primer_name, file(params.reference_db)) }

        GENERATE_PRIMER_DBS(primer_coords_ch)

        collected_coords = PRIMER_MAPPING.out.primer_coords
        collected_maps   = PRIMER_MAPPING.out.primer_coords

    } else {
        error "ERROR: analysis_mode must be 'study' or 'primer', got: ${params.analysis_mode}"
    }

    PREPARE_SIMS(
        collected_coords,
        file(params.reference_db)
    )

    ANALYSE_REAL(
        indexed_studies_ch,
        collected_coords,
        PREPARE_SIMS.out.consensus_info
    )

    sim_tasks_ch = PREPARE_SIMS.out.sim_tasks_csv
        .splitCsv(header: true)
        .map { row ->
            tuple(row.task_id, row.simulation_seed, row.num_steps, row.amplicon_length)
        }

    RUN_SIMULATION(
        sim_tasks_ch,
        PREPARE_SIMS.out.sim_reference_subset,
        collected_coords,
        PREPARE_SIMS.out.consensus_info
    )

    all_real_results_ch = ANALYSE_REAL.out.results_rds.collect()
    all_sim_results_ch  = RUN_SIMULATION.out.results_rds.collect()

    AGGREGATE(
        all_real_results_ch,
        all_sim_results_ch,
        collected_coords,
        PREPARE_SIMS.out.consensus_info
    )

    GENERATE_VERDICTS(AGGREGATE.out.master_verdict_table)

    GENERATE_REPORT(
        AGGREGATE.out.aggregated_dir,
        GENERATE_VERDICTS.out.verdict_data
    )

    GENERATE_INDEX(
        AGGREGATE.out.aggregated_dir,
        GENERATE_VERDICTS.out.verdict_data
    )

    if (params.auto_trim) {
        log.warn "auto_trim=true: Trimming with KEEP-only verdicts. " +
                 "Review output/aggregated_data/verdict_data_all_levels.csv first."
        TRIM_SEQUENCES(
            GENERATE_VERDICTS.out.verdict_data,
            COLLECT_MAPS.out.study_coords,
            PREPARE_SIMS.out.consensus_info,
            manifest_ch.collect()
        )
        MERGE_DATASETS(
            TRIM_SEQUENCES.out.standardized_fastas.collect(),
            TRIM_SEQUENCES.out.trim_summary,
            manifest_ch.collect()
        )
        if (params.run_validation) {
            VALIDATE(
                MERGE_DATASETS.out.combined_fasta,
                MERGE_DATASETS.out.combined_otu,
                MERGE_DATASETS.out.combined_tax,
                MERGE_DATASETS.out.combined_meta.ifEmpty(file('NO_META')),
                MERGE_DATASETS.out.asv_mapping.ifEmpty(file('NO_MAPPING')),
                AGGREGATE.out.aggregated_dir
            )
        }
    } else {
        log.info """
        ============================================================
        Pipeline complete. Review verdicts then run:
          Rscript R/11_select_studies.R
          nextflow run main.nf --entry STANDARDIZE -resume
        ============================================================
        """.stripIndent()
    }
}

workflow STANDARDIZE {

    selected_file  = file("${params.outdir}/aggregated_data/selected_studies_for_trim.txt",
                          checkIfExists: true)
    consensus_info = file("${params.outdir}/intermediate/consensusregioninfo.csv",
                          checkIfExists: true)
    study_coords   = file("${params.outdir}/intermediate/study_alignment_coords.csv",
                          checkIfExists: true)
    aggregated_dir = file("${params.outdir}/aggregated_data", checkIfExists: true)

    manifest_ch = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.study_name, row) }

    TRIM_SEQUENCES(
        selected_file,
        study_coords,
        consensus_info,
        manifest_ch.collect()
    )

    MERGE_DATASETS(
        TRIM_SEQUENCES.out.standardized_fastas.collect(),
        TRIM_SEQUENCES.out.trim_summary,
        manifest_ch.collect()
    )

    if (params.run_validation) {
        VALIDATE(
            MERGE_DATASETS.out.combined_fasta,
            MERGE_DATASETS.out.combined_otu,
            MERGE_DATASETS.out.combined_tax,
            MERGE_DATASETS.out.combined_meta.ifEmpty(file('NO_META')),
            MERGE_DATASETS.out.asv_mapping.ifEmpty(file('NO_MAPPING')),
            aggregated_dir
        )
    }
}
