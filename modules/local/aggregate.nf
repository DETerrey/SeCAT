process AGGREGATE {
    tag "aggregating all results"
    label 'mem_64g'
    publishDir(path: { "${params.outdir}/aggregated_data" }, mode: 'copy', pattern: "aggregated_data/*")
    // Duplicate user-facing verdict + stats tables into final_outputs/3_verdicts/ (flat files at task root)
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verdicts_dir}" }, mode: 'copy', pattern: "master_verdict_table.csv")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verdicts_dir}" }, mode: 'copy', pattern: "simulation_baseline_statistics.csv")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verdicts_dir}" }, mode: 'copy', pattern: "simulation_retention_curves.csv")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verdicts_dir}" }, mode: 'copy', pattern: "excluded_studies.csv")
    // Per-study + combined taxon-impact tables into final_outputs/supporting/taxon_impact/
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_taxon_impact_dir}" }, mode: 'copy', pattern: "master_taxon_impact_summary.csv")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_taxon_impact_dir}" }, mode: 'copy', pattern: "taxon_impact_*.csv")

    input:
    path real_results
    path sim_results
    path study_coords
    path consensus_info

    output:
    path "master_verdict_table.csv",           emit: master_verdict_table
    path "simulation_baseline_statistics.csv", emit: baseline_stats,   optional: true
    path "simulation_retention_curves.csv",    emit: retention_curves, optional: true
    path "excluded_studies.csv",               emit: excluded_studies, optional: true
    path "master_taxon_impact_summary.csv",    emit: taxon_impact_summary, optional: true
    path "taxon_impact_*.csv",                 emit: taxon_impact_studies, optional: true
    path "aggregated_data",                                    emit: aggregated_dir

    script:
    """
    # v5.0.0 detection-threshold wiring fix: AGGREGATE now honours params.yaml
    # (distance_cutoff_threshold, null_model_p_threshold, null_model_min_consecutive).
    # This comment changes the task hash so -resume re-runs AGGREGATE after the R fix.
    # cache-bust: excluded_studies.csv pre-assessment audit table
    mkdir -p output/intermediate output/real_data_results output/simulation_results aggregated_data
    cp ${study_coords}   output/intermediate/study_alignment_coords.csv 2>/dev/null || true
    cp ${consensus_info} output/intermediate/consensusregioninfo.csv    2>/dev/null || true
    for f in *_results.rds; do
        [[ "\$f" == *__seed_* ]] && continue
        study=\$(basename \$f _results.rds)
        mkdir -p output/real_data_results/\${study}
        cp \$f output/real_data_results/\${study}/\${study}_results.rds
    done
    for f in *__seed_*__results.rds; do
        fname=\$(basename \$f)
        task_id=\$(echo \$fname | sed 's/__seed_.*//')
        seed=\$(echo \$fname | sed 's/.*__seed_//' | sed 's/__results.rds//')
        mkdir -p output/simulation_results/\${task_id}/seed_\${seed}
        cp \$f output/simulation_results/\${task_id}/seed_\${seed}/results.rds
    done
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_CHANGEPOINT_METHOD="${params.changepoint_penalty_method}"
    export SECAT_CHANGEPOINT_MULTIPLIER="${params.changepoint_penalty_multiplier}"
    export SECAT_CHANGEPOINT_SEED="${params.changepoint_seed}"
    export SECAT_CHANGEPOINT_SCAN_MIN="${params.changepoint_scan_min}"
    export SECAT_CHANGEPOINT_SCAN_MAX="${params.changepoint_scan_max}"
    export SECAT_CHANGEPOINT_SCAN_BY="${params.changepoint_scan_by}"
    export SECAT_CHANGEPOINT_BOOTSTRAP_N="${params.changepoint_bootstrap_n}"
    export SECAT_NULL_MODEL_P="${params.null_model_p_threshold}"
    export SECAT_NULL_MODEL_MIN_CONSEC="${params.null_model_min_consecutive}"
    export SECAT_DISTANCE_CUTOFF="${params.distance_cutoff_threshold}"
    export SECAT_BC_CEILING="${params.bc_ceiling_threshold}"
    export SECAT_RETENTION_FLOOR="${params.retention_floor_threshold}"
    export SECAT_VERDICT_TOLERANCE_BP="${params.verdict_consensus_tolerance_bp}"
    export SECAT_CONSENSUS_OPT_THRESHOLD="${params.consensus_optimization_threshold}"
    export SECAT_MIN_CONSENSUS_STUDIES="${params.min_consensus_studies}"
    export SECAT_NULL_MODEL_MIN_TRIM_BP="${params.null_model_min_trim_bp ?: 5}"
    export SECAT_DISTANCE_CUTOFF_MIN_BP="${params.distance_cutoff_min_trim_bp ?: 5}"
    export SECAT_MIN_CONSENSUS_COVERAGE="${params.min_consensus_coverage ?: 0.50}"
    export SECAT_EXCLUDE_UNSTABLE="${params.exclude_unstable_from_consensus == null ? true : params.exclude_unstable_from_consensus}"
    export SECAT_OUTDIR="./output"
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/07_aggregate.R
    cp output/aggregated_data/* aggregated_data/ 2>/dev/null || true
    # Stage user-facing tables at task root so publishDir can mirror them to 3_verdicts/
    cp aggregated_data/master_verdict_table.csv           ./ 2>/dev/null || true
    cp aggregated_data/simulation_baseline_statistics.csv ./ 2>/dev/null || true
    cp aggregated_data/simulation_retention_curves.csv    ./ 2>/dev/null || true
    cp aggregated_data/excluded_studies.csv               ./ 2>/dev/null || true
    cp aggregated_data/master_taxon_impact_summary.csv    ./ 2>/dev/null || true
    cp aggregated_data/taxon_impact_*.csv                 ./ 2>/dev/null || true
    """
}
