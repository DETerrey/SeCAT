process RUN_SIMULATION {
    errorStrategy { task.exitStatus in 137..140 || task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2
    tag "${task_id}__seed_${seed}"
    cache 'lenient'
    label 'mem_18g'
    publishDir(path: { "${params.outdir}/simulation_results/${task_id}/seed_${seed}" }, mode: 'copy', pattern: "results.rds")

    input:
    tuple val(task_id), val(seed), val(num_steps), val(amplicon_length)
    path sim_reference_subset
    path study_coords
    path consensus_info
    path null_structure

    output:
    path "${task_id}__seed_${seed}__results.rds", emit: results_rds

    script:
    """
    export SECAT_REFERENCE_DB="${sim_reference_subset}"
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_VSEARCH_PATH="vsearch"  # binary is on PATH inside container
    export SECAT_TRIM_STEP_MODE="${params.trim_step_mode}"
    export SECAT_TRIM_INCREMENT="${params.trim_increment}"
    export SECAT_DEFAULT_MAX_TRIM_STEPS="${params.default_max_trim_steps}"
    export SECAT_CONSENSUS_BUFFER_STEPS="${params.consensus_buffer_steps}"
    export SECAT_SIM_ABUNDANCE_MODEL="${params.sim_abundance_model}"
    export SECAT_SIM_LOGNORMAL_MU="${params.sim_lognormal_mu}"
    export SECAT_SIM_LOGNORMAL_SIGMA="${params.sim_lognormal_sigma}"
    export SECAT_SIM_ADD_PCR_BIAS="${params.sim_add_pcr_bias}"
    export SECAT_SIM_PCR_CYCLES="${params.sim_pcr_cycles}"
    export SECAT_SIM_PCR_GC_BIAS="${params.sim_pcr_gc_bias}"
    export SECAT_SIM_PCR_OPTIMAL_GC="${params.sim_pcr_optimal_gc}"
    export SECAT_SIM_ADD_ERRORS="${params.sim_add_errors}"
    export SECAT_SIM_ERROR_RATE="${params.sim_error_rate}"
    export SECAT_SIM_INDEL_RATE="${params.sim_indel_rate}"
    export SECAT_SIM_ERROR_POSITION_BIAS="${params.sim_error_position_bias}"
    export SECAT_SIM_ADD_CHIMERAS="${params.sim_add_chimeras}"
    export SECAT_SIM_CHIMERA_RATE="${params.sim_chimera_rate}"
    export SECAT_NULL_STRUCTURE="${null_structure}"
    export SECAT_SIM_MATCHED_NULL="${params.sim_structure_matched_null}"
    export SECAT_SIM_CLUSTER_ID="${params.sim_cluster_identity}"
    export SECAT_SIM_COMMUNITY_SIZE="${params.simulation_community_size}"
    export SECAT_OUTDIR="."
    mkdir -p intermediate
    cp ${sim_reference_subset} intermediate/simulation_reference_subset.fasta 2>/dev/null || true
    cp ${study_coords}   intermediate/study_alignment_coords.csv 2>/dev/null || true
    cp ${consensus_info} intermediate/consensusregioninfo.csv    2>/dev/null || true
    echo "task_id,num_steps,amplicon_length,simulation_seed" > intermediate/simulation_tasks.csv
    echo "${task_id},${num_steps},${amplicon_length},${seed}" >> intermediate/simulation_tasks.csv
    export SECAT_PROJECTDIR="${projectDir}"
    # cache-bust: null taxonomy inheritance fix
    # cache-bust: null taxonomy inheritance fix
    Rscript ${projectDir}/R/05_sim_worker.R "${task_id}" "${seed}"
    mv simulation_results/${task_id}/seed_${seed}/results.rds ./${task_id}__seed_${seed}__results.rds 2>/dev/null || true
    """
}
