process PREPARE_SIMS {
    tag "preparing simulation tasks"
    cache 'lenient'
    label 'mem_4g'
    publishDir(path: { "${params.outdir}/intermediate" }, mode: 'copy')
    // Consensus region info goes to final_outputs/coordinates/ for user interpretation
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_coordinates_dir}" }, mode: 'copy', pattern: "consensus_region_info.csv")

    input:
    path study_coords
    val reference_db
    path clean_manifest
    path aligned_fastas

    output:
    path "simulation_tasks.csv",              emit: sim_tasks_csv
    path "simulation_reference_subset.fasta", emit: sim_reference_subset
    path "consensus_region_info.csv",           emit: consensus_info
    path "null_structure.rds",                emit: null_structure, optional: true
    path "study_alignment_coords.csv",        emit: study_align_coords, optional: true

    script:
    """
    export SECAT_REFERENCE_DB="${reference_db}"
    export SECAT_MANIFEST="\$(realpath ${clean_manifest})"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_NUM_SIMULATIONS="${params.num_simulations}"
    export SECAT_SIM_MAX_SILVA_SUBSET="${params.sim_max_silva_subset}"
    export SECAT_SIM_USE_PREBUILT="${params.sim_use_prebuilt_subset}"
    export SECAT_TRIM_STEP_MODE="${params.trim_step_mode}"
    export SECAT_TRIM_INCREMENT="${params.trim_increment}"
    export SECAT_DEFAULT_MAX_TRIM_STEPS="${params.default_max_trim_steps}"
    export SECAT_CONSENSUS_BUFFER_STEPS="${params.consensus_buffer_steps}"
    export SECAT_CONSENSUS_OPT_THRESHOLD="${params.consensus_optimization_threshold}"
    export SECAT_MIN_CONSENSUS_STUDIES="${params.min_consensus_studies}"
    export SECAT_OUTDIR="."
    mkdir -p intermediate
    cp ${study_coords} ./intermediate/study_alignment_coords.csv
    export SECAT_PROJECTDIR="${projectDir}"
    export SECAT_EXCLUDE_UNSTABLE="${params.exclude_unstable_from_consensus}"
    export SECAT_SIM_MAX_TEMPLATE="${params.sim_max_template}"
    export SECAT_SIM_CLUSTER_ID="${params.sim_cluster_identity}"
    # cache-bust: null-structure measurement fix
    # cache-bust: includeTerminalGaps fix
    Rscript ${projectDir}/R/04_prepare_sims.R
    cp ./intermediate/simulation_tasks.csv              ./simulation_tasks.csv
    cp ./intermediate/simulation_reference_subset.fasta ./simulation_reference_subset.fasta
    cp ./intermediate/consensus_region_info.csv           ./consensus_region_info.csv
    cp ./intermediate/study_alignment_coords.csv        ./study_alignment_coords.csv 2>/dev/null || true
    """
}
