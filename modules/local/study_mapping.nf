process STUDY_MAPPING {
    tag "${study_name}"
    label 'mem_64g_cpu4'
    publishDir(path: { "${params.outdir}/intermediate/study_mapping_parts" }, mode: 'copy', pattern: "mapping_part_*.csv")
    publishDir(path: { "${params.outdir}/intermediate/asv_coordinates" }, mode: 'copy', pattern: "*_coords.csv")
    publishDir(path: { "${params.outdir}/intermediate/aligned_fastas" }, mode: 'copy', pattern: "*_aligned.fasta")
    // Aligned FASTAs also go to per_study_data/ so users can re-trim or do phylogenetic work
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_per_study_dir}/aligned_fastas" }, mode: 'copy', pattern: "*_aligned.fasta")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_per_study_dir}/asv_coordinates" }, mode: 'copy', pattern: "*_coords.csv")

    input:
    tuple val(task_id), val(study_name), val(manifest_row)
    val reference_db
    path clean_manifest

    output:
    path "mapping_part_${task_id}.csv",  emit: mapping_part
    path "${study_name}_coords.csv",     emit: asv_coords,    optional: true
    tuple val(study_name), path("${study_name}_aligned.fasta"),  emit: aligned_fasta
    
    script:
    def absOutdir = params.outdir.startsWith('/') ? params.outdir : "${projectDir}/${params.outdir}"
    """
    export SGE_TASK_ID=${task_id}
    export SECAT_REFERENCE_DB="${reference_db}"
    export SECAT_MANIFEST="\$(realpath ${clean_manifest})"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_ALIGNMENT_MODE="${params.reference_alignment_mode}"
    export SECAT_SUBSET_SIZE="${params.reference_subset_size}"
    export SECAT_USE_ALL_ASVS="${params.use_all_asvs}"
    export SECAT_VSEARCH_PATH="vsearch"
    export SECAT_ALIGNMENT_METHOD="${params.study_alignment_method}"
    export SECAT_OUTDIR="${absOutdir}"
    export SECAT_PROJECTDIR="${projectDir}"
    export SECAT_MAPPING_SUPPORT_BAND="${params.mapping_support_band}"
    export SECAT_MAPPING_SUPPORT_WARN="${params.mapping_support_warn}"
    export SECAT_MAPPING_SUPPORT_FAIL="${params.mapping_support_fail}"
    Rscript ${projectDir}/R/02_study_mapping.R

    # Symlink outputs into work directory for Nextflow to stage
    ln -sf ${absOutdir}/intermediate/study_mapping_parts/mapping_part_${task_id}.csv ./mapping_part_${task_id}.csv || true
    ln -sf ${absOutdir}/intermediate/asv_coordinates/${study_name}_coords.csv ./${study_name}_coords.csv || true
    ln -sf ${absOutdir}/intermediate/aligned_fastas/${study_name}_aligned.fasta ./${study_name}_aligned.fasta
    """
}
