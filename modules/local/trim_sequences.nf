process TRIM_SEQUENCES {
    tag "trimming sequences"
    label 'mem_16g'
    publishDir(path: { "${params.outdir}/standardized_datasets" }, mode: 'copy', pattern: "*_standardized.fasta")
    // Preserve the per-study trimmed FASTAs (standardized_datasets is deleted by cleanup)
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_per_study_dir}/standardized_fastas" }, mode: 'copy', pattern: "*_standardized.fasta")
    publishDir(path: { "${params.outdir}/aggregated_data" }, mode: 'copy', pattern: "trim_summary.csv")
    // Duplicate trim summary into final_outputs/verdicts/
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verdicts_dir}" }, mode: 'copy', pattern: "trim_summary.csv")
    input:
    path selected_file
    path study_coords
    path consensus_info
    path manifest_rows
    output:
    path "*_standardized.fasta", emit: standardized_fastas
    path "trim_summary.csv",     emit: trim_summary
    script:
    def absOutdir = params.outdir.startsWith('/') ? params.outdir : "${projectDir}/${params.outdir}"
    """
    mkdir -p aggregated_data intermediate/aligned_fastas standardized_datasets
    cp ${selected_file}   aggregated_data/selected_studies_for_trim.txt
    cp ${study_coords}    intermediate/study_alignment_coords.csv
    cp ${consensus_info}  intermediate/consensus_region_info.csv
    ln -s ${absOutdir}/intermediate/aligned_fastas/* intermediate/aligned_fastas/ 2>/dev/null || true
    export SECAT_MANIFEST="${absOutdir}/cleaned_data/secat_manifest_clean.tsv"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_MIN_CONSENSUS_COVERAGE="${params.min_consensus_coverage}"
    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    export SECAT_SELECTION_MODE="${params.selection_mode ?: 'file'}"
    export SECAT_SELECTION_FILE="${params.selection_file ?: ''}"
    Rscript ${projectDir}/R/12_trim_sequences.R
    cp standardized_datasets/*_standardized.fasta ./ 2>/dev/null || true
    cp standardized_datasets/trim_summary.csv     ./ 2>/dev/null || true
    """
}
