process TRIM_SEQUENCES {
    tag "trimming sequences"
    label 'mem_16g'
    publishDir "${params.outdir}/standardized_datasets", mode: 'copy', pattern: "*_standardized.fasta"
    publishDir "${params.outdir}/aggregated_data",       mode: 'copy', pattern: "trim_summary.csv"
    input:
    path selected_file
    path study_coords
    path consensus_info
    path manifest_rows
    output:
    path "*_standardized.fasta", emit: standardized_fastas
    path "trim_summary.csv",     emit: trim_summary
    script:
    """
    mkdir -p aggregated_data intermediate/aligned_fastas standardized_datasets
    cp ${selected_file}   aggregated_data/selected_studies_for_trim.txt
    cp ${study_coords}    intermediate/study_alignment_coords.csv
    cp ${consensus_info}  intermediate/consensusregioninfo.csv
    ln -s ${projectDir}/${params.outdir}/intermediate/aligned_fastas/* intermediate/aligned_fastas/ 2>/dev/null || true
    export SECAT_MANIFEST="${params.manifest}"
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
