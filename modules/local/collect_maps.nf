process COLLECT_MAPS {
    tag "collecting ${mapping_parts.size()} studies"
    label 'mem_4g'
    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path mapping_parts

    output:
    path "study_alignment_coords.csv",  emit: study_coords
    path "study_mapping_summary.csv",   emit: collected_maps

    script:
    """
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_OUTDIR="."
    Rscript ${projectDir}/R/03_collect_maps.R
    """
}
