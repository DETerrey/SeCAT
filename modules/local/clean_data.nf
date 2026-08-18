process CLEAN_DATA {
    tag "cleaning all studies"
    label 'mem_16g'
    publishDir(path: { "${params.outdir}/cleaned_data" }, mode: 'copy')

    input:
    path manifest

    output:
    path "secat_manifest_clean.tsv", emit: clean_manifest
    path "cleaned/**",               emit: cleaned_files, optional: true

    script:
    def absOutdir = params.outdir.startsWith('/') ? params.outdir : "${projectDir}/${params.outdir}"
    """
    export SECAT_MANIFEST="${manifest}"
    export SECAT_OUTDIR="${absOutdir}"
    export SECAT_PROJECTDIR="${projectDir}"
    # cache-bust: separator-agnostic sample ID matching
    export SECAT_MIN_SAMPLE_DEPTH="${params.min_sample_depth}"
    export SECAT_MIN_ASV_PREVALENCE="${params.min_asv_prevalence}"
    export SECAT_MIN_ASV_READS="${params.min_asv_reads}"
    # cache-bust: integrated depth+abundance QC
    Rscript ${projectDir}/R/00_clean_data.R

    # Symlink cleaned manifest into work directory
    ln -sf ${absOutdir}/cleaned_data/secat_manifest_clean.tsv ./secat_manifest_clean.tsv || true
    ln -sfn ${absOutdir}/cleaned_data/cleaned ./cleaned || true
    """
}

