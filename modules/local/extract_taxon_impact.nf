process EXTRACT_TAXON_IMPACT {
    tag "extracting per-taxon impact tables"
    cache 'lenient'
    label 'mem_8g'
    publishDir "${params.final_outputs}/${params.final_outputs_taxon_impact_dir}", mode: 'copy', pattern: "*.csv"
    publishDir "${params.final_outputs}/${params.final_outputs_taxon_impact_dir}", mode: 'copy', pattern: "README.txt"

    input:
    path real_rds_files

    output:
    path "*.csv",      emit: taxon_csvs, optional: true
    path "README.txt", emit: readme,     optional: true

    script:
    """
    # Stage every per-study .rds into the workdir so the R script picks them up
    for f in ${real_rds_files}; do
        ln -sf \$f .
    done
    Rscript ${projectDir}/R/15_extract_taxon_impact.R
    """
}
