process GENERATE_VERDICTS {
    tag "generating verdicts"
    cache 'lenient'
    label 'mem_4g'
    publishDir(path: { "${params.outdir}/aggregated_data" }, mode: 'copy')
    // Duplicate final verdicts into final_outputs/verdicts/
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verdicts_dir}" }, mode: 'copy', pattern: "verdict_data_all_levels.csv")

    input:
    path master_verdict_table

    output:
    path "verdict_data_all_levels.csv", emit: verdict_data

    script:
    """
    mkdir -p output/aggregated_data
    cp ${master_verdict_table} output/aggregated_data/master_verdict_table.csv
    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/10_gen_verdicts.R
    cp output/aggregated_data/verdict_data_all_levels.csv ./verdict_data_all_levels.csv
    """
}
