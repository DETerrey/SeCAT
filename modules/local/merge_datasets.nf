process MERGE_DATASETS {
    tag "merging datasets"
    label 'mem_30g'
    // ── Existing publishDir layout (preserved for -resume compatibility) ──
    publishDir(path: { "${params.outdir}/meta_analysis" }, mode: 'copy', pattern: "combined_*")
    publishDir(path: { "${params.outdir}/meta_analysis" }, mode: 'copy', pattern: "asv_mapping_final.tsv")
    publishDir(path: { "${params.outdir}/comparison" }, mode: 'copy', pattern: "pre_consensus/*.tsv")
    publishDir(path: { "${params.outdir}/comparison" }, mode: 'copy', pattern: "post_consensus/*.tsv")

    // ── Consolidated final outputs (user-facing) ────────────────────────────
    // These duplicate the published files into ${params.final_outputs}/ so
    // downstream users only need to look in one place.
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_tables_dir}" }, mode: 'copy', pattern: "combined_*")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_tables_dir}" }, mode: 'copy', pattern: "asv_mapping_final.tsv")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_comparison_dir}" }, mode: 'copy', pattern: "pre_consensus/*.tsv")
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_comparison_dir}" }, mode: 'copy', pattern: "post_consensus/*.tsv")

    input:
    path standardized_fastas
    path trim_summary
    path manifest_rows
    path selected_studies_file

    output:
    path "combined_sequences.fasta",        emit: combined_fasta
    path "combined_feature_table.tsv",      emit: combined_otu
    path "combined_taxonomy.tsv",           emit: combined_tax
    path "combined_metadata.tsv",           emit: combined_meta,        optional: true
    path "asv_mapping_final.tsv",           emit: asv_mapping,          optional: true
    path "pre_consensus/feature_table.tsv", emit: pre_otu
    path "pre_consensus/taxonomy.tsv",      emit: pre_tax
    path "pre_consensus/metadata.tsv",      emit: pre_meta,             optional: true
    path "post_consensus/feature_table.tsv",emit: post_otu
    path "post_consensus/taxonomy.tsv",     emit: post_tax
    path "post_consensus/metadata.tsv",     emit: post_meta,            optional: true
    path "post_consensus/sequences.fasta",  emit: post_fasta,           optional: true

    script:
    def absOutdir = params.outdir.startsWith('/') ? params.outdir : "${projectDir}/${params.outdir}"
    """
    mkdir -p standardized_datasets aggregated_data \
             meta_analysis intermediate \
             comparison/pre_consensus comparison/post_consensus \
             pre_consensus post_consensus

    for f in ${standardized_fastas}; do cp \$f standardized_datasets/; done
    cp ${trim_summary}           standardized_datasets/trim_summary.csv
    cp ${selected_studies_file}  aggregated_data/selected_studies_for_trim.txt

    ln -s ${absOutdir}/intermediate/aligned_fastas \
          intermediate/aligned_fastas 2>/dev/null || true

    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_VSEARCH_PATH="vsearch"
    export SECAT_MERGE_METHOD="${params.merge_method}"
    export SECAT_HARMONIZE_METADATA="${params.harmonize_metadata}"
    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    export SECAT_SELECTION_MODE="${params.selection_mode ?: 'file'}"
    export SECAT_SELECTION_FILE="${params.selection_file ?: ''}"

    ln -s ${projectDir}/R R 2>/dev/null || true
    # cache-bust: taxonomy column aliases (Assignation/per_ident)
    Rscript ${projectDir}/R/13_merge_datasets.R

    # Stage post-consensus outputs to work dir root for Nextflow to capture
    cp meta_analysis/combined_sequences.fasta      ./                       2>/dev/null || true
    cp meta_analysis/combined_feature_table.tsv    ./                       2>/dev/null || true
    cp meta_analysis/combined_taxonomy.tsv         ./                       2>/dev/null || true
    cp meta_analysis/combined_metadata.tsv         ./                       2>/dev/null || true
    cp meta_analysis/asv_mapping_final.tsv         ./                       2>/dev/null || true

    cp meta_analysis/combined_feature_table.tsv    post_consensus/feature_table.tsv  2>/dev/null || true
    cp meta_analysis/combined_taxonomy.tsv         post_consensus/taxonomy.tsv       2>/dev/null || true
    cp meta_analysis/combined_metadata.tsv         post_consensus/metadata.tsv       2>/dev/null || true
    cp meta_analysis/combined_sequences.fasta      post_consensus/sequences.fasta    2>/dev/null || true

    # Stage pre-consensus outputs
    cp comparison/pre_consensus/feature_table.tsv  pre_consensus/feature_table.tsv   2>/dev/null || true
    cp comparison/pre_consensus/taxonomy.tsv       pre_consensus/taxonomy.tsv        2>/dev/null || true
    cp comparison/pre_consensus/metadata.tsv       pre_consensus/metadata.tsv        2>/dev/null || true
    """
}
