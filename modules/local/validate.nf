process VALIDATE {
    tag "validation (multi-tier)"
    label 'mem_30g'
    publishDir(path: { "${params.outdir}/validation" }, mode: 'copy')
    // Duplicate verification stats into final_outputs/verification/
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verification_dir}" }, mode: 'copy', pattern: "outputs/**",          saveAs: { f -> f.replaceFirst('^outputs/', '') })
    publishDir(path: { "${params.final_outputs}/${params.final_outputs_verification_dir}" }, mode: 'copy', pattern: "logs/validation.log")

input:
    path post_otu,    stageAs: 'post_feature_table.tsv'
    path post_tax,    stageAs: 'post_taxonomy.tsv'
    path post_meta,   stageAs: 'post_metadata.tsv'
    path post_fasta,  stageAs: 'post_sequences.fasta'
    path pre_otu,     stageAs: 'pre_feature_table.tsv'
    path pre_tax,     stageAs: 'pre_taxonomy.tsv'
    path pre_meta,    stageAs: 'pre_metadata.tsv'
    path asv_mapping
    path aggregated_dir
    path clean_manifest,  stageAs: 'clean_manifest.tsv'
    path roster,          stageAs: 'roster.txt'
    path consensus_info,  stageAs: 'consensus_info.csv'
    path trim_summary_in, stageAs: 'trim_summary_in.csv'
    path sharing_in,      stageAs: 'cross_study_sharing_in.csv'

    output:
    path "outputs/**",          emit: validation_outputs
    path "logs/validation.log", emit: validation_log, optional: true

    script:
    """
    # v5.0.0 network change: VALIDATE Tier 3B now uses SPIEC-EASI (was Spearman).
    # This comment changes the task hash so -resume re-runs VALIDATE after the R change.
    mkdir -p logs outputs pre_consensus post_consensus \
             output/intermediate output/standardized_datasets

    # Stage post-consensus (trimmed MetaASV dataset)
    cp post_feature_table.tsv  post_consensus/feature_table.tsv
    cp post_taxonomy.tsv       post_consensus/taxonomy.tsv
    cp post_metadata.tsv       post_consensus/metadata.tsv   2>/dev/null || true
    cp post_sequences.fasta    post_consensus/sequences.fasta 2>/dev/null || true

    # Stage pre-consensus (original untrimmed dataset)
    cp pre_feature_table.tsv   pre_consensus/feature_table.tsv
    cp pre_taxonomy.tsv        pre_consensus/taxonomy.tsv
    cp pre_metadata.tsv        pre_consensus/metadata.tsv    2>/dev/null || true

    # Stage SeCAT intermediate outputs for coordinate verification (Tier 0C)
    cp -r ${aggregated_dir}/* output/ 2>/dev/null || true
    cp ${asv_mapping} asv_metaasv_map.tsv        2>/dev/null || true
    cp cross_study_sharing_in.csv cross_study_sharing.csv 2>/dev/null || true

    # --- Tier 0 staging (previously missing; Tier 0A-0C were silently skipped) ---
    # 0C inputs: consensus coordinates + trim summary at the paths the R expects
    cp consensus_info.csv  output/intermediate/consensus_region_info.csv   2>/dev/null || true
    cp trim_summary_in.csv output/standardized_datasets/trim_summary.csv 2>/dev/null || true
    # 0A/0B input: pre-consensus (cleaned, untrimmed) sequences, concatenated from
    # the cleaned per-study FASTAs in the clean manifest, restricted to rostered
    # studies where a roster is provided (falls back to all cleaned studies).
    ROSTER=roster.txt; [ -f "\$ROSTER" ] || ROSTER=/dev/null
    awk -F'\\t' 'NR==FNR { if (NF) r[\$1]=1; next }
         FNR==1 { for(i=1;i<=NF;i++){ if(\$i=="asv_fasta_path") c=i; if(\$i=="study_name") s=i }; next }
         c && \$c ~ /^\\// { if (!length(r) || r[\$s]) print \$c }' \
        "\$ROSTER" clean_manifest.tsv \
      | while read -r fa; do [ -f "\$fa" ] && cat "\$fa"; done > pre_consensus/sequences.fasta || true
    [ -s pre_consensus/sequences.fasta ] || rm -f pre_consensus/sequences.fasta

    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    export SECAT_VALIDATION_LEVELS="${params.validation_levels ?: "ASV,Genus,Family"}"

    # cache-bust: sharing baseline read from cross_study_sharing.csv

    Rscript ${projectDir}/R/validation_taxon.R \
        0 \
        FALSE \
        2>&1 | tee logs/validation.log
    """
}
