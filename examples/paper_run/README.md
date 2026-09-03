# Reproducing the SeCAT paper run

This folder preserves the exact configuration of the demonstration run in the
SeCAT paper (thirteen marine 16S rRNA datasets; SILVA 138.2 SSURef NR99
full-length aligned reference).

- params_paper.yaml  - effective parameters of the paper run (paths are
  placeholders; point manifest/reference_db/outdir/bind_paths at your system)
- secat_manifest.tsv - the thirteen-study input manifest (paths are
  placeholders; the source datasets are listed in Supplementary Table S1 of
  the paper and available from their original archive accessions)
- selection_roster.txt - the eleven studies carried into Phase 2

Run Phase 1, review 1_report/, then rerun with --step standardize -resume.
Use release v5.0.6 and the published container (ghcr.io/deterrey/secat:5.0.6).
