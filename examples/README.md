# SeCAT examples

- manifest_example.tsv — a minimal input manifest
- params_example.yaml — a clean parameter template with the shipped defaults
- paper_run/ — the exact configuration of the SeCAT paper run (params, manifest, roster)

## Manifest column reference

Tab-separated, **exactly 10 columns**, one row per study:

| # | Column | Required | Meaning |
|---|--------|----------|---------|
| 1 | study_name | yes | unique study identifier (used in all outputs) |
| 2 | primer_name | yes | primer-set label, e.g. 515F_806R, 341F_805R |
| 3 | initial_fwd_trim | yes (0 if none) | fixed 5' trim already applied upstream |
| 4 | initial_rev_trim | yes (0 if none) | fixed 3' trim already applied upstream |
| 5 | asv_fasta_path | yes | ASV representative sequences (FASTA) |
| 6 | asv_counts_path | yes | feature table, samples x ASVs |
| 7 | taxonomy_path | yes | per-ASV taxonomy table |
| 8 | metadata_path | yes | sample metadata table |
| 9 | metadata_variable | optional | metadata column to filter samples on |
| 10 | metadata_value | optional | value of that column to retain |

Earlier versions documented fwd_seq/rev_seq primer-sequence columns; the primer mode
is retired and those columns are no longer read.

The statistical defaults in both parameter files are those used in the SeCAT paper:
null-model p < 0.05 over 2 consecutive steps, empirical cutoff 0.20, changepoint
penalty by bootstrap CV scan 0.4–2.4, consensus tolerance 500 columns,
structure-matched simulation null.
