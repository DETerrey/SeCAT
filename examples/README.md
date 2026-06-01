# SeCAT input examples

This folder contains a complete, minimal set of example input files. **None of these point to real data** — paths are placeholders. Copy them to a working location, edit the paths/values, and pass the modified `params_example.yaml` to Nextflow.

The examples are intentionally tiny (2 studies, 3 samples each) so you can read them in a text editor and see every column. Real datasets have hundreds to tens of thousands of ASVs per study; the format is unchanged.

---

## What's here

| File | Role | When you edit it |
|---|---|---|
| `manifest_example.tsv` | The master input — one row per study. Tells SeCAT where every study's data lives. | Once, when you start a new analysis. |
| `params_example.yaml` | All tunable parameters. Sets paths, choice of primer/study mode, simulation depth, statistical thresholds. | Once at start; revisit if exploring parameter sensitivity. |
| `selection_roster_example.txt` | A curated list of studies to include in the trim/merge step. | Between the verdict step and the trim step (after manual review). |
| `study_A/feature_table_example.tsv` | What an ASV/OTU count table must look like (the file `asv_counts_path` should point at). | Never — you produce this from DADA2 / QIIME2 / Deblur. |
| `study_A/taxonomy_example.tsv` | What a taxonomy file must look like (the file `taxonomy_path` should point at). | Never — produced by the same upstream tool. |
| `study_A/sequences_example.fasta` | What an ASV FASTA must look like (the file `asv_fasta_path` should point at). | Never — produced upstream. |
| `study_A/metadata_example.tsv` | What a sample metadata file must look like (the file `metadata_path` should point at). | You curate this yourself per study. |

---

## How the files relate

```
params_example.yaml
   ├── manifest:        ──►  manifest_example.tsv
   │                              ├── row: Study_A ──► asv_fasta_path    ──► study_A/sequences_example.fasta
   │                              │                   ├── asv_counts_path ──► study_A/feature_table_example.tsv
   │                              │                   ├── taxonomy_path   ──► study_A/taxonomy_example.tsv
   │                              │                   └── metadata_path   ──► study_A/metadata_example.tsv
   │                              └── row: Study_B  ──► (same four file types)
   ├── reference_db:    ──►  /path/to/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta
   └── selection_file:  ──►  selection_roster_example.txt
```

---

## Format requirements (read this before you edit)

### Manifest (`manifest_example.tsv`)

Tab-separated. **Headers are case-sensitive and must match exactly.**

| Column | Required | Type | Notes |
|---|---|---|---|
| `study_name` | yes | string | Unique. Use `Author_Year` (e.g., `Beatty_2022`). No spaces, no slashes. |
| `primer_name` | yes | string | E.g., `515F_806R`. Used only as a label in reports. |
| `fwd_seq` | yes | DNA string | Forward primer sequence, IUPAC ambiguity codes allowed. |
| `rev_seq` | yes | DNA string | Reverse primer sequence. |
| `initial_fwd_trim` | yes | int | bp to trim from the 5' end *before* SeCAT does anything. Use 0 if your sequences are already primer-trimmed. |
| `initial_rev_trim` | yes | int | bp to trim from the 3' end before SeCAT processing. Use 0 if already trimmed. |
| `asv_fasta_path` | yes | absolute path | FASTA of ASV/OTU representative sequences. |
| `asv_counts_path` | yes | absolute path | Tab-separated feature table (ASVs × samples). See format below. |
| `taxonomy_path` | yes | absolute path | Tab-separated taxonomy assignments. |
| `metadata_path` | optional | absolute path | Tab-separated sample metadata. Leave blank or write `NA` if none. |
| `metadata_variable` | optional | string | Column name in metadata used to filter samples (e.g., `HolobiontArea`). |
| `metadata_value` | optional | string | Value(s) to keep. Leave blank to keep all samples. |

Paths are **absolute** and must be readable from the compute node (not just the login node). On JASMIN, that means anything under `/home/users/<you>/` or `/work/scratch-nopw2/` or `/work/xfc/` will be fine.

### Feature table (`study_A/feature_table_example.tsv`)

The first column is the ASV/OTU identifier. Each subsequent column is one sample. Values are integer counts (not relative abundances).

The first column header can be `#OTU ID`, `#Feature ID`, `OTU_ID`, `Feature_ID`, `taxonID` — SeCAT's loader handles all common variants. QIIME2 BIOM exports often produce `#OTU ID` with a leading `#` and the file may include a top comment line — that's also handled.

### Taxonomy (`study_A/taxonomy_example.tsv`)

Two-column TSV with header. Column 1 = ASV ID matching the feature table. Column 2 = the taxonomy string. Format is flexible:

- QIIME2 style: `d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; ...`
- SILVA style: `Bacteria;Proteobacteria;Gammaproteobacteria;...`
- With `Confidence` as a third column — also fine.

### Metadata (`study_A/metadata_example.tsv`)

Tab-separated. The first column is the sample identifier and **must match the column headers of the feature table exactly**. Subsequent columns are sample-level variables — anything you want to carry through to downstream analysis. SeCAT will harmonise common column names (`Lat` → `Latitude`, `seagrass_species` → `Host_Species`, etc.) if `harmonize_metadata: true` (the default).

### FASTA (`study_A/sequences_example.fasta`)

Standard FASTA. The header (after `>`) is the ASV ID and **must match the row names in the feature table and the taxonomy table**. Sequences should already be primer-trimmed (or you set `initial_fwd_trim` / `initial_rev_trim` to do it inside SeCAT).

---

## Common pitfalls

1. **Sample IDs don't match between feature table and metadata.** SeCAT will silently drop the mismatched samples. Always check the cleaned manifest after a run (`output/cleaned_data/secat_manifest_clean.tsv`).
2. **ASV IDs don't match between feature table, taxonomy, and FASTA.** Same — silent drop. Use `grep '>' my.fasta | head` and `head -1 my_feature_table.tsv` to spot-check.
3. **Mixing relative paths with absolute paths in the manifest.** Always use absolute paths.
4. **Excel mangling.** Never open a manifest or feature table in Excel — it will reformat sample IDs that look like dates, drop leading zeros, and silently corrupt sequence IDs. Use a plain text editor or LibreOffice with the right import options.
5. **CRLF line endings on Windows.** If you create the manifest on Windows, run `dos2unix` on it before passing to Nextflow.

---

## Verifying your inputs before launching

Run the healthcheck script from the project root to validate that all files referenced in the manifest exist:

```bash
./secat_healthcheck.sh
```

This is fast (no compute) and catches the most common configuration errors before you submit a 24-hour SLURM job.
