# ==============================================================================
# SCRIPT: scripts/00_clean_data.R (v2.2 — ROBUST NUMERIC LOADING)
# PURPOSE:
#   1. For each study in secat_manifest.tsv:
#      - Filter metadata to keep only target environments (OPTIONAL)
#      - Remove Chloroplast / Mitochondria ASVs
#      - Remove empty samples (pre-existing or post-filtering)
#      - Sync feature table and FASTA to those filters
#   2. Write all cleaned files into a per-study "clean" subdirectory.
#   3. Write an updated manifest: secat_manifest_clean.tsv
#
# v2.2 changes vs v2.1:
#   - load_feature_table() separated from load_table_robust() — feature tables
#     require guaranteed numeric count columns; metadata/taxonomy do not.
#     Previously, colClasses="character" + a heuristic coercion loop caused
#     float-encoded counts ("3653.0") to silently remain as character in some
#     studies (notably Kardish_2023), producing all-zero colSums and dropping
#     all samples.
#   - Full per-sample read count diagnostics printed before and after filtering.
#   - QIIME2 skip logic made explicit and validated.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(data.table)
  library(fs)
})

MANIFEST_IN  <- "secat_manifest.tsv"
MANIFEST_OUT <- "secat_manifest_clean.tsv"

TAXA_TO_REMOVE       <- "Chloroplast|Mitochondria"
MIN_SAMPLES          <- 3
MIN_READS_PER_SAMPLE <- 1

message("========================================================")
message("          SeCAT DATA CLEANING UTILITY v2.2")
message("        (Robust numeric loading + diagnostics)")
message("========================================================")
message(paste("Input manifest :", MANIFEST_IN))
message(paste("Output manifest:", MANIFEST_OUT))
message(paste("Drop taxa matching:", TAXA_TO_REMOVE))
message(paste("Min samples per study:", MIN_SAMPLES))
message(paste("Min reads per sample:", MIN_READS_PER_SAMPLE))
message("========================================================")

if (!file.exists(MANIFEST_IN)) stop("FATAL: Manifest not found: ", MANIFEST_IN)

manifest       <- read_tsv(MANIFEST_IN, show_col_types = FALSE)
manifest_clean <- manifest

# ==============================================================================
# HELPER: LOAD METADATA / TAXONOMY (character columns, safe for mixed types)
# ==============================================================================

load_table_robust <- function(filepath, first_col_name = NULL) {

  lines <- readLines(filepath, n = 5)

  # Skip only pure QIIME2 comment lines (e.g. "# Constructed from biom file").
  # The actual column header in QIIME2 files starts with "#OTU ID" or
  # "#Feature ID" — it contains tabs and must NOT be skipped; fread reads it
  # as the column header row. Pure comment lines have no tabs.
  skip <- 0
  qiime_headers <- NULL
  for (i in seq_along(lines)) {
    line <- lines[i]
    is_pure_comment <- startsWith(line, "#") && !grepl("\t", line)
    if (is_pure_comment) {
      skip <- i
      qiime_headers <- c(qiime_headers, line)
    } else {
      break
    }
  }

  dt <- fread(
    filepath,
    skip        = skip,
    check.names = FALSE,
    data.table  = FALSE,
    sep         = "\t",
    fill        = TRUE,
    quote       = "",
    strip.white = TRUE,
    blank.lines.skip = TRUE,
    colClasses  = "character"
  )

  colnames(dt) <- gsub("^#\\s*", "", colnames(dt))
  colnames(dt) <- trimws(colnames(dt))

  if (!is.null(first_col_name) && ncol(dt) >= 1) {
    colnames(dt)[1] <- first_col_name
  }

  # Remove fully empty columns (character-safe)
  empty_cols <- sapply(dt, function(x) all(is.na(x) | trimws(x) == ""))
  if (any(empty_cols)) dt <- dt[, !empty_cols, drop = FALSE]

  # Remove rows with empty/NA first column
  if (nrow(dt) > 0) {
    first <- trimws(dt[[1]])
    dt <- dt[!is.na(first) & first != "", ]
  }

  result <- as_tibble(dt)
  if (!is.null(qiime_headers)) attr(result, "qiime_headers") <- qiime_headers
  return(result)
}

# ==============================================================================
# HELPER: LOAD FEATURE TABLE (guarantees numeric count columns)
#
# Separated from load_table_robust because feature tables have a known
# structure: col 1 = ASV ID (character), cols 2..N = numeric counts.
# Using colClasses="character" then coercing is safer than relying on fread's
# type inference (which crashes on ambiguous date-like strings elsewhere).
# ==============================================================================

load_feature_table <- function(filepath) {

  lines <- readLines(filepath, n = 10)

  # Skip only pure QIIME2 comment lines (no tabs — not the header row).
  # "#OTU ID\tSRR..." is the actual header and must not be skipped.
  skip <- 0
  qiime_headers <- NULL
  for (i in seq_along(lines)) {
    line <- lines[i]
    is_pure_comment <- startsWith(line, "#") && !grepl("\t", line)
    if (is_pure_comment) {
      skip <- i
      qiime_headers <- c(qiime_headers, line)
    } else {
      break
    }
  }

  message(sprintf("       [FT loader] File: %s", basename(filepath)))
  message(sprintf("       [FT loader] Skipping %d header/comment line(s)", skip))
  if (skip > 0) message(sprintf("       [FT loader] Comment: %s", lines[1]))

  # Read with all-character to avoid any type inference crashes
  dt <- fread(
    filepath,
    skip        = skip,
    check.names = FALSE,
    data.table  = FALSE,
    sep         = "\t",
    fill        = TRUE,
    quote       = "",
    strip.white = TRUE,
    blank.lines.skip = TRUE,
    colClasses  = "character"
  )

  colnames(dt) <- gsub("^#\\s*", "", colnames(dt))
  colnames(dt) <- trimws(colnames(dt))

  # Validate structure: need at least 2 columns (ASV ID + ≥1 sample)
  if (ncol(dt) < 2) stop("Feature table has fewer than 2 columns after loading: ", filepath)

  colnames(dt)[1] <- "ASV_ID"

  # Remove fully empty rows/cols
  empty_cols <- sapply(dt, function(x) all(is.na(x) | trimws(x) == ""))
  if (any(empty_cols)) dt <- dt[, !empty_cols, drop = FALSE]
  if (nrow(dt) > 0) {
    first <- trimws(dt[[1]])
    dt <- dt[!is.na(first) & first != "", ]
  }

  # === CRITICAL: coerce count columns to numeric ===
  # All columns except ASV_ID must be numeric.
  # Counts may be stored as floats ("3653.0") from biom export.
  sample_cols <- names(dt)[-1]
  n_coerced   <- 0
  n_failed    <- 0

  for (col in sample_cols) {
    vals      <- dt[[col]]
    converted <- suppressWarnings(as.numeric(vals))
    na_before <- sum(is.na(vals))
    na_after  <- sum(is.na(converted))
    new_nas   <- na_after - na_before

    if (new_nas > length(vals) * 0.1) {
      # More than 10% new NAs introduced — coercion is mangling real data
      warning(sprintf(
        "[FT loader] Column '%s': %d new NAs after as.numeric() — keeping as character",
        col, new_nas
      ))
      n_failed <- n_failed + 1
    } else {
      dt[[col]] <- converted
      n_coerced <- n_coerced + 1
    }
  }

  message(sprintf("       [FT loader] Coerced %d/%d sample columns to numeric (%d failed)",
                  n_coerced, length(sample_cols), n_failed))

  # Diagnostic: show range of values to confirm counts are real
  all_counts <- unlist(dt[, -1], use.names = FALSE)
  all_counts <- suppressWarnings(as.numeric(all_counts))
  all_counts <- all_counts[!is.na(all_counts)]
  message(sprintf("       [FT loader] Count range: min=%.1f, max=%.1f, mean=%.1f",
                  min(all_counts), max(all_counts), mean(all_counts)))
  message(sprintf("       [FT loader] Zero entries: %d / %d (%.1f%%)",
                  sum(all_counts == 0), length(all_counts),
                  100 * sum(all_counts == 0) / length(all_counts)))

  result <- as_tibble(dt)
  if (!is.null(qiime_headers)) attr(result, "qiime_headers") <- qiime_headers
  return(result)
}

# ==============================================================================
# HELPER: MAKE CLEAN PATH
# ==============================================================================

make_clean_path <- function(path) {
  path       <- path[1]
  parent_dir <- fs::path_dir(path)
  fname      <- fs::path_file(path)
  fs::path(parent_dir, "clean", fname)
}

# ==============================================================================
# HELPER: AUTO-DETECT AND APPLY SAMPLE ID MAPPING (SAMD <-> Run accession)
# ==============================================================================

auto_detect_and_map_samples <- function(ft_sample_cols, meta_clean, sample_id_col,
                                        orig_metadata, study_name) {

  valid_samples <- meta_clean[[sample_id_col]]
  n_matched     <- sum(ft_sample_cols %in% valid_samples)
  match_pct     <- 100 * n_matched / length(valid_samples)

  message(sprintf("     → Sample matching: %d/%d metadata samples found in FT (%.1f%%)",
                  n_matched, length(valid_samples), match_pct))

  if (match_pct >= 90) {
    return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                mapping_applied = FALSE))
  }

  meta_has_samd       <- any(grepl("^SAMD", valid_samples))
  ft_has_run_accessions <- any(grepl("^(DRR|SRR|ERR)", ft_sample_cols))

  if (!meta_has_samd || !ft_has_run_accessions) {
    message("     ⚠️ Low match but no SAMD↔Run mismatch detected")
    return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                mapping_applied = FALSE))
  }

  message("     → Auto-detected SAMD↔Run ID mismatch. Attempting mapping...")

  orig_meta <- load_table_robust(orig_metadata, "SampleID")
  message(sprintf("       Metadata dimensions: %d rows × %d columns",
                  nrow(orig_meta), ncol(orig_meta)))

  # Strategy 1: column by name
  drr_col <- grep("run|accession|srr|drr|err|biosample.*run|experiment|library.*name",
                  names(orig_meta), value = TRUE, ignore.case = TRUE)[1]

  # Strategy 2: column by content
  if (is.na(drr_col)) {
    message("       → Searching columns for run accession codes...")
    for (col_name in names(orig_meta)[-1]) {
      values     <- as.character(orig_meta[[col_name]])
      n_run_codes <- sum(grepl("^(DRR|SRR|ERR)", values, ignore.case = TRUE), na.rm = TRUE)
      if (n_run_codes > length(valid_samples) * 0.3) {
        drr_col <- col_name
        message(sprintf("       ✓ Found run codes in '%s' (%d matches)", col_name, n_run_codes))
        break
      }
    }
  } else {
    message(sprintf("       ✓ Using column '%s' for mapping", drr_col))
  }

  # Strategy 3: sequential (last resort)
  if (is.na(drr_col)) {
    message("       ⚠️ No mapping column found. Attempting sequential mapping...")
    message("          WARNING: This assumes samples are in the same order!")

    if (nrow(meta_clean) <= length(ft_sample_cols)) {
      mapped_samples <- ft_sample_cols[1:nrow(meta_clean)]
      mapping_df     <- tibble(SAMD = meta_clean[[sample_id_col]], RUN = mapped_samples)
      message(sprintf("       ✓ Mapped %d samples sequentially", nrow(mapping_df)))
      message(sprintf("       Example: %s → %s", mapping_df$SAMD[1], mapping_df$RUN[1]))

      valid_samples <- mapping_df$RUN
      meta_clean <- meta_clean %>%
        mutate(original_id = .data[[sample_id_col]]) %>%
        left_join(mapping_df, by = c("original_id" = "SAMD")) %>%
        mutate(!!sample_id_col := coalesce(RUN, .data[[sample_id_col]])) %>%
        select(-RUN, -original_id)

      n_matched <- sum(ft_sample_cols %in% valid_samples)
      match_pct <- 100 * n_matched / length(valid_samples)
      message(sprintf("       → After sequential mapping: %d/%d matched (%.1f%%)",
                      n_matched, length(valid_samples), match_pct))
      return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                  mapping_applied = TRUE))
    } else {
      message("       ❌ Sequential mapping not possible (more meta samples than FT samples)")
      return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                  mapping_applied = FALSE))
    }
  }

  # Apply named-column mapping
  mapping_df <- orig_meta %>%
    select(SAMD = 1, RUN = all_of(drr_col)) %>%
    filter(!is.na(RUN), RUN != "", SAMD %in% valid_samples) %>%
    mutate(RUN = trimws(RUN))

  if (nrow(mapping_df) == 0) {
    message("       ❌ Mapping column found but no valid mappings extracted")
    return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                mapping_applied = FALSE))
  }

  message(sprintf("       ✓ Created mapping for %d samples", nrow(mapping_df)))
  valid_samples <- mapping_df$RUN
  meta_clean <- meta_clean %>%
    left_join(mapping_df, by = setNames("SAMD", sample_id_col)) %>%
    mutate(!!sample_id_col := coalesce(RUN, .data[[sample_id_col]])) %>%
    select(-RUN)

  n_matched <- sum(ft_sample_cols %in% valid_samples)
  match_pct <- 100 * n_matched / length(valid_samples)
  message(sprintf("       → After column mapping: %d/%d matched (%.1f%%)",
                  n_matched, length(valid_samples), match_pct))

  return(list(valid_samples = valid_samples, meta_clean = meta_clean,
              mapping_applied = TRUE))
}

# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================

for (i in seq_len(nrow(manifest))) {
  study      <- manifest[i, ]
  study_name <- study$study_name

  message(paste0("\n>>> Processing Study: ", study_name))

  # ── Metadata filtering parameters ─────────────────────────────────────────
  has_metadata_var   <- "metadata_variable" %in% names(study)
  has_metadata_val   <- "metadata_value"    %in% names(study)
  should_filter_metadata <- FALSE
  TARGET_ENVIRONMENTS    <- NULL
  target_col             <- NULL

  if (has_metadata_var && has_metadata_val) {
    metadata_var_raw <- study$metadata_variable
    metadata_val_raw <- study$metadata_value
    var_is_valid <- !is.na(metadata_var_raw) && nchar(trimws(metadata_var_raw)) > 0
    val_is_valid <- !is.na(metadata_val_raw) && nchar(trimws(metadata_val_raw)) > 0
    if (var_is_valid && val_is_valid) {
      should_filter_metadata <- TRUE
      target_col          <- trimws(metadata_var_raw)
      TARGET_ENVIRONMENTS <- trimws(unlist(strsplit(as.character(metadata_val_raw), ",")))
      message(paste("    Metadata filtering: ENABLED"))
      message(paste("    Target variable:", target_col))
      message(paste("    Target values:", paste(TARGET_ENVIRONMENTS, collapse = ", ")))
    } else {
      message("    Metadata filtering: DISABLED (no variable/value specified)")
    }
  } else {
    message("    Metadata filtering: DISABLED (columns not in manifest)")
  }

  # ── File existence check ───────────────────────────────────────────────────
  orig_counts   <- study$asv_counts_path
  orig_taxonomy <- study$taxonomy_path
  orig_metadata <- study$metadata_path
  orig_fasta    <- study$asv_fasta_path

  paths_orig <- c(counts = orig_counts, taxonomy = orig_taxonomy,
                  metadata = orig_metadata, fasta = orig_fasta)
  missing <- paths_orig[!file.exists(paths_orig)]

  if (length(missing) > 0) {
    message("  ⚠️ Skipping study due to missing files:")
    for (nm in names(missing)) message("     - ", nm, ": ", missing[[nm]])
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  clean_counts   <- make_clean_path(orig_counts)
  clean_taxonomy <- make_clean_path(orig_taxonomy)
  clean_metadata <- make_clean_path(orig_metadata)
  clean_fasta    <- make_clean_path(orig_fasta)

  for (p in unique(c(fs::path_dir(clean_counts), fs::path_dir(clean_taxonomy),
                     fs::path_dir(clean_metadata), fs::path_dir(clean_fasta)))) {
    fs::dir_create(p)
  }

  # ── [1/4] Metadata ────────────────────────────────────────────────────────
  message("  [1/4] Processing Metadata...")
  meta_df       <- load_table_robust(orig_metadata, "SampleID")
  sample_id_col <- "SampleID"

  if (should_filter_metadata && !is.null(target_col) && target_col %in% names(meta_df)) {
    meta_clean <- meta_df %>% filter(.data[[target_col]] %in% TARGET_ENVIRONMENTS)
    message(sprintf("     - Filtered by %s: kept %d / %d samples",
                    target_col, nrow(meta_clean), nrow(meta_df)))
  } else {
    if (should_filter_metadata && !is.null(target_col) && !target_col %in% names(meta_df)) {
      message(sprintf("  ⚠️ Column '%s' not found in metadata — keeping all samples", target_col))
    }
    meta_clean <- meta_df
    message(paste("     - No filtering applied. Kept all", nrow(meta_clean), "samples"))
  }

  if (nrow(meta_clean) < MIN_SAMPLES) {
    message(sprintf("  ⚠️ Only %d samples (minimum %d). Skipping.", nrow(meta_clean), MIN_SAMPLES))
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  # ── [2/4] Taxonomy ────────────────────────────────────────────────────────
  message("  [2/4] Filtering Taxonomy...")
  tax_df    <- load_table_robust(orig_taxonomy, "ASV_ID")
  tax_col   <- names(tax_df)[2]
  tax_clean <- tax_df %>%
    filter(!grepl(TAXA_TO_REMOVE, .data[[tax_col]], ignore.case = TRUE))
  message(sprintf("     - Removed %d ASVs matching Chloroplast/Mitochondria",
                  nrow(tax_df) - nrow(tax_clean)))
  write_tsv(tax_clean, clean_taxonomy)
  valid_asvs <- tax_clean[["ASV_ID"]]

  # ── [3/4] Feature Table ───────────────────────────────────────────────────
  message("  [3/4] Cleaning Feature Table...")

  ft_df         <- load_feature_table(orig_counts)
  qiime_headers <- attr(ft_df, "qiime_headers")
  ft_sample_cols <- names(ft_df)[-1]

  message(sprintf("     → Feature table loaded: %d ASVs × %d samples",
                  nrow(ft_df), length(ft_sample_cols)))

  # Per-sample read totals BEFORE any filtering
  sample_totals_orig <- colSums(ft_df[, -1], na.rm = TRUE)
  zero_orig          <- sum(sample_totals_orig == 0)
  message(sprintf("     → Pre-filter read totals: min=%.0f, median=%.0f, max=%.0f",
                  min(sample_totals_orig), median(sample_totals_orig), max(sample_totals_orig)))
  if (zero_orig > 0) {
    message(sprintf("     ⚠️ %d/%d samples (%.1f%%) have ZERO reads before filtering",
                    zero_orig, length(sample_totals_orig),
                    100 * zero_orig / length(sample_totals_orig)))
  }

  # Sample ID mapping
  mapping_result <- auto_detect_and_map_samples(
    ft_sample_cols, meta_clean, sample_id_col, orig_metadata, study_name)
  valid_samples <- mapping_result$valid_samples
  meta_clean    <- mapping_result$meta_clean

  n_matched <- sum(ft_sample_cols %in% valid_samples)
  match_pct <- 100 * n_matched / length(valid_samples)

  if (match_pct < 50) {
    message(sprintf("     ❌ Only %.1f%% of metadata samples matched in FT. Skipping.", match_pct))
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  # Filter to matched samples × valid ASVs
  cols_to_keep <- c("ASV_ID", intersect(ft_sample_cols, valid_samples))
  if (length(cols_to_keep) == 1) {
    message("     ❌ FATAL: No matching sample columns after mapping. Skipping.")
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  ft_clean <- ft_df %>%
    select(all_of(cols_to_keep)) %>%
    filter(.data[["ASV_ID"]] %in% valid_asvs)

  # Remove ASVs with zero reads across retained samples
  asv_totals <- rowSums(ft_clean[, -1], na.rm = TRUE)
  ft_clean   <- ft_clean[asv_totals > 0, ]
  message(sprintf("     → After ASV zero-filter: %d ASVs retained", nrow(ft_clean)))

  # Per-sample totals AFTER ASV filter — before sample removal
  sample_totals_post <- colSums(ft_clean[, -1], na.rm = TRUE)
  message(sprintf("     → Post-ASV-filter sample totals: min=%.0f, median=%.0f, max=%.0f",
                  min(sample_totals_post), median(sample_totals_post), max(sample_totals_post)))

  # Remove samples below read threshold
  empty_samples <- names(sample_totals_post)[sample_totals_post < MIN_READS_PER_SAMPLE]
  if (length(empty_samples) > 0) {
    message(sprintf("     ⚠️ Removing %d samples with < %d reads",
                    length(empty_samples), MIN_READS_PER_SAMPLE))
    message(sprintf("        Examples: %s", paste(head(empty_samples, 5), collapse = ", ")))
    # Show read totals for removed samples (diagnostic)
    removed_totals <- sort(sample_totals_post[empty_samples])
    message(sprintf("        Removed sample totals: %s",
                    paste(head(removed_totals, 10), collapse = ", ")))
    ft_clean <- ft_clean %>%
      select(all_of(c("ASV_ID", setdiff(names(ft_clean)[-1], empty_samples))))
  }

  message(sprintf("     → Final Feature Table: %d ASVs × %d samples",
                  nrow(ft_clean), ncol(ft_clean) - 1))

  final_valid_asvs <- ft_clean[["ASV_ID"]]
  final_ft_samples <- names(ft_clean)[-1]

  meta_clean <- meta_clean %>% filter(.data[[sample_id_col]] %in% final_ft_samples)
  message(sprintf("     → Aligned metadata: %d samples", nrow(meta_clean)))

  # Write metadata
  write_tsv(meta_clean, clean_metadata)

  # Write feature table (preserve QIIME2 header if present)
  if (!is.null(qiime_headers)) {
    writeLines(c(qiime_headers[1],
                 paste0("#", paste(colnames(ft_clean), collapse = "\t"))),
               clean_counts)
    write_tsv(ft_clean, clean_counts, append = TRUE, col_names = FALSE)
  } else {
    write_tsv(ft_clean, clean_counts, col_names = TRUE)
  }

  # ── [4/4] FASTA ───────────────────────────────────────────────────────────
  message("  [4/4] Cleaning FASTA...")
  seqs       <- readDNAStringSet(orig_fasta)
  seqs_clean <- seqs[names(seqs) %in% final_valid_asvs]
  message(sprintf("     - Removed %d sequences (%d retained)",
                  length(seqs) - length(seqs_clean), length(seqs_clean)))
  writeXStringSet(seqs_clean, clean_fasta)

  # ── Validation ────────────────────────────────────────────────────────────
  if (length(seqs_clean) == 0 || nrow(ft_clean) == 0 || ncol(ft_clean) <= 1) {
    message("  ⚠️ No valid ASVs or samples remaining. Excluding from manifest.")
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  manifest_clean$asv_counts_path[i] <- clean_counts
  manifest_clean$taxonomy_path[i]   <- clean_taxonomy
  manifest_clean$metadata_path[i]   <- clean_metadata
  manifest_clean$asv_fasta_path[i]  <- clean_fasta

  message("  ✅ Study cleaned and validated")
}

# ==============================================================================
# POST-CLEANING VALIDATION
# ==============================================================================

message("\n========================================================")
message("POST-CLEANING VALIDATION")
message("========================================================")

manifest_clean_final <- manifest_clean %>%
  filter(!is.na(asv_counts_path) & !is.na(asv_fasta_path))

validation_issues <- 0

for (i in seq_len(nrow(manifest_clean_final))) {
  study      <- manifest_clean_final[i, ]
  study_name <- study$study_name

  clean_ft   <- load_feature_table(study$asv_counts_path)
  clean_meta <- load_table_robust(study$metadata_path, "SampleID")

  ft_samples   <- names(clean_ft)[-1]
  meta_samples <- clean_meta$SampleID
  n_match      <- sum(ft_samples %in% meta_samples)
  match_pct    <- 100 * n_match / length(ft_samples)

  if (match_pct < 100) {
    message(sprintf("  %s: ⚠️ %d/%d FT samples in metadata (%.1f%%)",
                    study_name, n_match, length(ft_samples), match_pct))
    missing_from_meta <- setdiff(ft_samples, meta_samples)
    extra_in_meta     <- setdiff(meta_samples, ft_samples)
    if (length(missing_from_meta) > 0)
      message(sprintf("    Missing from metadata: %s", paste(head(missing_from_meta, 5), collapse = ", ")))
    if (length(extra_in_meta) > 0)
      message(sprintf("    Extra in metadata: %s",    paste(head(extra_in_meta, 5), collapse = ", ")))
    validation_issues <- validation_issues + 1
  } else {
    message(sprintf("  %s: ✅ Perfect match (%d samples)", study_name, n_match))
  }
}

write_tsv(manifest_clean_final, MANIFEST_OUT)

message("\n========================================================")
message(paste("Cleaning complete. Updated manifest:", MANIFEST_OUT))
message(paste("Studies in final manifest:", nrow(manifest_clean_final)))
if (validation_issues > 0) {
  message(sprintf("⚠️ %d studies have validation issues — review above", validation_issues))
} else {
  message("✅ All studies validated successfully!")
}
message("Use secat_manifest_clean.tsv for downstream MESAP runs.")
message("========================================================")
