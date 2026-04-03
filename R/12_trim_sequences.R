#!/usr/bin/env Rscript
# scripts/12_trim_sequences.R
# PURPOSE: Extracts consensus region from ALIGNED sequences

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(Biostrings)
})

# --- LOAD CONFIG ---
cat("\n=== STARTING STANDARDIZATION (From Aligned Sequences) ===\n")
cat("Loading configuration...\n")

config_file <- here("secat_config.R")
if (file.exists(config_file)) {
  source(config_file)
  cat(sprintf("  ✓ Config loaded: MIN_REQUIRED_LENGTH = %d bp\n", MIN_REQUIRED_LENGTH))
} else {
  cat("  ⚠️ Config file not found, using defaults\n")
  MIN_REQUIRED_LENGTH <- 50
  OUTPUT_DIR <- "output/standardized_datasets"
  ALIGNED_DIR <- "output/intermediate/aligned_fastas"
}

# Minimum fraction of sequences that must survive the length filter.
# Studies falling below this are flagged as FAIL_LOW_YIELD and skipped.
# Ettinger-type failures (0.07% yield) are caught here before the sprintf crash.
if (!exists("MIN_YIELD_RATE")) MIN_YIELD_RATE <- 0.50

# Ensure directories exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- 1. Load Inputs ---
selection_file <- here("output/aggregated_data/selected_studies_for_trim.txt")
if (!file.exists(selection_file)) stop("FATAL: Selection file not found.")
selected_studies <- readLines(selection_file)
cat(sprintf("Targeting %d selected studies...\n", length(selected_studies)))

# Load study coordinates
coords_file <- here("output/intermediate/study_alignment_coords.csv")
if (!file.exists(coords_file)) stop("FATAL: Coordinates file not found.")
coords <- read_csv(coords_file, show_col_types = FALSE)

# Load consensus region
consensus_file <- here("output/intermediate/consensusregioninfo.csv")
if (!file.exists(consensus_file)) stop("FATAL: Consensus region file not found.")
consensus_info <- read_csv(consensus_file, show_col_types = FALSE)

consensus_start <- consensus_info$ConsensusStart[1]
consensus_end   <- consensus_info$ConsensusEnd[1]

if (is.na(consensus_start) || is.na(consensus_end)) {
  stop("FATAL: Consensus coordinates missing or invalid.")
}

cat(sprintf("Global Consensus Target: %d - %d (Length: %d bp)\n\n",
            consensus_start, consensus_end, consensus_end - consensus_start + 1))

# --- 2. Process Each Study ---
trim_summary <- tibble(
  study_name         = character(),
  status             = character(),
  original_seqs      = integer(),
  trimmed_seqs       = integer(),
  aligned_length     = integer(),
  degapped_length_min = integer(),
  degapped_length_max = integer()
)

for (study in selected_studies) {
  cat(sprintf("Processing: %s\n", study))

  # Get study coordinates
  study_coords <- coords %>% filter(study_name == !!study)
  if (nrow(study_coords) == 0) {
    cat("  [SKIP] No coordinates found.\n\n")
    next
  }

  # Look for aligned FASTA
  aligned_fasta <- file.path(ALIGNED_DIR, paste0(study, "_aligned.fasta"))

  if (!file.exists(aligned_fasta)) {
    cat(sprintf("  [FAIL] Aligned FASTA not found: %s\n", aligned_fasta))
    cat("  -> You need to run Phase 2 (alignment) first!\n\n")
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_NO_ALIGNMENT",
      original_seqs = NA_integer_, trimmed_seqs = 0L,
      aligned_length = NA_integer_, degapped_length_min = NA_integer_,
      degapped_length_max = NA_integer_
    ))
    next
  }

  # Load ALIGNED sequences (~43k positions wide)
  aligned_seqs    <- readDNAStringSet(aligned_fasta)
  original_count  <- length(aligned_seqs)
  cat(sprintf("  -> Loaded %d aligned sequences\n", original_count))

  alignment_widths <- unique(width(aligned_seqs))
  cat(sprintf("  -> Alignment width(s): %s\n", paste(alignment_widths, collapse = ", ")))

  study_ref_start <- study_coords$ref_start[1]
  study_ref_end   <- study_coords$ref_end[1]

  cat(sprintf("  -> Study maps to SILVA: %d - %d\n", study_ref_start, study_ref_end))
  cat(sprintf("  -> Consensus region:     %d - %d\n", consensus_start, consensus_end))

  # Check geometric overlap
  if (study_ref_end < consensus_start || study_ref_start > consensus_end) {
    cat("  [FAIL] Study does not overlap consensus region.\n\n")
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_NO_OVERLAP",
      original_seqs = original_count, trimmed_seqs = 0L,
      aligned_length = alignment_widths[1], degapped_length_min = NA_integer_,
      degapped_length_max = NA_integer_
    ))
    next
  }

  # Calculate effective overlap
  effective_start <- max(consensus_start, study_ref_start)
  effective_end   <- min(consensus_end,   study_ref_end)
  aligned_len     <- effective_end - effective_start + 1L

  cat(sprintf("  -> Effective overlap:    %d - %d\n", effective_start, effective_end))
  cat(sprintf("  -> Extracting alignment columns: %d - %d\n", effective_start, effective_end))

  # Extract consensus region (keeping gaps) then degap
  trimmed_aligned  <- subseq(aligned_seqs, start = effective_start, end = effective_end)
  trimmed_degapped <- DNAStringSet(gsub("-", "", as.character(trimmed_aligned)))

  degapped_lengths <- width(trimmed_degapped)
  min_len <- min(degapped_lengths)
  max_len <- max(degapped_lengths)
  cat(sprintf("  -> Degapped length range: %d - %d bp\n", min_len, max_len))

  # Apply minimum length filter
  cat(sprintf("  -> Applying minimum length filter: %d bp\n", MIN_REQUIRED_LENGTH))
  valid_seqs    <- trimmed_degapped[degapped_lengths >= MIN_REQUIRED_LENGTH]
  dropped_count <- original_count - length(valid_seqs)

  if (dropped_count > 0) {
    cat(sprintf("  -> Dropped %d sequences (< %d bp after degapping)\n",
                dropped_count, MIN_REQUIRED_LENGTH))
    dropped_lengths <- degapped_lengths[degapped_lengths < MIN_REQUIRED_LENGTH]
    if (length(dropped_lengths) > 0) {
      cat(sprintf("     Dropped length range: %d - %d bp (median: %.0f bp)\n",
                  min(dropped_lengths), max(dropped_lengths), median(dropped_lengths)))
    }
  }

  # --- All sequences too short ---
  if (length(valid_seqs) == 0L) {
    cat("  [FAIL] All sequences too short after trimming.\n\n")
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_TOO_SHORT",
      original_seqs = original_count, trimmed_seqs = 0L,
      aligned_length = aligned_len,
      degapped_length_min = min_len, degapped_length_max = max_len
    ))
    next
  }

  # --- Low-yield check ---
  # Catches studies whose sequences do not genuinely cover the consensus
  # region (e.g. Ettinger_2017: only 12/17675 = 0.07% survived).
  # These pass the geometric overlap test but are biologically unusable.
  yield_rate <- length(valid_seqs) / original_count
  if (yield_rate < MIN_YIELD_RATE) {
    cat(sprintf("  [FAIL] Only %.1f%% of sequences passed length filter (minimum: %.0f%%).\n",
                100 * yield_rate, 100 * MIN_YIELD_RATE))
    cat("         Sequences do not adequately cover the consensus region.\n")
    cat("         This study will be excluded from the merged dataset.\n\n")
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_LOW_YIELD",
      original_seqs = original_count, trimmed_seqs = length(valid_seqs),
      aligned_length = aligned_len,
      degapped_length_min = as.integer(min(width(valid_seqs))),
      degapped_length_max = as.integer(max(width(valid_seqs)))
    ))
    next
  }

  # --- Save output ---
  output_fasta <- file.path(OUTPUT_DIR, paste0(study, "_standardized.fasta"))
  writeXStringSet(valid_seqs, output_fasta, width = 80)

  final_lengths <- width(valid_seqs)
  cat(sprintf("  [OK] Saved %d sequences (yield: %.1f%%)\n",
              length(valid_seqs), 100 * yield_rate))
  cat(sprintf("  -> Length range: %d - %d bp (median: %.0f bp)\n",
              min(final_lengths), max(final_lengths), median(final_lengths)))
  cat(sprintf("  -> Output: %s\n\n", output_fasta))

  trim_summary <- bind_rows(trim_summary, tibble(
    study_name = study, status = "SUCCESS",
    original_seqs = original_count, trimmed_seqs = length(valid_seqs),
    aligned_length = aligned_len,
    degapped_length_min = as.integer(min(final_lengths)),
    degapped_length_max = as.integer(max(final_lengths))
  ))
}

# --- 3. Save Summary ---
summary_file <- file.path(OUTPUT_DIR, "trim_summary.csv")
write_csv(trim_summary, summary_file)
cat(sprintf("✓ Trimming summary saved: %s\n", summary_file))

# --- 4. Print Summary ---
cat("\n================================================================================\n")
cat("                        TRIMMING SUMMARY\n")
cat("================================================================================\n\n")
print(trim_summary, n = Inf)

n_success  <- sum(trim_summary$status == "SUCCESS",        na.rm = TRUE)
n_low_yield <- sum(trim_summary$status == "FAIL_LOW_YIELD", na.rm = TRUE)
n_failed   <- sum(startsWith(trim_summary$status, "FAIL"), na.rm = TRUE) - n_low_yield

cat(sprintf("\n✓ Successfully trimmed : %d / %d studies\n", n_success, length(selected_studies)))
if (n_low_yield > 0) {
  low_yield_studies <- trim_summary$study_name[trim_summary$status == "FAIL_LOW_YIELD"]
  cat(sprintf("⚠ Low yield (excluded): %d study/studies — %s\n",
              n_low_yield, paste(low_yield_studies, collapse = ", ")))
  cat("  These studies map geometrically to the consensus region but their\n")
  cat("  sequences have insufficient coverage at those alignment positions.\n")
  cat("  Review aligned FASTAs and consider excluding from future runs.\n")
}
if (n_failed > 0) {
  cat(sprintf("✗ Other failures       : %d study/studies\n", n_failed))
}

# --- 5. Launch Merger ---
if (n_success > 0) {
  cat("\n--- Launching Merger... ---\n")
  system("Rscript scripts/13_merge_datasets.R")
} else {
  cat("\n❌ No studies were successfully trimmed. Skipping merge.\n")
}

cat("\n=== STANDARDIZATION COMPLETE ===\n")
