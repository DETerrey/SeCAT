#===============================================================================
# FILENAME:   R/secat_mapping.R
# PIPELINE:   SeCAT (Sequence Consensus Analysis Tool)
# VERSION:    1.0.0
# AUTHOR:     [Author Name]
#
# PURPOSE:
#   Determines the genomic coordinates (start/end) for a given study relative
#   to a reference database.
#
# STRATEGIES:
#   1. STUDY MODE (Recommended):
#      Uses Profile Alignment (via the DECIPHER package) to map a subset of
#      actual ASV sequences to a reference alignment (e.g., SILVA). This
#      empirically determines the region covered by the study.
#
#   2. PRIMER MODE (Legacy/Fallback):
#      Parses theoretical coordinates directly from primer names (e.g.,
#      extracting "515" and "806" from "515F_806R"). This relies on standardized
#      naming conventions and does not validate the actual sequences.
#
# INPUTS:
#   - Study Info (List): Contains file paths and metadata.
#   - Reference DB (FASTA): A pre-aligned reference database.
#   - Config (List): Parameters for alignment and subsampling.
#
# OUTPUTS:
#   - Summary List containing:
#     * summary: Data frame with study-level consensus coordinates.
#     * coords: Data frame with per-ASV coordinates.
#     * alignment: The DNAStringSet of aligned ASVs.
#===============================================================================

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(tibble))

#-------------------------------------------------------------------------------
# Function:   calculate_mode
#-------------------------------------------------------------------------------
# Description:
#   Calculates the statistical mode. (Duplicate utility).
#-------------------------------------------------------------------------------
calculate_mode <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#-------------------------------------------------------------------------------
# Function:   parse_primer_positions (Legacy Helper)
#-------------------------------------------------------------------------------
# Description:
#   Extracts numeric coordinates from standard primer names (e.g., "515F_806R").
#-------------------------------------------------------------------------------
parse_primer_positions <- function(primer_name) {
    # Allow separators like "_", "-", or "." and extra prefixes
    parts <- strsplit(primer_name, "[-_.]")[[1]]
    parts <- parts[nzchar(parts)]

    # Find the first forward- and reverse-like tokens
    fwd_token <- parts[grep("F$", parts, ignore.case = TRUE)[1]]
    rev_token <- parts[grep("R$", parts, ignore.case = TRUE)[1]]

    get_pos <- function(x) {
        if (is.na(x)) return(NA_integer_)
        m <- regmatches(x, regexpr("[0-9]+", x))
        if (length(m) > 0) as.integer(m) else NA_integer_
    }

    fwd_pos <- get_pos(fwd_token)
    rev_pos <- get_pos(rev_token)

    list(start = fwd_pos, end = rev_pos)
}

#-------------------------------------------------------------------------------
# Function:   calculate_consensus_coordinates
#-------------------------------------------------------------------------------
# Description:
#   Aggregates individual ASV coordinates into a single consensus window for
#   the study.
#
# Methods:
#   - "modal" (Default): The most common start/end positions. Best for
#     amplicon data where most reads should start at the exact primer site.
#   - "percentile": Uses 10th/90th percentiles. Robust to outliers/artifacts.
#   - "minimum": Extends from min start to max end. Maximizes coverage but
#     riskier with noisy data.
#-------------------------------------------------------------------------------
calculate_consensus_coordinates <- function(all_ref_starts, all_ref_ends, method = "modal") {
  if (method == "modal") {
    start <- calculate_mode(all_ref_starts)
    end <- calculate_mode(all_ref_ends)
  } else if (method == "percentile") {
    start <- quantile(all_ref_starts, 0.10, type = 1)
    end <- quantile(all_ref_ends, 0.90, type = 1)
  } else if (method == "minimum") {
    start <- min(all_ref_starts)
    end <- max(all_ref_ends)
  } else {
    warning(paste("Unknown method:", method, ". Defaulting to 'modal'"))
    start <- calculate_mode(all_ref_starts)
    end <- calculate_mode(all_ref_ends)
  }
  return(list(start = start, end = end))
}

#-------------------------------------------------------------------------------
# Function:   map_study_to_reference
#-------------------------------------------------------------------------------
# Description:
#   The core mapping function. Coordinates the entire workflow of loading data,
#   aligning (or parsing), and returning coordinates.
#
# Scientific Context:
#   Accurate coordinate mapping is essential for meta-analysis. If Study A
#   claims to be "V4" but actually covers V3-V4 due to a protocol modification,
#   relying on metadata alone would lead to invalid biological comparisons.
#   This function allows us to verify the *actual* coverage empirically.
#
# Algorithm (Study Mode):
#   1. Subsampling: Randomly selects N ASVs (default 1000) to represent the study.
#      This is computationally efficient and statistically sufficient to find the
#      consensus window.
#   2. Reference Prep: Loads a pre-aligned reference (e.g., SILVA), optionally
#      subsetting it to speed up profile alignment.
#   3. Alignment: Uses `DECIPHER::AlignProfiles`. This treats the reference as
#      a fixed profile and aligns the ASVs to it, preserving the reference's
#      coordinate system.
#   4. Coordinate Extraction: Converts the alignment gaps (`-`) into numerical
#      start/end indices relative to the reference columns.
#
# Parameters:
#   @param study_info [list] - Metadata including `asv_fasta_path`.
#   @param reference_db_path [character] - Path to reference alignment.
#   @param config [list] - Configuration options (e.g., `ASV_SAMPLE_SIZE`).
#
# Returns:
#   @return [list] - Detailed results object with summary stats and the raw alignment.
#-------------------------------------------------------------------------------
map_study_to_reference <- function(study_info, reference_db_path, config) {

  message(paste("Processing study:", study_info$study_name))

  # Determine Mode from Config (Default to "study" if missing)
  mode <- if(!is.null(config$ANALYSIS_MODE)) config$ANALYSIS_MODE else "study"
  message(paste("  - Analysis Mode:", mode))

  # --- LEGACY PATH: PRIMER MODE ---
  # Trusts the primer name string to define coordinates.
  if (mode == "primer") {
      message("  - Using theoretical primer coordinates (Legacy Mode).")

      # 1. Parse primer positions
      primer_positions <- parse_primer_positions(study_info$primer_name)
      if (is.na(primer_positions$start) || is.na(primer_positions$end)) {
        warning(paste("Could not parse positions from primer name:", study_info$primer_name))
        return(NULL)
      }

      # 2. Return theoretical coordinates wrapped in list structure
      # We assume these apply directly to the reference genome/alignment
      summary_df <- tibble::tibble(
        study_name = study_info$study_name,
        primer_name = study_info$primer_name,
        ref_start_method = "PRIMER_THEORETICAL",
        ref_start = primer_positions$start,
        ref_end = primer_positions$end,
        analysis_amplicon_length = primer_positions$end - primer_positions$start + 1,
        original_start_theoretical = primer_positions$start,
        initial_fwd_trim = study_info$initial_fwd_trim,
        initial_rev_trim = study_info$initial_rev_trim,
        num_asvs_mapped = NA # Not mapping ASVs
      )

      return(list(
        summary = summary_df,
        coords = NULL,
        alignment = NULL
      ))
  }

  # --- MODERN PATH: STUDY/ALIGNMENT MODE (DECIPHER) ---
  # Empirically aligns ASVs to the reference profile.

  # 1. Check Dependencies
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    stop("Package 'DECIPHER' is required but not installed.")
  }

  # 2. Validate Study ASVs
  if (!file.exists(study_info$asv_fasta_path)) {
    warning(paste("FASTA not found for", study_info$study_name, ". Skipping."))
    return(NULL)
  }

  asv_sequences <- Biostrings::readDNAStringSet(study_info$asv_fasta_path)
  if (length(asv_sequences) == 0) {
    warning(paste("FASTA is empty for", study_info$study_name, ". Skipping."))
    return(NULL)
  }
  message(paste("  - Loaded", length(asv_sequences), "ASV sequences."))

  # 3. ASV Subsampling Logic
  # Aligning millions of ASVs is slow; a random subset defines the window accurately.
  if (config$USE_ALL_ASVS_FOR_MAFFT) {
     message("  - Using ALL ASVs for alignment (High Compute Mode).")
  } else {
     sample_size <- min(config$ASV_SAMPLE_SIZE, length(asv_sequences))
     if (length(asv_sequences) > sample_size) {
       set.seed(42)
       asv_sample_indices <- sample(seq_along(asv_sequences), sample_size)
       asv_sequences <- asv_sequences[asv_sample_indices]
       message(paste("  - Subsampled to", length(asv_sequences), "ASVs for mapping."))
     }
  }

  # 4. Robust Reference Loading
  message("  - Loading Reference Database...")
  # Using readBStringSet allows flexibility (works for RNA or DNA files)
  ref_raw_bstring <- Biostrings::readBStringSet(reference_db_path)

  # Optionally subset the reference to speed up Profile creation
  if (config$REFERENCE_ALIGNMENT_MODE == "full") {
      message("  - Using FULL reference database as profile.")
      ref_subset_bstring <- ref_raw_bstring
  } else {
      subset_size <- config$REFERENCE_SUBSET_SIZE
      actual_size <- min(subset_size, length(ref_raw_bstring))
      message(paste("  - Using SUBSET of reference database (", actual_size, "sequences) as profile."))
      set.seed(123)
      ref_indices <- sample(seq_along(ref_raw_bstring), actual_size)
      ref_subset_bstring <- ref_raw_bstring[ref_indices]
  }

  message("  - Standardizing reference profile format...")
  # Ensure standard DNA alphabet (T not U) and gap characters (-)
  clean_seqs <- function(bstring) {
      seqs_char <- as.character(bstring)
      seqs_char <- gsub(".", "-", seqs_char, fixed = TRUE)
      seqs_char <- gsub("U", "T", seqs_char, ignore.case = TRUE)
      return(seqs_char)
  }
  cleaned_chars <- clean_seqs(ref_subset_bstring)
  ref_profile <- Biostrings::DNAStringSet(cleaned_chars)

  # 5. Align ASVs to Profile
  message("  - Aligning ASVs to SILVA Profile using DECIPHER...")
  # Degap ASVs to treat them as unaligned queries
  asv_clean <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(asv_sequences)))
  names(asv_clean) <- names(asv_sequences)

  # Step 1: Align ASVs to each other (temporarily)
  asv_profile <- DECIPHER::AlignSeqs(asv_clean, verbose = FALSE)
  # Step 2: Align the ASV profile to the Reference profile
  merged_alignment <- DECIPHER::AlignProfiles(pattern = asv_profile, subject = ref_profile)

  # Filter to keep only the study sequences (exclude reference profile sequences)
  aligned_asvs <- merged_alignment[names(asv_clean)]

  # 6. Extract Coordinates from Alignment
  message("  - Extracting coordinates...")
  # Helper to find first and last non-gap character
  get_coords <- function(seq_str) {
    chars <- strsplit(seq_str, "")[[1]]
    non_gaps <- which(chars != "-")
    if (length(non_gaps) == 0) return(c(NA, NA))
    return(c(min(non_gaps), max(non_gaps)))
  }

  seq_strs <- as.character(aligned_asvs)
  coords_list <- lapply(seq_strs, get_coords)

  all_ref_starts <- sapply(coords_list, `[`, 1)
  all_ref_ends   <- sapply(coords_list, `[`, 2)

  valid_idx <- !is.na(all_ref_starts) & !is.na(all_ref_ends)
  all_ref_starts <- all_ref_starts[valid_idx]
  all_ref_ends <- all_ref_ends[valid_idx]

  if (length(all_ref_starts) == 0) {
    warning("No valid mapping coordinates found.")
    return(NULL)
  }

  # 7. Calculate Consensus Region for this Study
  consensus <- calculate_consensus_coordinates(all_ref_starts, all_ref_ends, method = config$STUDY_ALIGNMENT_METHOD)
  modal_start <- consensus$start
  modal_end <- consensus$end
  analysis_length <- modal_end - modal_start + 1

  # 8. Build Output Summary
  output_summary <- tibble::tibble(
    study_name = study_info$study_name,
    primer_name = study_info$primer_name,
    ref_start_method = paste0("DECIPHER_SILVA_", config$STUDY_ALIGNMENT_METHOD),
    ref_start = modal_start,
    ref_end = modal_end,
    analysis_amplicon_length = analysis_length,
    original_start_theoretical = NA,
    initial_fwd_trim = study_info$initial_fwd_trim,
    initial_rev_trim = study_info$initial_rev_trim,
    num_asvs_mapped = length(all_ref_starts)
  )

  # 9. Build Detailed Coordinate Table (NEW)
  # Preserves per-ASV data for downstream debugging or detailed analysis
  detailed_coords <- tibble::tibble(
    study_name = study_info$study_name,
    asv_id = names(asv_clean)[valid_idx],
    start = all_ref_starts,
    end = all_ref_ends
  )

  message(paste("  - Done. SILVA Column Region:", modal_start, "-", modal_end))

  return(list(
    summary = output_summary,
    coords = detailed_coords,
    alignment = aligned_asvs
  ))
}

#-------------------------------------------------------------------------------
# Function:  relax_consensus_coords
#-------------------------------------------------------------------------------
# Description:
#  Relax invalid consensus coordinates with deterministic fallback
#' Ensures Start < End for ALL downstream consumers (viz, aggregation, verdicts)
#' @param clique_result List from find_largest_overlapping_clique()
#' @param all_starts,all_ends Original study coordinate vectors  
#' @param study_names All study names
#' @return Fixed clique_result with valid start/end
#' @export
#   Calculates the statistical mode. (Duplicate utility).
#-------------------------------------------------------------------------------
relax_consensus_coords <- function(clique_result, all_starts, all_ends, study_names, study_name = "consensus") {
  cons_start <- clique_result$start
  cons_end <- clique_result$end
  
  if (!is.na(cons_start) && !is.na(cons_end) && cons_start >= cons_end) {
    message(sprintf("[WARN] Invalid Consensus Detected (Start: %d >= End: %d) for %s. Applying smart relaxation...", 
                    cons_start, cons_end, study_name))
    
    # Use ALL study bounds as conservative fallback (not just clique)
    study_min <- min(all_starts, na.rm = TRUE)
    study_max <- max(all_ends, na.rm = TRUE)
    
    # Ensure minimum 50bp window matching min_overlap
    relaxed_start <- min(study_min, study_max - 50)
    relaxed_end <- max(study_max, study_min + 50)
    
    clique_result$start <- relaxed_start
    clique_result$end <- relaxed_end
    clique_result$n_studies <- length(study_names)  # Flag as "relaxed consensus"
    
    message(sprintf("  -> Relaxed to study bounds [%d,%d] (length: %d bp)", 
                    relaxed_start, relaxed_end, relaxed_end - relaxed_start))
  }
  return(clique_result)
}
