# ==============================================================================
# SCRIPT:   scripts/04_prepare_sims.R
# PIPELINE: SeCAT (Sequence Consensus Analysis Tool)
# PHASE:    Phase 3 Setup: Simulation Planning
# VERSION:  13 (Robust Base R Primer Logic)
# AUTHOR:   [Author Name]
#
# PURPOSE:
#   Generates the "Task List" for the simulation phase. This script defines
#   exactly which simulations need to be run (one per study/primer x replicates).
#
# KEY OPERATIONS:
#   1. Reference Subsetting: Creates a random subset of the full reference
#      database to serve as the "Ground Truth" community for simulations.
#      Uses memory-efficient Reservoir Sampling.
#   2. Task Calculation: Determines the maximum number of trim steps possible
#      for each amplicon (Study or Primer) before length drops below minimum.
#   3. Task Expansion: Multiplies each task by the number of desired replicates
#      (e.g., 100 seeds) to create the full job manifest.
#
# INPUTS:
#   - Intermediate Outputs from Phase 1/2 (Coordinate CSVs).
#   - Full Reference DB [REFERENCE_DB_PATH].
#
# OUTPUTS:
#   - output/intermediate/simulation_tasks.csv (The master job list).
#   - output/intermediate/simulation_reference_subset.fasta (The ground truth).
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(Biostrings))

# Simple logger
log_and_flush <- function(message) {
  cat(paste(Sys.time(), "|", message, "\n"))
  flush.console()
}

# ------------------------------------------------------------------------------
# Helper: Reservoir Sampling for FASTA (Optimized)
# ------------------------------------------------------------------------------
# Why this exists: Loading a 50GB SILVA FASTA into R just to pick N sequences
# crashes most nodes. This function scans the file line-by-line (O(N)) and
# keeps a random sample in memory (O(k)), never exceeding modest RAM usage.
#
# Returns raw character vector to avoid "U" vs "T" validation errors common
# when Biostrings tries to validate RNA sequences as DNA before conversion.
# ------------------------------------------------------------------------------

sample_fasta_reservoir <- function(filepath, n_samples = 100, seed = 123) {
    set.seed(seed)

    log_and_flush(sprintf("  -> Starting reservoir sampling: targeting %d sequences", n_samples))

    con <- file(filepath, "r")
    on.exit(close(con))

    reservoir_seqs <- character(n_samples)
    reservoir_headers <- character(n_samples)
    n_filled <- 0
    total_seen <- 0

    current_header <- NULL
    current_seq_buffer <- character()

    process_sequence <- function(header, seq_fragments) {
        # Combine buffer into string
        full_seq <- paste(seq_fragments, collapse="")

        # Reservoir Logic:
        # 1. Fill the reservoir until full (first n_samples).
        # 2. For every subsequent item i, replace an existing item with probability k/i.
        if (total_seen < n_samples) {
            n_filled <<- n_filled + 1
            reservoir_headers[n_filled] <<- header
            reservoir_seqs[n_filled] <<- full_seq
        } else {
            j <- sample.int(total_seen + 1, 1)
            if (j <= n_samples) {
                reservoir_headers[j] <<- header
                reservoir_seqs[j] <<- full_seq
            }
        }
    }

    # INCREASED CHUNK SIZE for speed (from 2k to 50k lines per read)
    chunk_size <- 50000

    while(TRUE) {
        lines <- readLines(con, n = chunk_size)
        if(length(lines) == 0) break

        # Optimization: grep positions of headers in this chunk to avoid loop overhead
        header_indices <- grep("^>", lines)

        if (length(header_indices) == 0) {
            # Whole chunk is sequence data
            if (!is.null(current_header)) {
                current_seq_buffer <- c(current_seq_buffer, lines)
            }
            next
        }

        # Process mixed chunk (headers + sequences)
        prev_idx <- 1
        for (idx in header_indices) {
            # Content before this header belongs to previous seq
            if (idx > 1) {
                if (!is.null(current_header)) {
                    current_seq_buffer <- c(current_seq_buffer, lines[prev_idx:(idx-1)])
                }
            } else if (idx == 1 && !is.null(current_header)) {
                # Header is first line, flush previous if exists
                process_sequence(current_header, current_seq_buffer)
                total_seen <- total_seen + 1
            }

            if (idx > 1 && !is.null(current_header)) {
                process_sequence(current_header, current_seq_buffer)
                total_seen <- total_seen + 1
            }

            # Start new sequence
            current_header <- sub("^>", "", lines[idx])
            current_seq_buffer <- character()
            prev_idx <- idx + 1
        }

        # Trailing content in chunk
        if (prev_idx <= length(lines)) {
            current_seq_buffer <- c(current_seq_buffer, lines[prev_idx:length(lines)])
        }

        # Log less frequently to avoid I/O delay
        if (total_seen %% 50000 == 0) {
            cat(sprintf("\r    Scanning sequence %d...", total_seen))
        }
    }

    # Handle last sequence in file
    if (!is.null(current_header)) {
        process_sequence(current_header, current_seq_buffer)
        total_seen <- total_seen + 1
    }
    cat("\n")

    # Return RAW LIST (Character Vector) to safe-guard against RNA/DNA mismatch
    valid_indices <- which(reservoir_headers != "")

    log_and_flush(sprintf("  -> Sampling complete: %d sequences selected from %d total",
                          length(valid_indices), total_seen))

    return(list(
        seqs = reservoir_seqs[valid_indices],
        headers = reservoir_headers[valid_indices],
        total_scanned = total_seen
    ))
}

# ------------------------------------------------------------------------------
# Main Logic
# ------------------------------------------------------------------------------
main <- function() {
    log_and_flush("--- SCRIPT STARTED: prepare_simulation_tasks.R V13 ---")

    source(here::here("secat_config.R"))
    log_and_flush(paste("ANALYSIS_MODE detected:", ANALYSIS_MODE))

    # Load Configs (with defaults)
    TRIM_INCREMENT <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 10
    MIN_FINAL_SEQUENCE_LENGTH <- if (exists("MIN_FINAL_SEQUENCE_LENGTH")) MIN_FINAL_SEQUENCE_LENGTH else 100
    MAX_ABSOLUTE_TRIM_STEPS <- if (exists("MAX_ABSOLUTE_TRIM_STEPS")) MAX_ABSOLUTE_TRIM_STEPS else 2000
    NUM_SIMULATIONS_PER_PRIMER <- if (exists("NUM_SIMULATIONS_PER_PRIMER")) NUM_SIMULATIONS_PER_PRIMER else 100
    STEP_MODE <- if(exists("TRIM_STEP_MODE")) TRIM_STEP_MODE else "scaled"
    DEFAULT_STEPS <- if(exists("DEFAULT_MAX_TRIM_STEPS")) DEFAULT_MAX_TRIM_STEPS else 40
    BUFFER_STEPS <- if(exists("CONSENSUS_BUFFER_STEPS")) CONSENSUS_BUFFER_STEPS else 10

    # NEW: Get subset size from config
    SUBSET_SIZE <- if(exists("SIMULATION_MAX_SILVA_SUBSET")) SIMULATION_MAX_SILVA_SUBSET else 10000

    manifest <- readr::read_tsv(here::here(SECAT_MANIFEST_PATH), show_col_types = FALSE)
    simulation_tasks <- tibble()

    # =========================================================
    # STEP 1: Generate Reference Subset (RESERVOIR SAMPLING)
    # =========================================================
    # Only strictly necessary if we are in Study Mode OR if Primer Mode actually needs it later.
    # The original script ran this for Study Mode. We will keep logic consistent.
    
    if (TRUE) {
        log_and_flush("\n=== GENERATING SILVA SUBSET FOR SIMULATIONS ===")
        subset_path <- here::here("output/intermediate/simulation_reference_subset.fasta")

        if (!file.exists(subset_path)) {
            log_and_flush(sprintf("  -> Performing Reservoir Sampling (%d seqs) from Reference DB...", SUBSET_SIZE))
            if (!file.exists(REFERENCE_DB_PATH)) stop("FATAL: Reference DB not found at: ", REFERENCE_DB_PATH)

            tryCatch({
                # Run reservoir sampler (Returns Raw Character List)
                sample_res <- sample_fasta_reservoir(REFERENCE_DB_PATH, n_samples = SUBSET_SIZE, seed = 42)

                # Clean chars (U -> T) BEFORE creating DNAStringSet
                # This ensures RNA DBs (like SILVA) are converted to DNA for VSEARCH
                raw_seqs <- sample_res$seqs
                cleaned_seqs <- gsub(".", "-", raw_seqs, fixed = TRUE)
                cleaned_seqs <- gsub("U", "T", cleaned_seqs, ignore.case = TRUE)

                # Now safely create DNAStringSet
                ref_clean <- Biostrings::DNAStringSet(cleaned_seqs)
                names(ref_clean) <- sample_res$headers

                Biostrings::writeXStringSet(ref_clean, subset_path)

                # Report file size and statistics
                file_size_mb <- file.info(subset_path)$size / 1024^2
                avg_length <- mean(Biostrings::width(ref_clean))

                log_and_flush(paste("  -> Saved simulation reference subset to:", subset_path))
                log_and_flush(sprintf("  -> Subset contains %d sequences (%.1f MB)",
                                     length(ref_clean), file_size_mb))
                log_and_flush(sprintf("  -> Average sequence length: %.0f bp", avg_length))

            }, error = function(e) {
                stop(paste("Failed to generate reference subset:", e$message))
            })
        } else {
            # Verify existing subset
            existing_seqs <- Biostrings::readDNAStringSet(subset_path)
            file_size_mb <- file.info(subset_path)$size / 1024^2

            log_and_flush("  -> Simulation reference subset already exists. Skipping generation.")
            log_and_flush(sprintf("  -> Using existing subset: %d sequences (%.1f MB)",
                                 length(existing_seqs), file_size_mb))

            if (length(existing_seqs) < SUBSET_SIZE * 0.9) {
                warning(sprintf("Existing subset has only %d sequences but config specifies %d. Consider regenerating.",
                               length(existing_seqs), SUBSET_SIZE))
            }
        }
    }

    # =========================================================
    # STEP 2: Task Generation Logic
    # =========================================================
    if (ANALYSIS_MODE == "primer") {
        log_and_flush("\n=== Executing in 'primer' mode. ===")
        coords_path <- here::here("output/intermediate/primer_coords_phase1_output.csv")
        primer_coords <- readr::read_csv(coords_path, show_col_types = FALSE)
        all_primers <- unique(primer_coords$primer_name)

        # Helper for Primer Absolute Mode: Estimates amplicon length from FASTA
        # FIXED: Uses Base R subsetting to avoid dplyr scoping errors
        calc_abs <- function(p_name_val) {
             # Explicit subsetting - safe against NSE issues
             studies <- manifest[manifest$primer_name == p_name_val, ]
             
             if(nrow(studies)==0) return(30)
             max_len <- 0
             for(i in 1:nrow(studies)) {
                 fp <- studies$asv_fasta_path[i]
                 if(file.exists(fp)) {
                     try({ s <- readDNAStringSet(fp); max_len <- max(max_len, max(width(s))) }, silent=TRUE)
                 }
             }
             if(max_len==0) return(30)
             # Formula: (Length - Min_Required) / Increment = Max Steps
             return(min(floor(max(0, max_len - MIN_FINAL_SEQUENCE_LENGTH)/TRIM_INCREMENT), MAX_ABSOLUTE_TRIM_STEPS))
        }

        # Generate Steps Table
        primer_steps_table <- purrr::map_dfr(all_primers, function(curr_primer) {
            steps <- if (STEP_MODE == "scaled") DEFAULT_STEPS + BUFFER_STEPS else calc_abs(curr_primer)
            tibble(primer_name = curr_primer, max_trim_steps = steps)
        })

        readr::write_csv(primer_steps_table, here::here("output/intermediate/primer_max_trim_steps.csv"))

        # Create Task Matrix: Primers x Replicates (Base R Version)
        # FIXED: Replaces dplyr joins with Base R merge/expand.grid to prevent "object not found" errors
        
        # 1. Prepare Lengths
        primer_coords$amplicon_length <- primer_coords$primer_end - primer_coords$primer_start + 1
        unique_lengths <- unique(primer_coords[, c("primer_name", "amplicon_length")])
        
        # 2. Merge Steps + Lengths
        merged_tasks <- merge(primer_steps_table, unique_lengths, by="primer_name")
        colnames(merged_tasks)[colnames(merged_tasks) == "primer_name"] <- "task_id"
        colnames(merged_tasks)[colnames(merged_tasks) == "max_trim_steps"] <- "num_steps"
        
        # 3. Expand Replicates
        seeds <- 1:NUM_SIMULATIONS_PER_PRIMER
        simulation_tasks <- expand.grid(task_id = merged_tasks$task_id, simulation_seed = seeds, stringsAsFactors = FALSE)
        
        # 4. Finalize
        simulation_tasks <- merge(simulation_tasks, merged_tasks, by="task_id")
        simulation_tasks <- simulation_tasks[, c("task_id", "num_steps", "amplicon_length", "simulation_seed")]
        simulation_tasks <- as_tibble(simulation_tasks)

    } else if (ANALYSIS_MODE == "study") {
        log_and_flush("\n=== Executing in 'study' mode. ===")
        coords_path <- here::here("output/intermediate/study_alignment_coords.csv")
        study_coords <- readr::read_csv(coords_path, show_col_types = FALSE)

        # Create Task Matrix: Studies x Replicates
        simulation_tasks <- study_coords %>%
            dplyr::mutate(
                task_id = study_name,
                num_steps = if (STEP_MODE == "scaled") {
                     DEFAULT_STEPS + BUFFER_STEPS
                } else {
                    # Empirical step calculation based on alignment length
                    steps_calc <- floor((analysis_amplicon_length - MIN_FINAL_SEQUENCE_LENGTH) / TRIM_INCREMENT)
                    pmin(pmax(0, steps_calc), MAX_ABSOLUTE_TRIM_STEPS)
                }
            ) %>%
            dplyr::select(task_id, num_steps, amplicon_length = analysis_amplicon_length) %>%
            tidyr::crossing(simulation_seed = 1:NUM_SIMULATIONS_PER_PRIMER)

        # === NEW: Generate and Save Consensus Region Info ===
        log_and_flush("--- Calculating and saving consensus region info ---")

        # Load the secat_consensus.R functions
        source(here::here("R/secat_consensus.R"))

        # Calculate global consensus
        consres <- find_largest_overlapping_clique(
          starts = study_coords$ref_start,
          ends = study_coords$ref_end,
          study_names = study_coords$study_name,
          min_overlap = 50
        )

        # Create consensus info dataframe
        consensus_info <- tibble::tibble(
          ConsensusStart = if(!is.null(consres$start)) consres$start else NA,
          ConsensusEnd = if(!is.null(consres$end)) consres$end else NA,
          ConsensusLength = if(!is.null(consres$start)) (consres$end - consres$start) else 0,
          NumStudiesInConsensus = if(!is.null(consres$n_studies)) consres$n_studies else 0,
          IncludedStudies = if(!is.null(consres$included_studies)) paste(consres$included_studies, collapse = ";") else "",
          OutlierStudies = if(!is.null(consres$excluded_studies)) paste(consres$excluded_studies, collapse = ";") else ""
        )

        # Write to file
        consensus_path <- here::here("output/intermediate/consensusregioninfo.csv")
        readr::write_csv(consensus_info, consensus_path)
        log_and_flush(sprintf("✓ Saved consensus region info to %s", consensus_path))

        # Log details
        if (!is.null(consres$start)) {
            log_and_flush(sprintf("  Consensus: %d-%d bp (%d studies)",
                                 consres$start, consres$end, consres$n_studies))
            log_and_flush(sprintf("  Outliers: %s",
                                 ifelse(length(consres$excluded_studies) > 0,
                                        paste(consres$excluded_studies, collapse = ", "),
                                        "None")))
        } else {
            log_and_flush("  WARNING: No valid consensus clique could be found.")
        }
    }

    # Save Master Task List (CRITICAL: Must run for BOTH modes)
    output_path <- here::here("output/intermediate/simulation_tasks.csv")
    readr::write_csv(simulation_tasks, output_path)
    log_and_flush(sprintf("Saved %d simulation tasks to %s", nrow(simulation_tasks), output_path))
    log_and_flush(sprintf("  -> %d unique tasks × %d replicates each",
                         length(unique(simulation_tasks$task_id)),
                         NUM_SIMULATIONS_PER_PRIMER))
}
# Run Main
tryCatch(main(), error = function(e) {
    log_and_flush(paste("ERROR:", conditionMessage(e)))
    quit(save = "no", status = 1)
})
