#!/usr/bin/env Rscript
# =====================================================================
# WORKER SCRIPT: generate_primer_db.R (V4 - Memory Optimized Streaming)
#
# PURPOSE:
# Takes a single task ID from SGE and generates the corresponding
# primer-specific database using chunked/streaming processing to
# minimize memory usage (~2GB instead of 20GB+).
# =====================================================================

# --- MEMORY OPTIMIZATION: Force garbage collection ---
gc(reset = TRUE, full = TRUE)

# --- Load Libraries & Config ---
suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
})

# Load config WITHOUT here:: to avoid loading extra packages
if (file.exists("secat_config.R")) {
  source("secat_config.R")
} else {
  source(here::here("secat_config.R"))  # Fallback
}

# Explicit GC after config load
gc(full = TRUE)

# --- Get Task ID from SGE ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("FATAL: SGE_TASK_ID not provided.")
task_id <- as.integer(args[1])
message(paste("Starting worker for Task ID:", task_id))

# --- Define Input/Output Paths ---
PRIMER_COORDS_CSV <- "output/intermediate/primer_coords_phase1_output.csv"
OUTPUT_DIR <- "output/primer_databases"

# --- Read Primer Info for This Task ---
if (!file.exists(PRIMER_COORDS_CSV)) {
  stop("FATAL: Primer coordinates file not found: ", PRIMER_COORDS_CSV)
}

primer_sets <- read_csv(PRIMER_COORDS_CSV, show_col_types = FALSE)
if (task_id > nrow(primer_sets)) {
  stop(paste("FATAL: Task ID", task_id, "is out of bounds."))
}

primer_info <- primer_sets[task_id, ]

primer_name <- primer_info$primer_name
start_pos <- primer_info$primer_start
end_pos <- primer_info$primer_end

output_fasta <- file.path(OUTPUT_DIR, paste0("db_", primer_name, ".fasta"))
message(paste("  -> Target Primer:", primer_name, "| Coords:", start_pos, "-", end_pos))

# Clean up
rm(primer_sets, primer_info)
gc(full = TRUE)

# --- Check Database Exists ---
if (!file.exists(REFERENCE_DB_PATH)) {
  stop("FATAL: Reference database not found at: ", REFERENCE_DB_PATH)
}

# --- Memory-Efficient Streaming Extraction ---
message("Using memory-efficient chunked processing...")

# Create FASTA index (fast, ~1 second, minimal memory)
message("  -> Building FASTA index...")
fai <- Biostrings::fasta.index(REFERENCE_DB_PATH)
n_seqs <- nrow(fai)
message(paste("     Found", n_seqs, "sequences in reference database"))

# Validate first sequence length
message("  -> Checking sequence format...")
first_seq <- readBStringSet(REFERENCE_DB_PATH, nrec = 1)
full_len <- width(first_seq)[1]
message(paste("     Reference alignment length:", full_len, "bp"))

if (full_len < end_pos) {
  stop(paste("FATAL: Sequences (", full_len, "bp) are shorter than primer end position (", end_pos, ")!"))
}

rm(first_seq)
gc()

# Define chunk size (50k seqs = ~500MB RAM per chunk)
chunk_size <- 50000
n_chunks <- ceiling(n_seqs / chunk_size)

# Prepare output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Clear output file if exists
if (file.exists(output_fasta)) {
  message("  -> Removing existing output file...")
  file.remove(output_fasta)
}

message(paste("Processing", n_seqs, "sequences in", n_chunks, "chunks of", chunk_size, "..."))

# --- Streaming Processing Loop ---
total_extracted <- 0
start_time <- Sys.time()

for (chunk_id in 1:n_chunks) {
  # Calculate sequence range for this chunk
  start_idx <- (chunk_id - 1) * chunk_size + 1
  end_idx <- min(chunk_id * chunk_size, n_seqs)
  chunk_length <- end_idx - start_idx + 1
  
  # Read chunk using BStringSet (handles gaps/alignment characters)
  chunk_seqs <- tryCatch({
    readBStringSet(REFERENCE_DB_PATH, nrec = chunk_length, skip = start_idx - 1)
  }, error = function(e) {
    message(paste("WARNING: Error reading chunk", chunk_id, ":", conditionMessage(e)))
    return(NULL)
  })
  
  if (is.null(chunk_seqs) || length(chunk_seqs) == 0) {
    message(paste("  Skipping empty chunk", chunk_id))
    next
  }
  
  # Filter: Keep only sequences long enough
  all_widths <- width(chunk_seqs)
  valid_indices <- which(all_widths >= end_pos)
  
  if (length(valid_indices) > 0) {
    valid_seqs <- chunk_seqs[valid_indices]
    
    # Extract amplicon region (primer coordinates)
    amplicons <- subseq(valid_seqs, start = start_pos, end = end_pos)
    
    # Write to output file (append mode)
    writeXStringSet(amplicons, output_fasta, append = TRUE)
    
    total_extracted <- total_extracted + length(amplicons)
    
    # Clean up chunk
    rm(valid_seqs, amplicons)
  }
  
  # Progress report every 5 chunks
  if (chunk_id %% 5 == 0 || chunk_id == n_chunks) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- chunk_id / elapsed
    eta_secs <- (n_chunks - chunk_id) / rate
    
    message(sprintf("  Progress: %d/%d chunks (%.1f%%) | Extracted: %d seqs | ETA: %.1f min",
                    chunk_id, n_chunks, 100 * chunk_id / n_chunks, 
                    total_extracted, eta_secs / 60))
  }
  
  # Cleanup
  rm(chunk_seqs, all_widths, valid_indices)
  gc()
}

# --- Final Report ---
elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

message("========================================")
message(paste("Successfully wrote:", output_fasta))
message(paste("Total extracted sequences:", total_extracted))
message(paste("Processing time:", round(elapsed_total / 60, 2), "minutes"))
message(paste("Extraction rate:", round(n_seqs / elapsed_total), "sequences/sec"))
message("========================================")

# Final cleanup
rm(fai)
gc()
