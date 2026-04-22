#!/usr/bin/env Rscript
#===============================================================================
# SCRIPT:   scripts/01_primer_mapping.R
# PIPELINE: SeCAT (Sequence Consensus Analysis Tool)
# PHASE:    Phase 1: Primer Coordinate Mapping
# VERSION:  3.0 (Config-driven)
# AUTHOR:   [Author Name]
#
# PURPOSE:
#   Determines the physical genomic coordinates for every primer set listed in
#   the master manifest by aligning them to the full reference database.
#
# BIOLOGICAL RATIONALE:
#   To perform a meta-analysis across studies with different primer sets (e.g.,
#   V3-V4 vs V4), we must first identify the "Consensus Region"—the shared
#   intersection of all amplicons. This script maps every primer to a canonical
#   reference (e.g., SILVA/GreenGenes) to find these start/end positions.
#   This defines the maximum possible window for the subsequent analysis.
#
# INPUTS:
#   - Manifest: [SECAT_MANIFEST_PATH] (TSV with `primer_name`, `fwd_seq`, `rev_seq`)
#   - Reference DB: [REFERENCE_DB_PATH] (FASTA, e.g., SILVA 138)
#
# OUTPUTS:
#   - File: output/intermediate/primer_coords_phase1_output.csv
#   - Format: CSV with columns `primer_name`, `primer_start`, `primer_end`
#
# DEPENDENCIES:
#   - R/secat_consensus.R: Contains `find_consensus_region()` logic.
#   - R/secat_config.R: Global paths and constants.
#
#===============================================================================

#===============================================================================
# SECTION 1: INITIALIZATION
#===============================================================================

# --- Load Libraries & Config ---
# Suppress startup messages for cleaner log files in HPC environments
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Biostrings))

# Load global configuration (paths, constants)
source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R")) 

# Load core utility engines
source("R/secat_utils.R")     # General utilities
source("R/secat_consensus.R") # Contains the primer alignment logic

#===============================================================================
# SECTION 2: VALIDATION & DATA LOADING
#===============================================================================

# --- Main ---
message("--- Phase 1: Finding Primer Coordinates using Reference DB ---")

# Check that the reference database exists before attempting to load it
# This prevents obscure errors inside the alignment function
if (!file.exists(REFERENCE_DB_PATH)) {
  stop("FATAL: The main reference database was not found at the path specified in R: ", REFERENCE_DB_PATH)
}

# Load the master manifest containing the primer definitions
message(paste("Loading manifest from:", SECAT_MANIFEST_PATH))
manifest <- read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)

#===============================================================================
# SECTION 3: PRIMER ALIGNMENT
#===============================================================================

# Find primer coordinates using the full reference database
# This calls `find_consensus_region` from R/secat_consensus.R, which:
# 1. Loads the reference DB
# 2. Aligns each unique primer sequence (allowing mismatches)
# 3. Returns the modal start/end positions
primer_coords <- find_consensus_region(manifest, REFERENCE_DB_PATH, VSEARCH_PATH)

# Remove duplicates - keep only unique primer coordinates
# Scientific Note: Multiple studies in the manifest likely use the same
# standard primer sets (e.g., 515F/806R). We only need one coordinate entry
# per unique primer pair.
primer_coords <- primer_coords %>%
  distinct(primer_name, .keep_all = TRUE)

#===============================================================================
# SECTION 4: OUTPUT
#===============================================================================

# Define the output directory and file
output_dir <- "output/intermediate"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_path <- file.path(output_dir, "primer_coords_phase1_output.csv")

# Save the results for downstream use (e.g., Phase 2 simulations)
write_csv(primer_coords, output_path)

message(paste("--- Phase 1 Complete. Primer coordinates saved to:", output_path, "---"))
