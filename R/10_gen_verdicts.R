#!/usr/bin/env Rscript
# scripts/10_gen_verdicts.R
# PURPOSE: Prepare raw verdict data for interactive selection

suppressPackageStartupMessages({
  suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
  library(here)
})

cat("--- Generating Verdict Data ---\n")

# Load raw verdict data from aggregation
master_verdict_file <- here("output/aggregated_data/master_verdict_table.csv")

if (!file.exists(master_verdict_file)) {
  stop("FATAL: master_verdict_table.csv not found. Run aggregation first.")
}

verdicts <- read_csv(master_verdict_file, show_col_types = FALSE)

cat(sprintf("Loaded %d verdict rows covering %d studies and %d levels.\n", 
            nrow(verdicts),
            length(unique(verdicts$Study)),
            length(unique(verdicts$Level))))

# Save as intermediate for the selector to use
write_csv(verdicts, here("output/aggregated_data/verdict_data_all_levels.csv"))

cat("✓ Verdict data saved for interactive selection.\n")
