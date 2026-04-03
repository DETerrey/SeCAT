
#!/usr/bin/env Rscript
# =====================================================================
# WRAPPER SCRIPT: render_wrapper.R (Definitive)
#
# PURPOSE:
# This script acts as a robust bridge between the SGE shell script and
# RMarkdown. It takes a study name as a simple command-line argument
# and passes it correctly to rmarkdown::render, solving all parameter
# passing issues.
# =====================================================================

# Load the 'here' library to build robust, absolute paths
suppressPackageStartupMessages(library(here))

# Get the study name from the first command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("FATAL: No study_name was provided to the render wrapper script.")
}
study_name_param <- args[1]

# Print a message to the log to confirm we received the parameter
message(paste("Wrapper script received study_name:", study_name_param))

# --- Call rmarkdown::render with the correct, absolute paths ---
# This is now done safely inside a pure R environment.
rmarkdown::render(
  input = here("scripts", "generate_report.Rmd"),
  params = list(study_name = study_name_param),
  output_file = here("output", "reports", paste0(study_name_param, "_report.html"))
)

message("Report rendering complete.")

