# =============================================================================
# Script:  00_setup.R
# Purpose: Shared setup for the review/ folder. Source this at the top of
#          every other script in scripts/ so that package loading, the
#          random seed, and the session log are handled consistently.
# Author:  Antho
# Project: TE_Project / review (reviewer response, Genome Biology submission)
# =============================================================================

# --- Reproducibility ---------------------------------------------------------
# Fix the seed once, here, so every downstream script inherits the same
# random state for anything stochastic (permutations, bootstraps, etc.).
set.seed(20260824)

# --- Paths --------------------------------------------------------------------
# All paths are relative to the review/ project root (open review.Rproj in
# RStudio, or setwd() to this folder, before sourcing this script).
dir_data    <- "data"
dir_outputs <- "outputs"
dir_logs    <- "logs"

dir.create(dir_data,    showWarnings = FALSE, recursive = TRUE)
dir.create(dir_outputs, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_logs,    showWarnings = FALSE, recursive = TRUE)

# --- Packages -------------------------------------------------------------
# Keep this list explicit and minimal; add packages here as scripts need
# them rather than loading them ad hoc inside individual scripts.
required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "tibble",
  "stringr",
  "ggplot2"
)
# Script-specific dependencies are loaded inside the script that needs them
# (e.g. terra in 02_candidate_TEs.R, Biostrings in
# 05_tfbs_enrichment.R) so that a missing optional package does
# not block the other scripts.

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# --- Session log ------------------------------------------------------------
# Call this at the END of each script (not here) with a script-specific tag,
# e.g. log_session("02_candidate_TEs"), so every run is traceable to the
# exact package/R versions used to produce its outputs.
log_session <- function(tag) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(dir_logs, sprintf("%s_%s_sessionInfo.txt", tag, ts))
  writeLines(capture.output(sessionInfo()), log_file)
  message("Session info written to: ", log_file)
}
