# ============================================================
# 00_setup.R
# ============================================================
#
# PURPOSE
#   Central setup file sourced once in master_runner.R after cfg is assembled.
#   Loads packages, builds paths from config, creates directories,
#   and provides shared utility functions for the pipeline.
#
# FOLDER STRUCTURE
#   Sequences_2026/
#     raw_data/
#       <data_folder>/        # raw oTree exports, set in 00_run_cfg.R
#     analysis/               # git repo, R project root
#       data/
#         <data_folder>/
#           etl/              # ETL outputs: sequences.csv, participants.csv, etc.
#           output/           # derived quantities: model summaries, tables
#           models/           # Stan fit objects and auxiliary RDS files
#           figures/          # paper figures
#       scripts/
#       stan/
#
# USAGE
#   Sourced once in master_runner.R after cfg is assembled.
#   All scripts access paths via global path_* variables defined here.
#   Switching datasets: change cfg$run$data_folder in 00_run_cfg.R.
#
# DEPENDENCIES
#   Packages listed in pkgs below; installed automatically if missing.
#
# ============================================================

options(stringsAsFactors = FALSE)
options(width = 120)

# ============================================================
# PACKAGES
# ============================================================
pkgs <- c(
  # Core data manipulation
  "data.table", "dplyr", "tidyr", "stringr", "lubridate", "psych",
  # Project structure
  "here",
  # Bayesian modelling
  "rstan",
  # Figures
  "ggplot2", "ggdist", "patchwork", "scales",
  # Utilities
  "parallel"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# PATHS
# ============================================================

# Repository root is one level above analysis/
path_repo <- here::here("..")

# Active dataset folder from config (e.g. "main" or "pilot")
stopifnot(!is.null(cfg$run$data_folder), nzchar(cfg$run$data_folder))
.ds <- as.character(cfg$run$data_folder)

# Data directories
path_raw <- file.path(here::here(".."), "raw_data", cfg$run$data_folder)
path_src <- file.path(here::here("data"), cfg$run$data_folder, "etl")
path_out <- file.path(here::here("data"), cfg$run$data_folder, "output")
path_mod <- file.path(here::here("data"), cfg$run$data_folder, "models")
path_fig <- file.path(here::here("data"), cfg$run$data_folder, "figures")

# ============================================================
# CREATE DIRECTORIES
# ============================================================
# Note: path_raw is not created here -- it must already exist with raw data
for (.p in c(path_src, path_out, path_mod, path_fig)) {
  dir.create(.p, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# HELPERS
# ============================================================

# ---- Console messages ----
msg <- function(...) cat(..., "\n")

# ---- Config validation ----
validate_cfg <- function(cfg) {
  
  stopifnot(is.list(cfg))
  
  # run
  stopifnot(
    !is.null(cfg$run),
    nzchar(as.character(cfg$run$data_folder)),
    is.numeric(cfg$run$seed),
    length(cfg$run$seed) == 1L,
    !is.null(cfg$run$treatment)
  )
  
  # design
  stopifnot(
    !is.null(cfg$design),
    !is.null(cfg$design$seq),
    !is.null(cfg$design$seq$treatments)
  )
  
  tr <- unique(as.character(cfg$run$treatment))
  stopifnot(
    length(tr) > 0L,
    all(nzchar(tr)),
    all(tr %in% names(cfg$design$seq$treatments))
  )
  
  # model
  stopifnot(!is.null(cfg$model))
  
  invisible(TRUE)
}

# ---- Skip existing artifacts ----
should_skip <- function(paths, cfg, type = c("model", "output"), label = "artifact") {
  
  type <- match.arg(type)
  stopifnot(is.character(paths), length(paths) >= 1L)
  
  overwrite_flag <- switch(
    type,
    model  = isTRUE(cfg$run$overwrite_models),
    output = isTRUE(cfg$run$overwrite_outputs)
  )
  
  exists_any <- any(file.exists(paths))
  
  if (exists_any && !overwrite_flag) {
    hint <- switch(
      type,
      model  = "Set cfg$run$overwrite_models = TRUE to regenerate.",
      output = "Set cfg$run$overwrite_outputs = TRUE to regenerate."
    )
    warning(
      label, " already exists. Skipping.\n",
      paste0("  ", paths, collapse = "\n"),
      "\n", hint
    )
    return(TRUE)
  }
  
  FALSE
}

# ---- Treatment multiplier lookup ----
get_mult <- function(design, tr) {
  
  stopifnot(is.list(design), !is.null(design$seq), !is.null(design$seq$treatments))
  
  tr <- as.character(tr)
  if (length(tr) != 1L || !nzchar(tr)) stop("Treatment must be a single non-empty string.")
  
  if (!tr %in% names(design$seq$treatments)) {
    stop(
      "Unknown treatment '", tr, "'. Available in design$seq$treatments: ",
      paste(names(design$seq$treatments), collapse = ", ")
    )
  }
  
  m <- design$seq$treatments[[tr]]
  
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m <= 1) {
    stop("Invalid multiplier for treatment '", tr, "': ", paste(m, collapse = ", "))
  }
  
  m
}