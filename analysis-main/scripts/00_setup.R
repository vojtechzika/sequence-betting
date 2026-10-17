# ============================================================
# 00_setup.R
# ============================================================
#
# PURPOSE
#   Central setup file sourced once in run-analysis.R after cfg is assembled.
#   Loads packages, builds paths from config, creates directories,
#   and provides shared utility functions for the pipeline.
#
# FOLDER STRUCTURE
#   <project root>/             # R project root (sequences.Rproj)
#     raw-data/otree/           # oTree exports (.csv); every .csv in the folder is read
#     raw-data/manual/          # manual input files (.csv)
#     data/
#       source/                 # tables extracted unchanged from the raw export
#       derived/                # tables computed from source/ and read by later scripts
#       models/                 # one subfolder per Stan model
#       output/                 # results meant for reading (csv/, tables/, figures/, reports/)
#     scripts/
#     selected-models.csv       # model folder used for each model key
#
# USAGE
#   Sourced once in run-analysis.R after cfg is assembled.
#   All scripts access paths via global path_* variables defined here.
#   The project analyses one dataset: to analyse another, copy the project
#   and replace the files in raw-data/otree/.
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
# Attached globally
pkgs_attach <- c(
  # Data manipulation and project structure
  "data.table", "here",
  # Bayesian modelling
  "rstan",
  # Figures
  "ggplot2", "ggdist", "patchwork", "scales",
  # Utilities
  "parallel"
)

# Installed, but called via pkg:: or loaded inside the scripts that use them
pkgs_use <- c(
  "psych",                                          # LOT-R reliability (21)
  "ggrepel", "gridExtra",                            # figures
  "mclust", "sandwich", "lmtest",                   # mechanism checks, RT (93, 99)
  "lme4", "lmerTest", "performance",                # trial-level GLMMs (95b, 99)
  "marginaleffects", "stargazer"                    # 95b
)

pkgs <- c(pkgs_attach, pkgs_use)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs_attach, library, character.only = TRUE))

# ============================================================
# PATHS
# ============================================================

# Raw data folder with two subfolders: otree/ holds the oTree exports (every
# .csv in it is read) and manual/ the manually prepared input files
# (hand-coded debrief strategies, out-of-sample sequences)
path_raw    <- here::here("raw-data", "otree")
path_manual <- here::here("raw-data", "manual")

# Generated data
.data_dir <- here::here("data")

path_src    <- file.path(.data_dir, "source")         # tables extracted from the raw export
path_mod    <- file.path(.data_dir, "models")         # one subfolder per model
path_der    <- file.path(.data_dir, "derived")        # tables computed from source/ and read by later scripts
path_out    <- file.path(.data_dir, "output")         # results meant for reading
path_csv    <- file.path(path_out, "csv")
path_tab    <- file.path(path_out, "tables") # LaTeX tables only (.tex); csv files go to path_csv
path_fig    <- file.path(path_out, "figures")
path_rep    <- file.path(path_out, "reports")      # text reports on the run (environment, function tests)

# ============================================================
# CREATE DIRECTORIES
# ============================================================
# Note: path_raw is not created here -- it must already exist with raw data
for (.p in c(path_src, path_out, path_mod, path_fig, path_tab, path_der, path_csv, path_rep)) {
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


# ============================================================
# environment
# Writes the computational environment of a run and the raw data files it
# was run on (name, size, modification time, MD5 checksum) to a text file.
# OUTPUT: path_rep/run-environment.txt
# ============================================================
write_environment <- function(cfg) {
  
  f_out <- file.path(path_rep, "run-environment.txt")
  
  # Raw data files the run read (every .csv in raw-data/otree/ and raw-data/manual/)
  raw_files <- c(list.files(path_raw, pattern = "\\.csv$", full.names = TRUE),
                 list.files(path_manual, pattern = "\\.csv$", full.names = TRUE))
  raw_info  <- if (length(raw_files) == 0L) {
    "(none found)"
  } else {
    fi <- file.info(raw_files)
    vapply(seq_along(raw_files), function(i) {
      paste0(basename(raw_files[i]),
             "  |  ", format(fi$size[i], big.mark = ",", scientific = FALSE), " bytes",
             "  |  modified ", format(fi$mtime[i], "%Y-%m-%d %H:%M:%S"),
             "  |  MD5 ", unname(tools::md5sum(raw_files[i])))
    }, character(1))
  }
  
  git_hash <- tryCatch(
    trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  git_dirty <- tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L,
    error = function(e) NA, warning = function(w) NA
  )
  
  lines <- c(
    "RUN ENVIRONMENT",
    paste0("Written:        ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Seed:           ", cfg$run$seed),
    paste0("Treatments:     ", paste(cfg$run$treatment, collapse = ", ")),
    paste0("Git commit:     ", git_hash, if (isTRUE(git_dirty)) "  (uncommitted changes)" else ""),
    "",
    "RAW DATA (raw-data/otree/, raw-data/manual/)",
    raw_info,
    "",
    paste0("R:              ", R.version.string),
    paste0("rstan:          ", as.character(utils::packageVersion("rstan"))),
    paste0("StanHeaders:    ", as.character(utils::packageVersion("StanHeaders"))),
    paste0("Stan (stanc):   ", tryCatch(rstan::stan_version(), error = function(e) NA_character_)),
    paste0("Platform:       ", R.version$platform),
    "",
    "STAN SAMPLER SETTINGS (cfg$model$stan)",
    capture.output(str(cfg$model$stan)),
    "",
    "SIMULATION SETTINGS (cfg$model$simulation)",
    capture.output(str(cfg$model$simulation)),
    "",
    "PACKAGE VERSIONS",
    vapply(pkgs, function(p) paste0(p, ": ", as.character(utils::packageVersion(p))), character(1)),
    "",
    "SESSION INFO",
    capture.output(sessionInfo())
  )
  
  writeLines(lines, f_out)
  msg("Saved: ", f_out)
  invisible(f_out)
}
