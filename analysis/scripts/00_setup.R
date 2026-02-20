options(stringsAsFactors = FALSE)
options(width = 120)

# ---- Packages (minimal, add later when needed) ----
pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "here", "lubridate")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Paths ----
# analysis/ is the R project root; repo root is one level up
path_repo <- here::here("..")

# Raw + clean roots
path_raw   <- file.path(path_repo, "data", "raw")
path_clean <- file.path(path_repo, "data", "clean")

# Dataset-specific roots
path_raw_ds   <- function(ds) file.path(path_raw, ds)
path_clean_ds <- function(ds) file.path(path_clean, ds)

# Analysis artifacts (dataset-specific, kept under data/clean/<ds>/)
path_out_ds <- function(ds) file.path(path_clean_ds(ds), "output")
path_fig_ds <- function(ds) file.path(path_clean_ds(ds), "figures")
path_mod_ds <- function(ds) file.path(path_clean_ds(ds), "models")

# ---- Create dirs ----
dir.create(path_raw,   showWarnings = FALSE, recursive = TRUE)
dir.create(path_clean, showWarnings = FALSE, recursive = TRUE)

for (ds in c("pilot", "main")) {
  dir.create(path_raw_ds(ds),   showWarnings = FALSE, recursive = TRUE)
  dir.create(path_clean_ds(ds), showWarnings = FALSE, recursive = TRUE)
  
  dir.create(path_out_ds(ds), showWarnings = FALSE, recursive = TRUE)
  dir.create(path_fig_ds(ds), showWarnings = FALSE, recursive = TRUE)
  dir.create(path_mod_ds(ds), showWarnings = FALSE, recursive = TRUE)
}

# ---- Helpers ----
msg <- function(...) cat(..., "\n")

# ----------------------------
# Resolve treatment plan
# ----------------------------
# treatment options:
#   "all"  = pooled only (no filtering)
#   "each" = run separately for each treatment present
#   "both" = pooled + separately for each treatment present
#   c("m25","m19") = specified subset (no pooled)
resolve_treatments <- function(dataset, treatment) {
  
  f_merged <- file.path(path_clean_ds(dataset), "merged.csv")
  if (!file.exists(f_merged)) {
    stop(
      "merged.csv not found at: ", f_merged, "\n",
      "Treatment resolution requires merged.csv. Call resolve_treatments() only after ETL has created it."
    )
  }
  
  dt <- data.table::fread(f_merged, select = "participant.treatment")
  available <- sort(unique(as.character(dt$participant.treatment)))
  available <- available[!is.na(available) & nzchar(available)]
  
  if (length(available) == 0) {
    stop("No non-missing participant.treatment values found in merged.csv.")
  }
  
  if (identical(treatment, "all"))  return(list(pooled = TRUE,  by = character(0), available = available))
  if (identical(treatment, "both")) return(list(pooled = TRUE,  by = available,    available = available))
  if (identical(treatment, "each")) return(list(pooled = FALSE, by = available,    available = available))
  
  if (!is.character(treatment)) {
    stop("Invalid treatment specification. Use 'all', 'each', 'both', or a character vector like c('m25','m19').")
  }
  
  missing <- setdiff(treatment, available)
  if (length(missing) > 0) {
    stop(
      "Unknown treatment(s): ", paste(missing, collapse = ", "),
      "\nAvailable: ", paste(available, collapse = ", ")
    )
  }
  
  list(pooled = FALSE, by = treatment, available = available)
}


# treatment-specific multiplier lookup (design spec)
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