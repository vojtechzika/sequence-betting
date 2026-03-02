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

# ============================================================
# HELPERS
# ============================================================

### Console Messages
msg <- function(...) cat(..., "\n")

### Validate config files
validate_cfg <- function(cfg) {
  
  stopifnot(is.list(cfg))
  
  # ---- run ----
  stopifnot(
    !is.null(cfg$run),
    nzchar(as.character(cfg$run$dataset)),
    is.numeric(cfg$run$seed),
    length(cfg$run$seed) == 1L,
    !is.null(cfg$run$treatment)
  )
  
  # ---- design ----
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
  
  # ---- model ----
  stopifnot(!is.null(cfg$model))
  
  invisible(TRUE)
}


### rewritting existing files


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