options(stringsAsFactors = FALSE)
options(width = 120)

# ---- Packages (minimal, add later when needed) ----
pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "here", "lubridate")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Paths ----
# analysis/ is the R project root; repo root is one level up
path_repo  <- here::here("..")

path_raw   <- file.path(path_repo, "data", "raw")
path_clean <- file.path(path_repo, "data", "clean")
path_clean_ds <- function(ds) file.path(path_clean, ds)

# analysis artifacts (dataset-specific, kept under data/clean/<ds>/)
path_out_ds <- function(ds) file.path(path_clean_ds(ds), "output")
path_fig_ds <- function(ds) file.path(path_clean_ds(ds), "figures")
path_mod_ds <- function(ds) file.path(path_clean_ds(ds), "models")

# ---- Create dirs ----
dir.create(path_raw, showWarnings = FALSE, recursive = TRUE)
dir.create(path_clean, showWarnings = FALSE, recursive = TRUE)

dir.create(path_clean_ds("pilot"), showWarnings = FALSE, recursive = TRUE)
dir.create(path_clean_ds("main"),  showWarnings = FALSE, recursive = TRUE)

dir.create(path_out_ds("pilot"), showWarnings = FALSE, recursive = TRUE)
dir.create(path_fig_ds("pilot"), showWarnings = FALSE, recursive = TRUE)
dir.create(path_mod_ds("pilot"), showWarnings = FALSE, recursive = TRUE)

dir.create(path_out_ds("main"), showWarnings = FALSE, recursive = TRUE)
dir.create(path_fig_ds("main"), showWarnings = FALSE, recursive = TRUE)
dir.create(path_mod_ds("main"), showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----
msg <- function(...) cat(..., "\n")