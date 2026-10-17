# ============================================================
# run-analysis.R
# ============================================================
#
# PURPOSE
#   Top-level entry point for the analysis pipeline. Sources all
#   configuration files, assembles the cfg object, initializes the
#   environment, and runs the pipeline stages in order.
#
# USAGE
#   Open sequences.Rproj in RStudio and run this file.
#   Before running, ensure raw data is in place:
#     raw-data/otree/   (oTree exports, .csv)
#     raw-data/manual/  (manual input files, .csv)
#
# PIPELINE STAGES
#   1. Processing -- reads raw oTree exports, structures and cleans data
#   2. Modelling  -- fits the Stan models (risk preference r, a*, drift,
#                    analysis samples, RQ models)
#   3. Reporting  -- audits the selected models, produces tables and figures
#
# OUTPUT
#   All outputs are written to data/
#     source/   -- tables extracted from the raw export
#     derived/  -- tables computed from source/ and read by later scripts
#     models/   -- one subfolder per Stan model
#     output/   -- csv/, tables/ and figures/ meant for reading
#
# CONFIGURATION
#   settings.R                 -- run controls (treatments, seed, overwrite) and
#                                 Stan sampling settings and priors
#   scripts/config/01_design.R -- experimental design and inferential thresholds
#   selected-models.csv        -- model folder used for each model key
#
# DEPENDENCIES
#   R >= 4.2, rstan >= 2.21, data.table, ggplot2, ggdist, patchwork
#   See 00_setup.R for full package list.
#
# ============================================================

# ---- Settings and design ----
source(here::here("settings.R"))
source(here::here("scripts", "config", "01_design.R"))

# ---- Pure functions ----
source(here::here("scripts", "config", "02_eu-functions.R"))

# ---- Model helpers ----
source(here::here("scripts", "modelling", "_model-helpers.R"))

# ---- Bundle cfg ----
cfg <- list(
  run        = settings$run,
  design     = design_cfg(),
  model      = settings$model,
  synthesis  = settings$synthesis,
  discussion = settings$discussion
)

# ---- Setup (requires cfg for paths) ----
source(here::here("scripts", "00_setup.R"))
validate_cfg(cfg)


# ---- Function tests ----
# Validates the expected-utility functions at the design values; writes
# data/output/reports/eu-functions-test.txt and stops if a check fails
source(here::here("scripts", "config", "03_eu-functions-test.R"))
eu_functions_test()

# ---- Raw data check ----
if (!dir.exists(path_raw)) {
  stop(
    "\nRaw data folder not found: ", path_raw,
    "\nCreate raw-data/otree/ and place the oTree export files there before running.\n"
  )
}
if (length(list.files(path_raw)) == 0) {
  stop(
    "\nRaw data folder is empty: ", path_raw,
    "\nPlace oTree export files there before running.\n"
  )
}

# ---- Environment record ----
write_environment(cfg)

# ---- Load stages ----
source(here::here("scripts", "01_run-processing.R"))
source(here::here("scripts", "02_run-modelling.R"))
source(here::here("scripts", "03_run-reporting.R"))

# ---- Run pipeline ----
run_processing(cfg)
run_modelling(cfg)
run_reporting(cfg)

msg("\nRun completed.\n")