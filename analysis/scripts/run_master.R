# ============================================================
# master_runner.R
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
#     ../raw_data/<data_folder>/   (set data_folder in scripts/config/00_run.R)
#
# PIPELINE STAGES
#   1. ETL     -- reads raw oTree exports, structures and cleans data
#   2. Indices -- computes participant-level quantities (r_i, a*, drift)
#   3. Analysis -- fits Stan models, produces tables and figures
#
# OUTPUT
#   All outputs are written to analysis/data/<data_folder>/
#     etl/      -- cleaned input files
#     output/   -- CSV tables and model summaries
#     models/   -- Stan fit objects (RDS)
#     figures/  -- paper figures (PNG)
#
# CONFIGURATION
#   scripts/config/00_run.R    -- execution controls (dataset, seed, overwrite)
#   scripts/config/01_design.R -- experimental design and inferential thresholds
#   scripts/config/02_models.R -- Stan sampling settings
#
# DEPENDENCIES
#   R >= 4.2, rstan >= 2.21, data.table, ggplot2, ggdist, patchwork
#   See 00_setup.R for full package list.
#
# ============================================================

# ---- Configs ----
source(here::here("scripts", "config", "00_run.R"))
source(here::here("scripts", "config", "01_design.R"))
source(here::here("scripts", "config", "02_models.R"))

# ---- Pure functions ----
source(here::here("scripts", "config", "03_eu_functions.R"))

# ---- Bundle cfg ----
cfg <- list(
  run    = run_cfg(),
  design = design_cfg(),
  model  = model_cfg()
)

# ---- Setup (requires cfg for paths) ----
source(here::here("scripts", "00_setup.R"))
validate_cfg(cfg)


# ---- Raw data check ----
if (!dir.exists(path_raw)) {
  stop(
    "\nRaw data folder not found: ", path_raw,
    "\nCreate the folder and place oTree export files there before running.\n"
  )
}
if (length(list.files(path_raw)) == 0) {
  stop(
    "\nRaw data folder is empty: ", path_raw,
    "\nPlace oTree export files there before running.\n"
  )
}


# ---- Load stages ----
source(here::here("scripts", "01_stage_etl.R"))
source(here::here("scripts", "02_stage_indices.R"))
source(here::here("scripts", "03_stage_analysis.R"))

# ---- Run pipeline ----
run_etl(cfg)
run_indices(cfg)
run_analysis(cfg)

msg("\nMASTER completed for data_folder:", cfg$run$data_folder, "\n")