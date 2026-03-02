#############################
# ---- MASTER RUNNER ----   #
#############################

source(here::here("scripts", "00_setup.R"))

# ---- Configs ----
source(here::here("scripts", "config", "00_run.R"))
source(here::here("scripts", "config", "01_design.R"))
source(here::here("scripts", "config", "02_models.R"))

# ---- Pure functions ----
source(here::here("scripts", "config", "03_eu_functions.R"))

# --- Bundle into a single cfg passed everywhere ---
cfg <- list(
  run    = run_cfg(),
  design = design_cfg(),
  model  = model_cfg()
)

validate_cfg(cfg)

# ---- Load stages ----
source(here::here("scripts", "01_stage_etl.R"))
source(here::here("scripts", "02_stage_indices.R"))
source(here::here("scripts", "03_stage_analysis.R"))

# ---- Run pipeline ----

run_etl(cfg) # extract data from oTree source and populate raw tables
run_indices(cfg) # Continue with indices the analyses will need
run_analysis(cfg) # Run analyses


msg("\nMASTER completed for dataset:", cfg$dataset, "\n")