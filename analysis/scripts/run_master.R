#############################
# ---- MASTER RUNNER ----   #
#############################

source(here::here("scripts", "00_setup.R"))

# --- Load configs (functions) ---
source(here::here("scripts", "config", "00_run.R"))     # run_cfg()
source(here::here("scripts", "config", "01_design.R"))  # design_cfg()
source(here::here("scripts", "config", "02_models.R"))  # model_cfg() + model functions

# --- Instantiate configs (objects) ---
run    <- run_cfg()
design <- design_cfg()
model  <- model_cfg()

# --- Bundle into a single cfg passed everywhere ---
cfg <- list(
  run    = run,
  design = design,
  model  = model
)

# ---- Load stages ----
source(here::here("scripts", "01_stage_etl.R"))
source(here::here("scripts", "02_stage_indices.R"))
source(here::here("scripts", "03_stage_analysis.R"))

# ---- Run pipeline ----

# First create csvs
run_etl(cfg) # datasets must be created first

# Define plan on what subsets of the data to run analyses on (e.g., pooled, each treatment separately, or a custom subset of treatments)
cfg$plan <- resolve_treatments(cfg$run$dataset, cfg$run$treatment) # must trail ETL, because it depends on the dataset being processed and the treatments that are present in it

# Then create indices, which depend on the datasets
run_indices(cfg)

# Finally run analyses, which depend on the datasets and indices
run_analysis(cfg)


msg("\nMASTER completed for dataset:", cfg$dataset, "\n")