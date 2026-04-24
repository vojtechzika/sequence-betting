# ============================================================
# 00_run_cfg.R
# Execution-level controls only
# ============================================================
run_cfg <- function() {
  list(
    
    # DATA FOLDER
    # Must match the name of the raw data folder in ../raw_data/
    # Workflow: (1) create ../raw_data/<name>/ and place oTree exports there,
    #           (2) set data_folder to the same name here.
    data_folder = "mar2026",
    
    # Treatments to analyze (first = confirmatory)
    treatment = c("m25", "m19"),
    
    # Restrict RQ2/RQ3 to MPL-consistent participants only
    consistent_only = TRUE,
    
    # Overwrite behavior
    overwrite_outputs = TRUE,
    overwrite_models  = TRUE,
    
    # Global seed for reproducibility
    seed = 12345
  )
}