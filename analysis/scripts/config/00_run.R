# ============================================================
# 00_run_cfg.R
# Execution-level controls only
# ============================================================

run_cfg <- function() {
  list(
    # Dataset to analyse: "pilot" or "main"
    dataset   = "main",
    
    # Treatments to include (first = confirmatory)
    treatment = c("m25", "m19"),
    
    # Overwrite behavior
    overwrite_outputs = FALSE,
    overwrite_models  = FALSE,
    
    # Global seed for reproducibility
    seed = 12345
  )
}