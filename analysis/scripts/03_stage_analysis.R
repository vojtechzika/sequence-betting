############################
# ---- Analysis Stage ----
############################

run_analysis <- function(cfg) {
  
  stopifnot(is.list(cfg), "dataset" %in% names(cfg), "plan" %in% names(cfg))
  dataset <- cfg$dataset
  plan <- cfg$plan
  
  # If you want stage files runnable standalone, keep this:
  source(here::here("scripts", "00_setup.R"))
  
  # --- Load Sources --- #
  source(here::here("scripts", "analysis", "01_descriptive_participants.R"))
  source(here::here("scripts", "analysis", "02_descriptive_sequences.R"))
  
  source(here::here("scripts", "analysis", "11_rq1_sequences.R"))
  source(here::here("scripts", "analysis", "12_rq1_participants.R"))
  
  # ----------------------------
  # 00 Space: Descriptives
  # ----------------------------
  descriptive_participants(cfg)
  descriptive_sequences(cfg)
  
  # ----------------------------
  # 10 Space: RQ1
  # ----------------------------
  rq1_sequences(cfg)
  rq1_participants(cfg)
  
  # ----------------------------
  # 20 Space: RQ2 (placeholder)
  # ----------------------------
  # rq2_sequences(cfg)
  # rq2_participants(cfg)
  
  msg("\nAnalysis phase completed for:", dataset, "\n")
}