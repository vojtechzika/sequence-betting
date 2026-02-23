############################
# ---- Analysis Stage ----
############################

run_analysis <- function(cfg) {
  
  stopifnot(
    is.list(cfg),
    !is.null(cfg$run),
    !is.null(cfg$run$dataset),
    !is.null(cfg$plan)
  )
  
  ds   <- cfg$run$dataset
  plan <- cfg$plan  # required for downstream funcs
  
  # --- Load Sources --- #
  source(here::here("scripts", "analysis", "01_descriptive_participants.R"))
  source(here::here("scripts", "analysis", "02_descriptive_sequences.R"))
  
  source(here::here("scripts", "analysis", "11_rq1_stan.R"))
  source(here::here("scripts", "analysis", "12_rq1_sequences.R"))
  source(here::here("scripts", "analysis", "13_rq1_participants.R"))
  
  source(here::here("scripts", "analysis", "21_rq2_stan.R"))
  source(here::here("scripts", "analysis", "22_rq2_sequences.R"))
  source(here::here("scripts", "analysis", "23_rq2_participants.R"))
  
  # ----------------------------
  # 00 Space: Descriptives
  # ----------------------------
  descriptive_participants(cfg)
  descriptive_sequences(cfg)
  
  # ----------------------------
  # 10 Space: RQ1
  # ----------------------------
  rq1_stan(cfg)
  rq1_sequences(cfg)
  rq1_participants(cfg)
  
  # ----------------------------
  # 20 Space: RQ2 (placeholder)
  # ----------------------------
   rq2_stan(cfg)
   rq2_sequences(cfg)
   rq2_participants(cfg)
   
  
  msg("\nAnalysis phase completed for:", ds, "\n")
}