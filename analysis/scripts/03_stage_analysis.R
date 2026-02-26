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
  
  source(here::here("scripts", "analysis", "31_rq3_stan.R"))
  source(here::here("scripts", "analysis", "32_rq3_sequences.R"))
  source(here::here("scripts", "analysis", "33_rq3_participants.R"))
  
  source(here::here("scripts", "analysis", "41_rq4_stan.R"))
  source(here::here("scripts", "analysis", "42_rq4_sequences.R"))
  source(here::here("scripts", "analysis", "43_rq4_participants.R"))
  
  source(here::here("scripts", "analysis", "51_ex1_1_anchors.R"))
  source(here::here("scripts", "analysis", "52_ex1_1_sequences.R"))
  source(here::here("scripts", "analysis", "53_ex1_1_participants.R"))
  source(here::here("scripts", "analysis", "54_ex1_2_stan_ghi.R"))
  source(here::here("scripts", "analysis", "55_ex1_3_ghi_vs_pure.R"))

  source(here::here("scripts", "analysis", "61_ex2_stan.R"))
  source(here::here("scripts", "analysis", "62_ex2_summary.R"))
  source(here::here("scripts", "analysis", "63_ex2_participants.R"))
  
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
  # 20 Space: RQ2 
  # ----------------------------
   rq2_stan(cfg)
   rq2_sequences(cfg)
   rq2_participants(cfg)
   
   # ----------------------------
   # 30 Space: RQ3 
   # ----------------------------
   rq3_stan(cfg)
   rq3_sequences(cfg)
   rq3_participants(cfg)
   
   # ----------------------------
   # 40 Space: RQ4
   # ----------------------------
   rq4_stan(cfg)
   rq4_sequences(cfg)
   rq4_participants(cfg)
   
   # ----------------------------
   # 50 Space: EX1
   # ----------------------------
   ex1_1_anchors(cfg)
   ex1_1_sequences(cfg)
   ex1_1_participants(cfg)
   ex1_2_stan_ghi(cfg)
   ex1_3_ghi_vs_pure(cfg)
   
   # ----------------------------
   # 60 Space: EX2
   # ----------------------------
   ex2_stan(cfg)
   ex2_summary(cfg)
   ex2_participants(cfg)
   
   
  
  msg("\nAnalysis phase completed for:", ds, "\n")
}