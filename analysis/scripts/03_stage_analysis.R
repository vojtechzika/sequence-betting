############################
# ---- Analysis Stage ----
############################

run_analysis <- function(cfg) {

  # --- Load Sources --- #
  source(here::here("scripts", "analysis", "01_descriptive_participants.R"))
  source(here::here("scripts", "analysis", "02_descriptive_sequences.R"))
  source(here::here("scripts", "analysis", "03_r_raincloud.R"))
  source(here::here("scripts", "analysis", "04_raw_bets_and_stakes_by_trt.R"))
  source(here::here("scripts", "analysis", "05_raw_stakes_sequences.R"))
  source(here::here("scripts", "analysis", "06_desc_times_and_optimism_by_trt.R"))
  
  source(here::here("scripts", "analysis", "11_rq1_stan.R"))
  source(here::here("scripts", "analysis", "12_rq1_diagnostics.R"))
  source(here::here("scripts", "analysis", "13_rq1_tables.R")) 
  source(here::here("scripts", "analysis", "14_rq1_figures.R")) 

  source(here::here("scripts", "analysis", "21_rq2_stan.R"))
  source(here::here("scripts", "analysis", "22_rq2_diagnostics.R"))
  source(here::here("scripts", "analysis", "23_rq2_tables.R"))
  source(here::here("scripts", "analysis", "24_rq2_figures.R")) 
  source(here::here("scripts", "analysis", "25_rq2_stake_vs_r.R")) 
  
  source(here::here("scripts", "analysis", "31_rq3_stan.R"))
  source(here::here("scripts", "analysis", "32_rq3_diagnostics.R"))
  source(here::here("scripts", "analysis", "33_rq3_tables.R"))
  source(here::here("scripts", "analysis", "34_rq3_figures.R"))
  source(here::here("scripts", "analysis", "39_rq3_floor_robustness.R"))
  
  #these two could be merged into one
  source(here::here("scripts", "analysis", "35_rq3_plot_forest.R")) 
  source(here::here("scripts", "analysis", "36_rq3_plot_raincloud.R"))
  
 
  source(here::here("scripts", "analysis", "41_rq4_stan.R"))
  source(here::here("scripts", "analysis", "42_rq4_diagnostics.R"))
  source(here::here("scripts", "analysis", "43_rq4_tables.R"))
  source(here::here("scripts", "analysis", "44_rq4_figures.R"))
   
  source(here::here("scripts", "analysis", "51_ex1_1_anchors.R"))
  source(here::here("scripts", "analysis", "52_ex1_1_diagnostics.R"))
  source(here::here("scripts", "analysis", "53_ex1_1_tables.R"))
  source(here::here("scripts", "analysis", "54_ex1_2_stan_ghi.R"))
  source(here::here("scripts", "analysis", "55_ex1_2_diagnostics.R"))
  source(here::here("scripts", "analysis", "56_ex1_3_ghi_vs_pure.R"))
  source(here::here("scripts", "analysis", "57_ex1_desc_pure_summary.R"))
  source(here::here("scripts", "analysis", "58_ex1_desc_pure_vs_covariates.R"))
  source(here::here("scripts", "analysis", "59_ex1_desc_pure_vs_mirror.R"))
   
  source(here::here("scripts", "analysis", "61_ex2_stan.R"))
  source(here::here("scripts", "analysis", "62_ex2_diagnostics.R"))
  source(here::here("scripts", "analysis", "63_ex2_tables.R"))
  
  source(here::here("scripts", "analysis", "71_ex3_contrasts.R"))
  #source(here::here("scripts", "analysis", "72_ex3_rq2_implied_p.R"))
  #source(here::here("scripts", "analysis", "73_ex3_implied_r.R"))
  source(here::here("scripts", "analysis", "74_ex3_figures.R"))

  source(here::here("scripts", "analysis", "81_ex4_similarity.R"))
  source(here::here("scripts", "analysis", "82_ex4_canonical_rules.R"))
  
  # SYNTHESIS
  source(here::here("scripts", "analysis", "91_synthesis_sequences.R"))
  source(here::here("scripts", "analysis", "92_synthesis_mirror_plot.R"))
  source(here::here("scripts", "analysis", "93_synthesis_mechanism_checks.R"))
  source(here::here("scripts", "analysis", "94_synthesis_stan.R"))
  source(here::here("scripts", "analysis", "95_synthesis_stan_outputs.R"))
  source(here::here("scripts", "analysis", "96_synthesis_debrief_validation.R"))
  source(here::here("scripts", "analysis", "97_synthesis_out_of_sample.R"))
  source(here::here("scripts", "analysis", "98_synthesis_cost_vs_direction.R"))
  source(here::here("scripts", "analysis", "99_synthesis_reaction_times.R"))
  
  # DISCUSSION
  source(here::here("scripts", "analysis", "101_discussion.R"))
  
  
  
  # ----------------------------
  # 00 Space: Descriptives
  # ----------------------------
  descriptive_participants(cfg)
  descriptive_sequences(cfg)
  r_raincloud(cfg)
  raw_bets_and_stakes_by_trt(cfg)
  raw_stakes_sequences(cfg)
  desc_times_and_optimism_by_trt(cfg)

  
  # ----------------------------
  # 10 Space: RQ1
  # ----------------------------
  rq1_stan(cfg)                      
  rq1_diagnostics(cfg)              
  rq1_tables(cfg)
  rq1_figures(cfg)
  
  # ----------------------------
  # 20 Space: RQ2 
  # ----------------------------
   rq2_stan(cfg)
   rq2_diagnostics(cfg)
   rq2_tables(cfg)
   rq2_figures(cfg)
   rq2_stake_vs_r(cfg)
   
   # ----------------------------
   # 30 Space: RQ3 
   # ----------------------------
   rq3_stan(cfg)
   rq3_diagnostics(cfg)
   rq3_tables(cfg)
   rq3_figures(cfg)
   
   rq1_3_forest(cfg)
   rq1_3_raincloud(cfg)
   
   rq3_floor_robustness(cfg)
   
   # ----------------------------
   # 40 Space: RQ4
   # ----------------------------
   rq4_stan(cfg)
   rq4_diagnostics(cfg)
   rq4_tables(cfg)
   rq4_figures(cfg)
   
   # ----------------------------
   # 50 Space: EX1
   # ----------------------------
   ex1_1_anchors(cfg)
   ex1_1_diagnostics(cfg)
   ex1_1_tables(cfg)
   ex1_2_stan_ghi(cfg)
   ex1_2_diagnostics(cfg)
   ex1_3_ghi_vs_pure(cfg)
   ex1_pure_summary(cfg)
   ex1_pure_vs_covariates(cfg)
   ex1_pure_vs_mirror(cfg)
   
   # ----------------------------
   # 60 Space: EX2
   # ----------------------------
   ex2_stan(cfg)
   ex2_diagnostics(cfg)
   ex2_tables(cfg)
   
   # ----------------------------
   # 70 Space: EX3
   # ----------------------------
   ex3_contrasts(cfg)
   #ex3_rq2_implied_p(cfg)
   #ex3_implied_r(cfg)
   ex3_figures(cfg)
   
   # ----------------------------
   # 80 Space: EX4
   # ----------------------------
   ex4_similarity(cfg)
   ex4_canonical_rules(cfg)
   
   # ----------------------------
   # 90 Space: Synthesis
   # ----------------------------
   synthesis_sequences(cfg)
   synthesis_mirror_plot(cfg)
   synthesis_mechanism_checks(cfg)
   synthesis_stan(cfg)
   synthesis_stan_outputs(cfg)
   synthesis_debrief_validation(cfg)
   synthesis_out_of_sample(cfg)
   synthesis_cost_vs_direction(cfg)
   synthesis_reaction_times(cfg)
   

   # ----------------------------
   # 100 Space: Discussion
   # ----------------------------
   discussion(cfg)


  
  msg("\nANALYSIS completed for data_folder:", cfg$run$data_folder, "\n")
}