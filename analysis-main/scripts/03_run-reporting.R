############################
# ---- Reporting Stage ----
############################

run_reporting <- function(cfg) {
  
  # ----------------------------
  # Helpers
  # ----------------------------
  source(here::here("scripts", "reporting", "_reporting-style.R"))
  
  
  # ----------------------------
  # 00 Space: Descriptives
  # ----------------------------
  source(here::here("scripts", "reporting", "00_audit.R"))
  audit(cfg)
  
  source(here::here("scripts", "reporting", "01_descriptive-participants.R"))
  descriptive_participants(cfg)
  
  source(here::here("scripts", "reporting", "02_descriptive-sequences.R"))
  descriptive_sequences(cfg)
  
  source(here::here("scripts", "reporting", "03_r-distribution.R"))
  r_distribution(cfg)
  
  source(here::here("scripts", "reporting", "04_r-switching-point-vs-bayes.R"))
  r_switching_point_vs_bayes(cfg)
  
  source(here::here("scripts", "reporting", "05_behavior-bets-and-stakes.R"))
  behavior_bets_and_stakes(cfg) 
  
  source(here::here("scripts", "reporting", "06_behavior-stakes-by-sequence.R"))
  behavior_stakes_by_sequence(cfg)
  
  source(here::here("scripts", "reporting", "07_response-times-and-optimism.R"))
  response_times_and_optimism(cfg)
  
  
  # ----------------------------
  # 10 Space: RQ1 (Betting)
  # ----------------------------
  
  source(here::here("scripts", "reporting", "11_rq1-tables.R"))
  rq1_tables(cfg)
  
  source(here::here("scripts", "reporting", "12_rq1-figures.R"))
  rq1_figures(cfg)
  

  
  # ----------------------------
  # 20 Space: RQ2 (Staking)
  # ----------------------------
  
  source(here::here("scripts", "reporting", "21_rq2-tables.R"))
  rq2_tables(cfg)
  
  source(here::here("scripts", "reporting", "22_rq2-figures.R"))
  rq2_figures(cfg)
  
  source(here::here("scripts", "reporting", "23_rq2-stake-vs-r.R"))
  rq2_stake_vs_r(cfg)
  
  source(here::here("scripts", "reporting", "24_rq2-model-comparison.R"))
  rq2_model_comparison(cfg)
  
  
  # ----------------------------
  # 30 Space: RQ3 (Welfare loss)
  # ----------------------------
  
  
  source(here::here("scripts", "reporting", "31_rq3-tables.R"))
  rq3_tables(cfg)
  
  source(here::here("scripts", "reporting", "32_rq3-figures.R"))
  rq3_figures(cfg)
  
  source(here::here("scripts", "reporting", "33_rq1-3-figures.R"))
  rq1_3_figures(cfg)
  
  # Sensitivity of RQ3 to the wealth floor (summary of the fits of 35_rq3-stan-floor-robustness.R)
  source(here::here("scripts", "reporting", "34_rq3-floor-robustness.R"))
  rq3_floor_robustness(cfg)
  
  # Comparison of the hurdle-Gamma, Gamma-only and Gaussian likelihoods
  source(here::here("scripts", "reporting", "35_rq3-likelihood-comparison.R"))
  rq3_likelihood_comparison(cfg)
  
  
  # ----------------------------
  # 40 Space: RQ4 (Side choice)
  # ----------------------------
  source(here::here("scripts", "reporting", "41_rq4-tables.R"))
  rq4_tables(cfg)
  
  source(here::here("scripts", "reporting", "42_rq4-figures.R"))
  rq4_figures(cfg)
  

  # ----------------------------
  # 50 Space: EX1 (anchor weights, pure sequences, regression summary)
  # ----------------------------
  source(here::here("scripts", "reporting", "51_ex1-weights-scatter.R"))
  ex1_weights_scatter(cfg)
  
  source(here::here("scripts", "reporting", "52_ex1-pure-sequences.R"))
  ex1_pure_sequences(cfg)
  
  source(here::here("scripts", "reporting", "53_ex1-regression-tables.R"))
  ex1_regression_tables(cfg)

  
  # ----------------------------
  # 60 Space: EX2 (optimism and response time)
  # ----------------------------
  source(here::here("scripts", "reporting", "61_ex2-tables.R"))
  ex2_tables(cfg)
  
   
  # ----------------------------
  # 70 Space: EX3 (treatment contrasts)
  # ----------------------------
  source(here::here("scripts", "reporting", "71_ex3-contrasts.R"))
  ex3_contrasts(cfg)
  
  source(here::here("scripts", "reporting", "72_ex3-figures.R"))
  ex3_figures(cfg)
  
  # ----------------------------
  # 80 Space: EX4 (canonical rules vs EX1 indices and RQ4 estimates)
  # ----------------------------
  source(here::here("scripts", "reporting", "81_ex4-similarity.R"))
  ex4_similarity(cfg)
  
  source(here::here("scripts", "reporting", "82_ex4-canonical-rules.R"))
  ex4_canonical_rules(cfg)
  
  # ----------------------------
  # 90 Space: Synthesis
  # ----------------------------
  source(here::here("scripts", "reporting", "91_synthesis-sequences.R"))
  synthesis_sequences(cfg)
  
  source(here::here("scripts", "reporting", "92_synthesis-mirror-plot.R"))
  synthesis_mirror_plot(cfg)
  
  source(here::here("scripts", "reporting", "93_synthesis-mechanism-checks.R"))
  synthesis_mechanism_checks(cfg)
  
  source(here::here("scripts", "reporting", "94_synthesis-joint-effects.R"))
  synthesis_joint_effects(cfg)
  
  source(here::here("scripts", "reporting", "95_synthesis-joint-heterogeneity.R"))
  synthesis_joint_heterogeneity(cfg)
  
  source(here::here("scripts", "reporting", "95b_synthesis-glmm-robustness.R"))
  synthesis_glmm_robustness(cfg)

  source(here::here("scripts", "reporting", "96_synthesis-debrief-validation.R"))
  synthesis_debrief_validation(cfg)

  source(here::here("scripts", "reporting", "97_synthesis-out-of-sample.R"))
  synthesis_out_of_sample(cfg)

  source(here::here("scripts", "reporting", "_synthesis-draws.R"))
  source(here::here("scripts", "reporting", "98a_synthesis-costdir-sequences.R"))
  synthesis_costdir_sequences(cfg)
  
  source(here::here("scripts", "reporting", "98b_synthesis-costdir-variance.R"))
  synthesis_costdir_variance(cfg)
  
  source(here::here("scripts", "reporting", "98c_synthesis-costdir-participants.R"))
  synthesis_costdir_participants(cfg)

  source(here::here("scripts", "reporting", "99_synthesis-reaction-times.R"))
  synthesis_reaction_times(cfg)

  # ----------------------------
  # 100 Space: Discussion
  # ----------------------------
  source(here::here("scripts", "reporting", "101_discussion-split-sample.R"))
  discussion_split_sample(cfg)
  

  msg("\nREPORTING completed.\n")
}