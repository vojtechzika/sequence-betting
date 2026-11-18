############################
# ---- Reporting Stage ----
############################

run_reporting <- function(cfg) {
  
  # ----------------------------
  # Helpers
  # ----------------------------
  source(here::here("scripts", "reporting", "_reporting-style.R"))
  
  # Sample helpers (analysis_trials) are defined in the modelling stage; load them
  # when reporting is run in a session where the modelling stage has not been run.
  if (!exists("analysis_trials", mode = "function"))
    source(here::here("scripts", "modelling", "09_analysis-sample.R"))
  
  
  # ----------------------------
  # 00 Space: Descriptives
  # ----------------------------
  source(here::here("scripts", "reporting", "00_audit.R"))
  audit(cfg)
  
  source(here::here("scripts", "reporting", "01_descriptive-participants.R"))
  descriptive_participants(cfg)
  
  source(here::here("scripts", "reporting", "02_descriptive-sequences.R"))
  descriptive_sequences(cfg)
  
  source(here::here("scripts", "reporting", "03_r-descriptive.R"))
  r_descriptive(cfg)
  
  source(here::here("scripts", "reporting", "04_r-distribution.R"))
  r_distribution(cfg)
  
  source(here::here("scripts", "reporting", "05_r-robustness.R"))
  r_robustness(cfg)
  
  source(here::here("scripts", "reporting", "06_betting-task-distribution.R"))
  betting_task_distribution(cfg)
  
  source(here::here("scripts", "reporting", "07_behavior-stakes-by-sequence.R"))
  behavior_stakes_by_sequence(cfg)
  
  source(here::here("scripts", "reporting", "08_response-time-distribution.R"))
  response_time_distribution(cfg)
  
  source(here::here("scripts", "reporting", "09_optimism-distribution.R"))
  optimism_distribution(cfg)
  
  
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
  
  source(here::here("scripts", "reporting", "25_rq2-floor-robustness.R"))
  rq2_floor_robustness(cfg)
  
  
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

  source(here::here("scripts", "reporting", "43_rq4-top-sequences-table.R"))
  rq4_top_sequences_table(cfg)
  

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

  source(here::here("scripts", "reporting", "62_ex2-figures.R"))
  ex2_figures(cfg)
  
   
  # ----------------------------
  # 70 Space: EX3 (treatment contrasts)
  # ----------------------------
  source(here::here("scripts", "reporting", "71_ex3-contrasts.R"))
  ex3_contrasts(cfg)
  
  # 72_ was retired: the figure is now produced by 33_rq1-3-figures.R
  
  source(here::here("scripts", "reporting", "73_ex3-summary.R"))
  ex3_summary(cfg)
  
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

  source(here::here("scripts", "reporting", "95c_synthesis-run-share-by-imbalance.R"))
  synthesis_run_share_by_imbalance(cfg)

  source(here::here("scripts", "reporting", "96_synthesis-debrief-validation.R"))
  synthesis_debrief_validation(cfg)
  source(here::here("scripts", "reporting", "96b_synthesis-debrief-descriptives.R"))
  synthesis_debrief_descriptives(cfg)

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
  # prior sensitivity (skips itself when no broad fits were requested)
  source(here::here("scripts", "reporting", "102_prior-sensitivity.R"))
  prior_sensitivity(cfg)

  # ----------------------------
  # 200 Space: tables of the paper, in order of appearance
  # (main text first, then the online appendix)
  # ----------------------------
  # descriptive statistics by treatment
  source(here::here("scripts", "reporting", "201_descriptive-table.R"))
  descriptive_table(cfg)

  source(here::here("scripts", "reporting", "202_model-summary-table.R"))
  model_summary_table(cfg)
  # joint model of the Synthesis: effects and cross-equation correlations
  source(here::here("scripts", "reporting", "203_synthesis-joint-table.R"))
  synthesis_joint_table(cfg)
  # individual heterogeneity of the Synthesis: participant-level heterogeneity
  source(here::here("scripts", "reporting", "204_synthesis-heterogeneity-table.R"))
  synthesis_heterogeneity_table(cfg)
  # personal run slope by imbalance class
  source(here::here("scripts", "reporting", "205_synthesis-run-by-imb-table.R"))
  run_by_imb_table(cfg)
  # mechanism checks of the Synthesis
  source(here::here("scripts", "reporting", "206_mechanism-checks-table.R"))
  mechanism_checks_table(cfg)
  # sequence-level estimates of RQ1 to RQ4
  source(here::here("scripts", "reporting", "207_sequence-estimates-table.R"))
  sequence_estimates_table(cfg)
  # deviations of the sequence-level estimates from the grand means
  source(here::here("scripts", "reporting", "208_sequence-deviations-table.R"))
  sequence_deviations_table(cfg)
  # top-5 groups of the deviations (not used in the current drafts)
  source(here::here("scripts", "reporting", "209_synthesis-top-sequences-table.R"))
  synthesis_top_sequences_table(cfg)
  # share of bets against the run by imbalance class and level
  source(here::here("scripts", "reporting", "210_synthesis-run-share-table.R"))
  synthesis_run_share_table(cfg)
  # closed debriefing items by slope class
  source(here::here("scripts", "reporting", "211_synthesis-debrief-likert-table.R"))
  synthesis_debrief_likert_table(cfg)
  # descriptive statistics of the debriefing questions
  source(here::here("scripts", "reporting", "212_synthesis-debrief-descriptives-table.R"))
  synthesis_debrief_descriptives_table(cfg)

  msg("\nREPORTING completed.\n")
}