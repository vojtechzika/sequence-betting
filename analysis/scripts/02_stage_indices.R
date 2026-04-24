############################
# ---- Indices Stage ----
############################

run_indices <- function(cfg) {
  
  source(here::here("scripts", "indices", "01_score_lotr.R"))
  source(here::here("scripts", "indices", "02_stan_r_from_mpl.R"))
  source(here::here("scripts", "indices", "03_cache_a_star_from_r_draws.R"))
  source(here::here("scripts", "indices", "04_pids_with_positive_a_star.R"))
  source(here::here("scripts", "indices", "05_classic_r.R"))
  source(here::here("scripts", "indices", "06_block_drift.R"))
  source(here::here("scripts", "indices", "07_button_order.R"))

  ### Run Stuff
  score_lotr(cfg)
  stan_r_from_mpl(cfg)
  cache_a_star_from_r_draws(cfg) # Produces a* cache for ALL treatments defined in design$seq$treatments
  pids_with_positive_a_star(cfg) # Flags normative betters and non-betters for downstream analyses
  classic_r(cfg) # r estimation as a sanity check
  block_drift_check(cfg)
  button_order_check(cfg)
  
  msg("\nINDICES completed for data_folder:", cfg$run$data_folder, "\n")
  
  invisible(TRUE)
}