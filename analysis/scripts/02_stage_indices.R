############################
# ---- Indices Stage ----
############################

run_indices <- function(cfg) {
  
  source(here::here("scripts", "indices", "01_build_master_sequences.R"))
  source(here::here("scripts", "indices", "02_score_lotr.R"))
  source(here::here("scripts", "indices", "03_stan_r_from_mpl.R"))
  source(here::here("scripts", "indices", "04_cache_a_star_from_r_draws.R"))
  source(here::here("scripts", "indices", "05_pids_with_positive_a_star.R"))
  
  
  ### Run Stuff
  build_master_sequences(cfg)
  score_lotr(cfg)
  stan_r_from_mpl(cfg)
  cache_a_star_from_r_draws(cfg) # Produces a* cache for ALL treatments defined in design$seq$treatments
  pids_with_positive_a_star(cfg) # Flags normative betters and non-betters for downstream analyses
  
  
  
  msg("\nIndices completed for:", cfg$run$dataset, "\n")
  
  invisible(TRUE)
}