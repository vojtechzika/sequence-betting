############################
# ---- Indices Stage ----
############################

run_indices <- function(cfg) {
  
  stopifnot(
    is.list(cfg),
    !is.null(cfg$run),
    !is.null(cfg$design),
    !is.null(cfg$model)
  )
  
  ds     <- cfg$run$dataset

  source(here::here("scripts", "indices", "01_stan_r_from_mpl.R"))
  source(here::here("scripts", "indices", "02_score_lotr.R"))
  #source(here::here("scripts", "indices", "03_score_response_times.R")) # not used in the current pipeline - time is computed directly in the analysis
  source(here::here("scripts", "indices", "04_cache_a_star_from_r_draws.R"))
  source(here::here("scripts", "indices", "05_pids_with_positive_a_star.R"))
  source(here::here("scripts", "indices", "09_build_master_sequences.R"))
  
  
  ### Run Stuff
  
  stan_r_from_mpl(cfg)
  score_lotr(cfg)
  cache_a_star_from_r_draws(cfg) # Produces a* cache for ALL treatments defined in design$seq$treatments
  pids_with_positive_a_star(cfg) # Flags normative betters and non-betters for downstream analyses
  
  build_master_sequences(cfg)
  
  msg("\nIndices phase completed for:", ds, "\n")
  
  invisible(TRUE)
}