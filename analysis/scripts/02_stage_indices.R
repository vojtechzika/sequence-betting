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
  design <- cfg$design
  model  <- cfg$model
  
  source(here::here("scripts", "indices", "01_stan_r_from_mpl.R"))
  source(here::here("scripts", "indices", "02_score_lotr.R"))
  source(here::here("scripts", "indices", "03_build_master_sequences.R"))
  source(here::here("scripts", "indices", "04_cache_a_star_from_r_draws.R"))
  
  stan_r_from_mpl(run = cfg$run, design = design, model = model)
  score_lotr(ds)
  build_master_sequences(ds)
  
  # Produces a* cache for ALL treatments defined in design$seq$treatments
  cache_a_star_from_r_draws(run = cfg$run, design = design, model = model)
  
  msg("\nIndices phase completed for:", ds, "\n")
  
  invisible(TRUE)
}