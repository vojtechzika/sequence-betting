# ============================================================
# scripts/indices/04_cache_a_star_from_r_draws.R
#
# Computes the expected-utility optimal stake a*(r_i; m)
# using posterior draws of risk preferences from the
# Holt–Laury model.
#
# For each treatment in design$seq$treatments:
#   - Reads mpl_r_draws_<tr>.rds
#   - Computes a* for every posterior draw and participant
#   - Saves a_star_draws_<tr>.rds
#
# The stake is obtained via discrete maximization over
# a ∈ {0,1,...,endowment} using CRRA utility.
#
# This script is treatment-consistent:
#   a_star_draws_<tr> are computed ONLY from the
#   corresponding mpl_r_draws_<tr>.
#
# No pooling across treatments occurs here.
# ============================================================

library(data.table)

cache_a_star_from_r_draws <- function(cfg) {
  
  stopifnot(
    is.list(cfg),
    !is.null(cfg$run),
    !is.null(cfg$design),
    !is.null(cfg$run$dataset)
  )
  
  ds     <- as.character(cfg$run$dataset)
  design <- cfg$design
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  e     <- design$seq$endowment
  xmin  <- design$seq$xmin
  p_win <- design$seq$coin_prob
  
  grid  <- 0:e
  r_tol <- get("R_TOL_DEFAULT", inherits = TRUE)
  
  # ============================================================
  # For each treatment: read matching r_draws_<tr>.rds
  # ============================================================
  
  for (tr in names(design$seq$treatments)) {
    
    f_r <- file.path(mod_dir, paste0("mpl_r_draws_", tr, ".rds"))
    stopifnot(file.exists(f_r))
    
    f_out <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    
    if (should_skip(f_out, cfg, "model",
                    paste0("a_star cache (", ds, "/", tr, ")"))) next
    
    hl <- readRDS(f_r)
    pid_levels <- hl$pid
    r_draws    <- hl$r_draws
    
    iters <- nrow(r_draws)
    N     <- ncol(r_draws)
    
    m <- design$seq$treatments[[tr]]
    
    a_star_draws <- matrix(NA_integer_, nrow=iters, ncol=N)
    
    for (i in 1:N) {
      a_star_draws[,i] <- vapply(
        r_draws[,i],
        FUN = a_star_discrete,
        FUN.VALUE = integer(1),
        m = m,
        e = e,
        p_win = p_win,
        xmin = xmin,
        r_tol = r_tol,
        grid = grid
      )
    }
    
    saveRDS(
      list(pid = pid_levels,
           a_star_draws = a_star_draws),
      f_out
    )
  }
}