# ============================================================
# 03_cache_a_star_from_r_draws.R
#
# PURPOSE
#   Computes the expected-utility optimal stake a*(r_i; m)
#   using posterior draws of risk preferences from the
#   Holt-Laury model. Saves a matrix of a* draws per
#   treatment for use in RQ2 and RQ3.
#
# INPUT
#   path_mod/mpl_r_draws_<tr>.rds  -- posterior draws of r_i
#
# OUTPUT
#   path_mod/a_star_draws_<tr>.rds -- matrix of a* draws (iters x N)
#
# NOTES
#   - Computed separately per treatment (no pooling)
#   - Uses discrete maximization over a in {0,...,endowment}
#   - If consistent_only = TRUE, uses consistent-only r draws
#   - r_tol and a_star_discrete defined in 03_eu_functions.R
# ============================================================

cache_a_star_from_r_draws <- function(cfg) {
  
  design <- cfg$design
  
  e     <- design$seq$endowment
  xmin  <- design$seq$xmin
  p_win <- design$seq$coin_prob
  grid  <- 0:e
  r_tol <- get("R_TOL_DEFAULT", inherits = TRUE)
  
  for (tr in names(design$seq$treatments)) {
    
    # Use consistent-only draws if requested
    tag <- if (isTRUE(cfg$run$consistent_only)) {
      paste0(tr, "_consistent")
    } else {
      tr
    }
    
    f_r   <- file.path(path_mod, paste0("mpl_r_draws_", tag, ".rds"))
    f_out <- file.path(path_mod, paste0("a_star_draws_", tr, ".rds"))
    
    stopifnot(file.exists(f_r))
    
    if (should_skip(
      paths = f_out,
      cfg   = cfg,
      type  = "model",
      label = paste0("a_star cache (", tr, ")")
    )) next
    
    hl         <- readRDS(f_r)
    pid_levels <- hl$pid
    r_draws    <- hl$r_draws
    iters      <- nrow(r_draws)
    N          <- ncol(r_draws)
    m          <- design$seq$treatments[[tr]]
    
    a_star_draws <- matrix(NA_integer_, nrow = iters, ncol = N)
    
    for (i in seq_len(N)) {
      a_star_draws[, i] <- vapply(
        r_draws[, i],
        FUN       = a_star_discrete,
        FUN.VALUE = integer(1),
        m         = m,
        e         = e,
        p_win     = p_win,
        xmin      = xmin,
        r_tol     = r_tol,
        grid      = grid
      )
    }
    
    saveRDS(
      list(pid = pid_levels, a_star_draws = a_star_draws),
      f_out
    )
    
    msg("Saved: ", f_out,
        " | tag: ", tag,
        " | iters: ", iters,
        " | N: ", N)
  }
}