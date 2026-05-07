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
#   path_mod/mpl_r_draws_<tr>_consistent.rds  -- consistent-only r draws
#                                                (primary when consistent_only = TRUE)
#   path_mod/mpl_r_draws_<tr>.rds             -- full-sample r draws
#                                                (primary when consistent_only = FALSE;
#                                                 also used for _full when consistent_only = TRUE)
#
# OUTPUT
#   path_mod/a_star_draws_<tr>.rds       -- primary (consistent-only or full depending on flag)
#   path_mod/a_star_draws_<tr>_full.rds  -- full-sample descriptive,
#                                           only produced if consistent_only = TRUE
#
# NOTES
#   - Computed separately per treatment (no pooling)
#   - Uses discrete maximization over a in {0,...,endowment}
#   - r_tol and a_star_discrete defined in 03_eu_functions.R
# ============================================================
cache_a_star_from_r_draws <- function(cfg) {
  
  design         <- cfg$design
  consistent_only <- isTRUE(cfg$run$consistent_only)
  
  e     <- design$seq$endowment
  xmin  <- design$seq$xmin
  p_win <- design$seq$coin_prob
  grid  <- 0:e
  r_tol <- get("R_TOL_DEFAULT", inherits = TRUE)
  
  for (tr in names(design$seq$treatments)) {
    
    # Build tag_map: r-draws input tag -> a_star output stem
    tag_map <- list()
    
    f_r_consistent <- file.path(path_mod, paste0("mpl_r_draws_", tr, "_consistent.rds"))
    f_r_full       <- file.path(path_mod, paste0("mpl_r_draws_", tr, ".rds"))
    
    if (consistent_only && file.exists(f_r_consistent)) {
      tag_map[[paste0(tr, "_consistent")]] <- paste0("a_star_draws_", tr)
    } else if (file.exists(f_r_full)) {
      tag_map[[tr]] <- paste0("a_star_draws_", tr)
    }
    
    if (consistent_only && file.exists(f_r_full)) {
      tag_map[[tr]] <- paste0("a_star_draws_", tr, "_full")
    }
    
    for (r_tag in names(tag_map)) {
      
      f_r   <- file.path(path_mod, paste0("mpl_r_draws_", r_tag, ".rds"))
      f_out <- file.path(path_mod, paste0(tag_map[[r_tag]], ".rds"))
      
      stopifnot(file.exists(f_r))
      
      if (should_skip(
        paths = f_out,
        cfg   = cfg,
        type  = "model",
        label = paste0("a_star cache (", r_tag, ")")
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
          " | tag: ", r_tag,
          " | iters: ", iters,
          " | N: ", N)
    }
  }
}