# scripts/indices/04_cache_a_star_from_r_draws.R

library(data.table)

cache_a_star_from_r_draws <- function(cfg) {
  
  stopifnot(
    is.list(cfg),
    !is.null(cfg$run),
    !is.null(cfg$design),
    !is.null(cfg$model),
    !is.null(cfg$run$dataset)
  )
  
  ds <- as.character(cfg$run$dataset)
  design  <- cfg$design
  model   <- cfg$model
  
  # ----------------------------
  # Inputs
  # ----------------------------
  f_draws <- file.path(path_mod_ds(ds), "mpl_r_draws.rds")
  stopifnot(file.exists(f_draws))
  
  hl <- readRDS(f_draws)
  stopifnot(is.list(hl), all(c("pid", "r_draws") %in% names(hl)))
  
  pid_levels <- as.character(hl$pid)
  r_draws <- hl$r_draws  # iters x N
  stopifnot(is.matrix(r_draws))
  
  iters <- nrow(r_draws)
  N <- ncol(r_draws)
  stopifnot(length(pid_levels) == N)
  
  # ----------------------------
  # Design + numeric params
  # ----------------------------
  stopifnot(!is.null(design$seq))
  e      <- as.integer(design$seq$endowment)
  xmin   <- as.numeric(design$seq$xmin)
  p_win  <- if (!is.null(design$seq$coin_prob)) as.numeric(design$seq$coin_prob) else 0.5
  
  stopifnot(length(e) == 1L, e > 0L)
  stopifnot(length(xmin) == 1L, is.finite(xmin), xmin > 0)
  stopifnot(length(p_win) == 1L, is.finite(p_win), p_win > 0, p_win < 1)
  
  stopifnot(!is.null(design$seq$treatments), length(design$seq$treatments) > 0)
  tr_names <- names(design$seq$treatments)
  stopifnot(length(tr_names) > 0, all(nzchar(tr_names)))
  
  # stake grid (integer ECU)
  grid <- 0:e
  
  # tolerance (constant; used in a_star_discrete)
  r_tol <- get("R_TOL_DEFAULT", inherits = TRUE)
  stopifnot(length(r_tol) == 1L, is.finite(r_tol), r_tol > 0)
  
  # ----------------------------
  # Compute + save per treatment
  # ----------------------------
  out_files <- character(0)
  
  for (tr in tr_names) {
    
    m <- design$seq$treatments[[tr]]
    if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m <= 1) {
      stop("Invalid multiplier for treatment '", tr, "': ", paste(m, collapse = ", "))
    }
    
    msg("Caching a* for treatment:", tr, "| multiplier:", m)
    
    a_star_draws <- matrix(NA_integer_, nrow = iters, ncol = N)
    for (i in 1:N) {
      a_star_draws[, i] <- vapply(
        r_draws[, i],
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
    
    out <- list(
      pid = pid_levels,
      a_star_draws = a_star_draws,
      params = list(
        dataset = ds,
        treatment = tr,
        mult = m,
        endowment = e,
        coin_prob = p_win,
        xmin = xmin,
        r_tol = r_tol
      )
    )
    
    f_out <- file.path(path_mod_ds(ds), paste0("a_star_draws_", tr, ".rds"))
    saveRDS(out, f_out)
    msg("Saved:", f_out)
    
    out_files <- c(out_files, f_out)
  }
  
  invisible(out_files)
}

# Example:
# run    <- run_cfg()
# design <- design_cfg()
# model  <- model_cfg()
# cache_a_star_from_r_draws(run, design, model)