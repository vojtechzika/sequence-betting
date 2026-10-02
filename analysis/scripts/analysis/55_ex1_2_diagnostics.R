# ============================================================
# 55_ex1_2_diagnostics.R
#
# PURPOSE
#   Convergence diagnostics for the EX1.2 Stan model
#   (GHI predictors: optimism, RT, risk parameter).
#
#   Checks the preregistered convergence criteria
#   (Rhat <= 1.01, ESS >= 400, zero divergent transitions)
#   for each treatment, mirroring the EX2 diagnostics.
#
# INPUT
#   path_mod/ex1_2_fit_<tr>.rds
#
# OUTPUT
#   path_out/ex1_2_diagnostics.csv
#     columns: dataset, treatment, N, Trep,
#              rhat_max, ess_min, divergences,
#              rhat_ok, ess_ok, divs_ok, converged
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex1_2_diagnostics <- function(cfg) {
  
  # ----------------------------------------------------------
  # 0. Configuration
  # ----------------------------------------------------------
  ds         <- as.character(cfg$run$data_folder)
  treatments <- unique(as.character(cfg$run$treatment))
  
  rhat_threshold <- 1.01
  ess_threshold  <- 400L
  
  f_out <- file.path(path_out, "ex1_2_diagnostics.csv")
  
  if (should_skip(f_out, cfg, "output", "EX1.2 diagnostics"))
    return(invisible(NULL))
  
  # ----------------------------------------------------------
  # 1. Loop over treatments
  # ----------------------------------------------------------
  rows <- list()
  
  for (tr in treatments) {
    
    f_fit <- file.path(path_mod, paste0("ex1_2_fit_", tr, ".rds"))
    
    if (!file.exists(f_fit)) {
      warning("EX1.2 diagnostics: fit not found for tr='", tr,
              "', skipping: ", f_fit)
      next
    }
    
    msg("EX1.2 diagnostics: checking tr='", tr, "'")
    
    fit  <- readRDS(f_fit)
    summ <- rstan::summary(fit)$summary
    
    rhat_max <- max(summ[, "Rhat"],  na.rm = TRUE)
    ess_min  <- min(summ[, "n_eff"], na.rm = TRUE)
    
    sp   <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divs <- as.integer(sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1L))))
    
    rhat_ok   <- rhat_max <= rhat_threshold
    ess_ok    <- ess_min  >= ess_threshold
    divs_ok   <- divs == 0L
    converged <- rhat_ok & ess_ok & divs_ok
    
    # Recover N and Trep from y_rep dimensions (matrix[Trep, N] in data block)
    y_rep_dims <- fit@sim$dims_oi[["y_rep"]]
    Trep <- if (!is.null(y_rep_dims) && length(y_rep_dims) >= 1L)
      as.integer(y_rep_dims[1L]) else NA_integer_
    N    <- if (!is.null(y_rep_dims) && length(y_rep_dims) >= 2L)
      as.integer(y_rep_dims[2L]) else NA_integer_
    
    msg("  Rhat_max=",    round(rhat_max, 4),
        " | ESS_min=",    round(ess_min,  0),
        " | divergences=", divs,
        " | converged=",  converged)
    
    if (!converged)
      warning("EX1.2 diagnostics: convergence criteria not met for tr='", tr, "'.")
    
    rows[[tr]] <- data.table(
      dataset     = ds,
      treatment   = tr,
      N           = N,
      Trep        = Trep,
      rhat_max    = round(rhat_max, 4),
      ess_min     = round(ess_min,  0),
      divergences = as.integer(divs),
      rhat_ok     = rhat_ok,
      ess_ok      = ess_ok,
      divs_ok     = divs_ok,
      converged   = converged
    )
  }
  
  # ----------------------------------------------------------
  # 2. Save
  # ----------------------------------------------------------
  if (length(rows) == 0L) {
    warning("EX1.2 diagnostics: no fits found; output not written.")
    return(invisible(NULL))
  }
  
  out <- rbindlist(rows)
  fwrite(out, f_out)
  msg("Saved: ", f_out)
  
  invisible(out)
}