# ============================================================
# 09_mpl_diagnostics.R
#
# PURPOSE
#   Reads existing MPL Stan fit objects and produces convergence
#   diagnostics CSV. Runs independently of 02_stan_r_from_mpl.R
#   so diagnostics can be regenerated without refitting models.
#
# INPUT
#   path_mod/mpl_fit_<tag>.rds  -- Stan fit objects
#
# OUTPUT
#   path_out/mpl_diagnostics.csv
# ============================================================

mpl_diagnostics <- function(cfg) {
  
  design <- cfg$design
  consistent_only <- isTRUE(cfg$run$consistent_only)
  
  all_tags <- c("pooled", names(design$seq$treatments))
  if (consistent_only) {
    all_tags <- c(all_tags, "pooled_consistent",
                  paste0(names(design$seq$treatments), "_consistent"))
  }
  
  f_out <- file.path(path_out, "mpl_diagnostics.csv")
  if (should_skip(f_out, cfg, "output", "MPL diagnostics")) return(invisible(NULL))
  
  rows <- list()
  
  for (tag in all_tags) {
    
    f_fit <- file.path(path_mod, paste0("mpl_fit_", tag, ".rds"))
    if (!file.exists(f_fit)) {
      msg("MPL diagnostics: fit not found for tag=", tag, " -- skipping")
      next
    }
    
    fit  <- readRDS(f_fit)
    summ <- summary(fit)$summary
    sp   <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    
    rhat_max    <- max(summ[, "Rhat"], na.rm = TRUE)
    ess_min     <- min(summ[, "n_eff"], na.rm = TRUE)
    divergences <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    
    msg("MPL diagnostics (", tag, "):",
        " Rhat_max=", round(rhat_max, 4),
        " ESS_min=",  round(ess_min, 0),
        " divergences=", divergences)
    
    rows[[tag]] <- data.table(
      tag         = tag,
      rhat_max    = round(rhat_max, 4),
      ess_min     = round(ess_min, 0),
      divergences = divergences
    )
  }
  
  if (length(rows) > 0) {
    fwrite(rbindlist(rows), f_out)
    msg("Saved: ", f_out)
  }
  
  invisible(rows)
}