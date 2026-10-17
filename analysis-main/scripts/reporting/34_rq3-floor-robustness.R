# ============================================================
# 34_rq3-floor-robustness.R
#
# PURPOSE
#   Summarises the sensitivity of the RQ3 grand-mean welfare loss to the
#   numerical wealth floor x_min used at the zero-wealth state (see
#   35_rq3-stan-floor-robustness.R for the construction of the fits). For each
#   floor the script reports the posterior of the grand-mean welfare loss, the
#   share of zero losses and the convergence statistics of the fit.
#
# FUNCTION
#   rq3_floor_robustness(cfg, tr, sample)
#     tr       treatment (default: first treatment of the run)
#     sample   "confirmatory" or "full"
#
# INPUT
#   model folder rq3-floor-robustness: fit-<tr>-<sample>-xmin<floor>.rds and
#   meta-<tr>-<sample>-xmin<floor>.rds for every floor
#
# OUTPUT
#   data/output/csv/rq3-<tr>-<sample>-floor-robustness.csv
#       one row per floor: posterior mean, median and 95% interval of the
#       grand-mean welfare loss, share of zero losses, trials used,
#       divergences, maximum R-hat, minimum effective sample size
#
# NOTES
#   - The model folder is read by name; it is not selected in selected-models.csv
#   - Output files are overwritten at every run
# ============================================================

rq3_floor_robustness <- function(cfg,
                                 tr     = unique(as.character(cfg$run$treatment))[1],
                                 sample = "confirmatory") {

  stopifnot(sample %in% c("confirmatory", "full"))
  dir <- file.path(path_mod, "rq3-floor-robustness")
  if (!dir.exists(dir)) stop("Model folder rq3-floor-robustness not found; run the modelling stage.", call. = FALSE)

  metas <- list.files(dir, pattern = paste0("^meta-", tr, "-", sample, "-xmin.*\\.rds$"))
  if (length(metas) == 0L) stop("No floor fits for ", tr, "/", sample, " in ", dir, call. = FALSE)

  out <- rbindlist(lapply(metas, function(f) {
    meta <- readRDS(file.path(dir, f))
    fit  <- readRDS(file.path(dir, sub("^meta-", "fit-", f)))

    grand <- rowMeans(rstan::extract(fit, pars = "mu_c")$mu_c)
    sm_s  <- rstan::summary(fit)$summary
    n_div <- sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                        function(p) sum(p[, "divergent__"]), numeric(1)))

    cbind(as.data.table(meta),
          data.table(grand_mean = mean(grand), grand_median = median(grand),
                     grand_q025 = unname(quantile(grand, 0.025)),
                     grand_q975 = unname(quantile(grand, 0.975)),
                     n_divergences = n_div, max_rhat = max(sm_s[, "Rhat"], na.rm = TRUE),
                     min_ess = min(sm_s[, "n_eff"], na.rm = TRUE)))
  }))
  setorder(out, xmin)

  f <- file.path(path_csv, paste0("rq3-", tr, "-", sample, "-floor-robustness.csv"))
  fwrite(out, f)
  msg("Saved: ", f)

  invisible(out)
}
