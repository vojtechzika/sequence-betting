# ============================================================
# 34_rq3-floor-robustness.R
#
# PURPOSE
#   Summarises the sensitivity of the RQ3 welfare loss to the numerical
#   wealth floor x_min used at the zero-wealth state (see
#   35_rq3-stan-floor-robustness.R for the construction of the fits). For each
#   floor the script reports the posterior of the grand-mean welfare loss, the
#   share of zero losses and the convergence statistics of the fit, and it
#   compares the sequence-level estimates mu_c of every floor with those of the
#   baseline floor (cfg$design$seq$xmin), which is refitted by the same code.
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
#       divergences, maximum R-hat, minimum effective sample size; relative to
#       the baseline floor, the largest absolute change of the posterior mean
#       of mu_c over sequences (in posterior standard deviations of the
#       baseline fit and as a share of the endowment), the Spearman correlation
#       of the sequence means, the number of sequences classified as strong
#       (posterior probability of a loss above rho of at least 0.95, rho =
#       cfg$design$rq3$rho[1]) and the number of sequences whose class
#       (strong, moderate, weak, neutral; cut-offs of 31_rq3-tables.R) differs
#       from the baseline
#
# NOTES
#   - The model folder is read by name; it is not selected in selected-models.csv
#   - Fits of one sample share the order of the sequences; the number of
#     sequences is checked
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

  rho <- as.numeric(cfg$design$rq3$rho)[1]

  # ---- Sequence classes from the probability of a loss above rho (cut-offs of 31_rq3-tables.R) ----
  seq_class <- function(L)
    cut(L, breaks = c(-Inf, 0.50, 0.80, 0.95, Inf),
        labels = c("neutral", "weak", "moderate", "strong"), right = FALSE)

  # ---- Fits, metadata and draws of mu_c (draws x sequences) for every floor ----
  fits <- lapply(metas, function(f) {
    fit <- readRDS(file.path(dir, sub("^meta-", "fit-", f)))
    list(meta = readRDS(file.path(dir, f)), fit = fit,
         mu = rstan::extract(fit, pars = "mu_c")$mu_c)
  })
  fits <- fits[order(vapply(fits, function(x) x$meta$xmin, numeric(1)))]

  is_base <- vapply(fits, function(x) isTRUE(x$meta$is_baseline), logical(1))
  if (sum(is_base) != 1L) stop("Exactly one baseline floor is expected for ", tr, "/", sample, call. = FALSE)
  ref     <- fits[[which(is_base)]]$mu
  m_ref   <- colMeans(ref)
  sd_ref  <- apply(ref, 2, sd)
  cls_ref <- seq_class(colMeans(ref > rho))

  out <- rbindlist(lapply(fits, function(x) {
    mu <- x$mu
    stopifnot(is.matrix(mu), ncol(mu) == ncol(ref))
    grand <- rowMeans(mu)
    m     <- colMeans(mu)
    cls   <- seq_class(colMeans(mu > rho))
    sm_s  <- rstan::summary(x$fit)$summary
    n_div <- sum(vapply(rstan::get_sampler_params(x$fit, inc_warmup = FALSE),
                        function(p) sum(p[, "divergent__"]), numeric(1)))

    cbind(as.data.table(x$meta),
          data.table(grand_mean = mean(grand), grand_median = median(grand),
                     grand_q025 = unname(quantile(grand, 0.025)),
                     grand_q975 = unname(quantile(grand, 0.975)),
                     n_divergences = n_div, max_rhat = max(sm_s[, "Rhat"], na.rm = TRUE),
                     min_ess = min(sm_s[, "n_eff"], na.rm = TRUE),
                     max_shift_sd = max(abs(m - m_ref) / sd_ref),
                     max_shift_endowment = max(abs(m - m_ref)),
                     spearman_baseline = cor(m, m_ref, method = "spearman"),
                     n_strong = sum(cls == "strong"),
                     n_label_changed = sum(cls != cls_ref),
                     rho = rho))
  }))
  setorder(out, xmin)

  f <- file.path(path_csv, paste0("rq3-", tr, "-", sample, "-floor-robustness.csv"))
  fwrite(out, f)
  msg("Saved: ", f)

  invisible(out)
}
