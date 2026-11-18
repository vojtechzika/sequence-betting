# ============================================================
# 25_rq2-floor-robustness.R
#
# PURPOSE
#   Summarises the robustness of the sequence-level stake deviation mu_a of
#   the Gaussian RQ2 model to the standardisation of the stake deviations
#   (preregistered robustness checks): the floor of the participant standard
#   deviation set to 1 and to 5 in place of the main value (cfg$design$rq2$sd_floor),
#   and the median absolute deviation (MAD) in place of the standard deviation
#   (fits built by 21_rq2-stan-gaussian.R). Each variant is compared with the
#   main fit of the same sample.
#
# FUNCTION
#   rq2_floor_robustness(cfg, tr, sample)
#     tr       treatment (default: first treatment of the run)
#     sample   sample of the variant fits (default "confirmatory")
#
# INPUT
#   model folder rq2-gaussian: fit-<tr>-<sample>.rds (main) and
#   fit-<tr>-<sample>-<variant>.rds, seq-levels-<tr>-<sample>[-<variant>].rds
#   for the variants floor1, floor5 and mad
#
# OUTPUT
#   data/output/csv/rq2-<tr>-<sample>-floor-robustness.csv
#       one row per variant (the main fit included): grand mean of mu_a,
#       largest absolute change of the posterior mean of mu_a relative to the
#       main fit, in posterior standard deviations of the main fit and as a
#       share of the endowment, Spearman correlation of the sequence means
#       with the main fit, and the numbers of sequences with posterior
#       probability of at least 0.95 of under- and of over-staking beyond
#       the main threshold rho (cfg$design$rq2$rho[1])
#
# NOTES
#   - The model folder is read by name; it is not selected in selected-models.csv
#   - Output files are overwritten at every run
# ============================================================

rq2_floor_robustness <- function(cfg,
                                 tr     = unique(as.character(cfg$run$treatment))[1],
                                 sample = "confirmatory") {

  folder   <- "rq2-gaussian"
  variants <- c("main", "floor1", "floor5", "mad")
  rho      <- as.numeric(cfg$design$rq2$rho)[1]

  # ---- Posterior draws of mu_a (draws x sequences), sequences in a common order ----
  read_mu <- function(variant) {
    sfx <- if (variant == "main") paste0(tr, "-", sample) else paste0(tr, "-", sample, "-", variant)
    f_fit <- file.path(path_mod, folder, paste0("fit-", sfx, ".rds"))
    if (!file.exists(f_fit)) return(NULL)
    lev <- as.character(readRDS(file.path(path_mod, folder, paste0("seq-levels-", sfx, ".rds"))))
    mu  <- rstan::extract(readRDS(f_fit), pars = "mu_a")$mu_a
    stopifnot(is.matrix(mu), ncol(mu) == length(lev))
    colnames(mu) <- lev
    mu
  }

  draws <- lapply(variants, read_mu)
  names(draws) <- variants
  if (is.null(draws$main)) stop("Main Gaussian fit for ", tr, "/", sample, " not found in ", folder, call. = FALSE)
  draws <- Filter(Negate(is.null), draws)

  ref   <- draws$main
  m_ref <- colMeans(ref)
  sd_ref <- apply(ref, 2, sd)

  out <- rbindlist(lapply(names(draws), function(v) {
    mu <- draws[[v]][, colnames(ref), drop = FALSE]
    m  <- colMeans(mu)
    U  <- colMeans(mu < -rho)
    O  <- colMeans(mu >  rho)
    data.table(
      variant = v, treatment = tr, sample = sample,
      grand_mean = mean(rowMeans(mu)),
      max_shift_sd = max(abs(m - m_ref) / sd_ref),
      max_shift_endowment = max(abs(m - m_ref)),
      spearman_main = cor(m, m_ref, method = "spearman"),
      n_under_strong = sum(U >= 0.95), n_over_strong = sum(O >= 0.95),
      rho = rho
    )
  }))

  f <- file.path(path_csv, paste0("rq2-", tr, "-", sample, "-floor-robustness.csv"))
  fwrite(out, f)
  msg("Saved: ", f)

  invisible(out)
}
