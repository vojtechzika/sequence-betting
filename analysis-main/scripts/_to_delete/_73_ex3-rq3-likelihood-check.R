# ============================================================
# 73_ex3-rq3-likelihood-check.R
#
# PURPOSE
#   Checks that, in the full FN sample, the hurdle-Gamma and Gamma-only RQ3
#   models give essentially the same welfare-loss estimates. This justifies
#   using the hurdle-Gamma model in both treatments for the EX3 contrast
#   (71_ex3-contrasts.R), where the FP arm requires the hurdle component.
#   The script reports the grand mean welfare loss of each model with its
#   95% credible interval, the Spearman correlation of the sequence means,
#   and the largest absolute difference between the sequence means.
#
# SAMPLE
#   The full sample of the first treatment in cfg$run$treatment (FN).
#
# INPUT
#   models "rq3-hurdle-gamma" and "rq3-gamma": fit-<tr>-full.rds (mu_c)
#
# OUTPUT (data/output/csv/)
#   ex3-rq3-likelihood-check.csv
#
# NOTES
#   - The two models are read by folder name, independent of the model
#     selected for the key "rq3"
#   - Output files are overwritten at every run
# ============================================================

ex3_rq3_likelihood_check <- function(cfg) {

  tr <- unique(as.character(cfg$run$treatment))[1L]
  get_mu <- function(folder)
    rstan::extract(readRDS(model_path(folder, paste0("fit-", tr, "-full.rds"))), pars = "mu_c")$mu_c

  h <- get_mu("rq3-hurdle-gamma")
  g <- get_mu("rq3-gamma")
  stopifnot(is.matrix(h), is.matrix(g), ncol(h) == ncol(g))

  gh <- rowMeans(h); gg <- rowMeans(g)
  qs <- function(x) quantile(x, c(0.025, 0.975), names = FALSE)
  mh <- colMeans(h); mg <- colMeans(g)

  out <- data.table(
    treatment = tr, sample = "full",
    grand_hurdle_median = median(gh), grand_hurdle_q025 = qs(gh)[1], grand_hurdle_q975 = qs(gh)[2],
    grand_gamma_median  = median(gg), grand_gamma_q025  = qs(gg)[1], grand_gamma_q975  = qs(gg)[2],
    spearman_seq = cor(mh, mg, method = "spearman"),
    max_abs_diff_seq = max(abs(mh - mg))
  )
  f <- file.path(path_csv, "ex3-rq3-likelihood-check.csv")
  fwrite(out, f);  msg("Saved: ", f)

  invisible(TRUE)
}
