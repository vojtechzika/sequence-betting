# ============================================================
# 73_ex3-summary.R
#
# PURPOSE
#   EX3: one-row summary per research question (betting, RQ1; stake deviation,
#   RQ2, fnconf-fnfull only; welfare loss, RQ3) and comparison of the contrasts computed in 71_ex3-contrasts.R:
#   fnfull-fpfull and fnconf-fnfull (first sample minus second
#   sample). For each quantity and comparison the summary gives the grand mean
#   of each sample, the contrast of the averages of the participant means, the
#   number of sequences whose contrast has a 95% credible interval above or
#   below zero, the number of sequences whose contrast exceeds the
#   practical-equivalence threshold with a posterior probability of at least
#   0.95, the same count and the posterior probability for the additional
#   preregistered threshold rho_alt (NA for the stake deviation), the range of
#   the sequence contrasts, and the rank correlation of the
#   sequence means across the two samples. Nothing is recomputed from the fits.
#
# INPUT (data/output/csv/)
#   ex3-rq<k>-sequence-contrasts-<tag>.csv      k = 1, 2, 3
#   ex3-rq<k>-participant-contrasts-<tag>.csv
#   ex3-rq<k>-spearman-<tag>.csv
#   rq1-<tr>-<sample>-model-summary.csv         grand_mean_bet
#   rq2-<tr>-<sample>-model-summary.csv         grand_mean_a
#   rq3-<tr>-<sample>-model-summary.csv         grand_mean_c
#
# OUTPUT
#   data/output/csv/ex3-summary.csv       one row per research question and comparison
#
# NOTES
#   - A sequence is counted as higher in the first (second) sample if the 95%
#     credible interval of its contrast lies above (below) zero
#   - The threshold is the practical-equivalence threshold rho of the
#     research question used in 71_ex3-contrasts.R
#   - Output files are overwritten at every run
# ============================================================

ex3_summary <- function(cfg) {

  tags <- c("fnfull-fpfull", "fnconf-fnfull")
  specs <- list(
    list(rq = "rq1", par = "grand_mean_bet", tags = tags),
    list(rq = "rq2", par = "grand_mean_a",   tags = "fnconf-fnfull"),
    list(rq = "rq3", par = "grand_mean_c",   tags = tags)
  )

  read_out <- function(f) {
    p <- file.path(path_csv, f)
    stopifnot(file.exists(p))
    fread(p, encoding = "UTF-8")
  }

  summ <- rbindlist(lapply(specs, function(sp) {
    rbindlist(lapply(sp$tags, function(tag) {
      sq <- read_out(paste0("ex3-", sp$rq, "-sequence-contrasts-", tag, ".csv"))
      pc <- read_out(paste0("ex3-", sp$rq, "-participant-contrasts-", tag, ".csv"))
      sm <- read_out(paste0("ex3-", sp$rq, "-spearman-", tag, ".csv"))
      gm <- function(tr, smp) {
        m <- read_out(paste0(sp$rq, "-", tr, "-", smp, "-model-summary.csv"))[parameter == sp$par]
        stopifnot(nrow(m) == 1L)
        m
      }
      g_a <- gm(pc$treatment_a, pc$sample_a); g_b <- gm(pc$treatment_b, pc$sample_b)

      data.table(
        rq                  = sp$rq, contrast = tag,
        n_participants_a    = pc$n_pid_a, n_participants_b = pc$n_pid_b,
        grand_mean_a        = g_a$median, grand_mean_a_q025 = g_a$q025, grand_mean_a_q975 = g_a$q975,
        grand_mean_b        = g_b$median, grand_mean_b_q025 = g_b$q025, grand_mean_b_q975 = g_b$q975,
        contrast_median     = pc$median, contrast_q025 = pc$q025, contrast_q975 = pc$q975,
        rho                 = pc$rho, p_abs_gt_rho = pc$p_abs_gt_rho,
        rho_alt             = pc$rho_alt, p_abs_gt_rho_alt = pc$p_abs_gt_rho_alt,
        n_sequences         = nrow(sq),
        n_seq_higher_a      = sum(sq$q025 > 0),
        n_seq_higher_b      = sum(sq$q975 < 0),
        n_seq_beyond_rho    = sum(sq$p_abs_gt_rho >= 0.95),
        n_seq_beyond_rho_alt = if (is.na(pc$rho_alt)) NA_integer_ else sum(sq$p_abs_gt_rho_alt >= 0.95),
        seq_contrast_min    = min(sq$median), seq_contrast_max = max(sq$median),
        spearman_median     = sm$median, spearman_q025 = sm$q025, spearman_q975 = sm$q975
      )
    }))
  }))

  f <- file.path(path_csv, "ex3-summary.csv")
  fwrite(summ, f)
  msg("Saved: ", f)

  invisible(TRUE)
}
