# ============================================================
# 09_optimism-distribution.R
#
# PURPOSE
#   Describes the raw distribution of dispositional optimism, measured by the
#   score of the Life Orientation Test-Revised (LOT-R, 0 to 24). The score is
#   shown as a histogram of the raw data, with one bar per score, for three
#   samples: all participants of the FN treatment, the confirmatory sample of
#   the FN treatment, and all participants of the FP treatment. A summary
#   table gives the distribution of the score in each sample, and two tests
#   compare the FN and the FP treatment. Because the FP sample has a more
#   balanced sex ratio than the FN sample, the mean score is also given for
#   women and men separately.
#
# INPUT
#   path_der/lotr-scored.csv              LOT-R score per participant
#   path_src/participants.csv             treatment and sex
#   model "mpl": scored-<tr>.csv          Holt-Laury consistency flag, per treatment
#
# OUTPUT
#   data/output/figures/optimism-distribution.png
#   data/output/csv/optimism-distribution-summary.csv
#                               one row per sample: number of participants,
#                               minimum, quartiles, maximum, mean and standard
#                               deviation
#   data/output/csv/optimism-distribution-tests.csv
#                               one row per test of the FN and the FP treatment
#                               (all participants): Wilcoxon rank-sum test and
#                               Kolmogorov-Smirnov test of the LOT-R score
#   data/output/csv/optimism-distribution-by-sex.csv
#                               one row per sample and sex (women, men): number
#                               of participants and mean and standard deviation
#                               of the score; participants of other or
#                               unreported sex are not shown
#
# NOTES
#   - The confirmatory sample consists of the FN participants with a consistent
#     Holt-Laury choice pattern. The FP treatment has no confirmatory sample
#   - Cronbach's alpha of the scale by treatment is in
#     data/output/csv/processing_lotr-alpha.csv (scripts/processing/21_score-lotr.R)
#   - The row labels of the figure are written vertically (FN Full, FN
#     Confirmatory, FP Full); treatment colours: _reporting-style.R
# ============================================================

optimism_distribution <- function(cfg) {

  style <- treatment_style(cfg)

  tr_fn <- style$treat[style$type == "FN"]
  tr_fp <- style$treat[style$type == "FP"]
  stopifnot(length(tr_fn) == 1L, length(tr_fp) == 1L)

  # ---- Scores, with treatment, sex and the consistency flag ----
  lotr <- fread(file.path(path_der, "lotr-scored.csv"), encoding = "UTF-8")
  lotr[, pid := as.character(pid)]
  par <- fread(file.path(path_src, "participants.csv"), encoding = "UTF-8")
  par[, `:=`(pid = as.character(pid), treat = as.character(treat))]
  par[, sex := fifelse(sex %in% c("F", "M"), sex, "O")]
  d <- merge(lotr[, .(pid, score = as.numeric(lotr_score))],
             par[, .(pid, treat, sex)], by = "pid")
  d <- d[!is.na(score)]

  cons <- rbindlist(lapply(style$treat, function(tr) {
    sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
    sc[, .(pid = as.character(pid), consistent = as.integer(inconsistent) == 0L)]
  }))
  d <- merge(d, cons, by = "pid", all.x = TRUE)

  # ---- Samples: label, treatment, restriction to consistent, colour ----
  samples <- data.table(
    sample     = c("FN, full", "FN, confirmatory", "FP, full"),
    treat      = c(tr_fn, tr_fn, tr_fp),
    consistent = c(FALSE, TRUE, FALSE),
    fill       = c(style[treat == tr_fn, colour], style[treat == tr_fn, colour_light],
                   style[treat == tr_fp, colour])
  )

  dat <- rbindlist(lapply(seq_len(nrow(samples)), function(k) {
    x <- d[treat == samples$treat[k]]
    if (samples$consistent[k]) x <- x[consistent == TRUE]
    x[, sample := samples$sample[k]]
    x
  }))
  dat[, sample := factor(sample, levels = samples$sample)]

  # ---- Summary table ----
  summary_tbl <- dat[, .(n = .N, min = min(score), q25 = unname(quantile(score, 0.25)),
                         median = median(score), q75 = unname(quantile(score, 0.75)),
                         max = max(score), mean = mean(score), sd = sd(score)),
                     by = sample]
  setorder(summary_tbl, sample)
  f_csv <- file.path(path_csv, "optimism-distribution-summary.csv")
  fwrite(summary_tbl, f_csv)
  msg("Saved: ", f_csv)

  # ---- Tests: FN versus FP, all participants ----
  a <- dat[sample == "FN, full", score]; b <- dat[sample == "FP, full", score]
  w  <- suppressWarnings(wilcox.test(a, b, exact = FALSE))
  ks <- suppressWarnings(ks.test(a, b, exact = FALSE))
  tests <- data.table(
    comparison = "FN full vs FP full", test = c("Wilcoxon rank-sum", "Kolmogorov-Smirnov"),
    n_fn = length(a), n_fp = length(b), mean_fn = mean(a), mean_fp = mean(b),
    statistic = c(unname(w$statistic), unname(ks$statistic)),
    p_value   = c(w$p.value, ks$p.value))
  f_tests <- file.path(path_csv, "optimism-distribution-tests.csv")
  fwrite(tests, f_tests)
  msg("Saved: ", f_tests)

  # ---- Mean score by sex ----
  by_sex <- dat[sex %in% c("F", "M"),
                .(n = .N, mean = mean(score), sd = if (.N >= 2L) sd(score) else NA_real_),
                by = .(sample, sex)]
  setorder(by_sex, sample, sex)
  f_sex <- file.path(path_csv, "optimism-distribution-by-sex.csv")
  fwrite(by_sex, f_sex)
  msg("Saved: ", f_sex)

  # ---- Histogram: one bar per score, one row per sample ----
  fills <- setNames(samples$fill, samples$sample)
  h <- dat[, .N, by = .(sample, score)]
  h[, share := N / sum(N), by = sample]
  gm <- dat[, .(m = mean(score), n = .N), by = sample]
  gm[, lab := sprintf("mean %.1f", m)]
  gm[, hj  := fifelse(m > 12, 1.08, -0.08)]

  fig <- ggplot(h, aes(x = score, y = share, fill = sample)) +
    geom_col(width = 1, colour = "white", linewidth = 0.25, alpha = 0.9) +
    geom_vline(data = gm, aes(xintercept = m), colour = "grey20", linewidth = 0.7) +
    geom_text(data = gm, aes(x = m, y = Inf, label = lab, hjust = hj), inherit.aes = FALSE,
              size = 3, colour = "grey20", vjust = 1.6) +
    geom_text(data = gm, aes(x = -Inf, y = Inf, label = paste0("n = ", n)), inherit.aes = FALSE,
              size = 3, colour = "grey35", hjust = -0.15, vjust = 1.6) +
    scale_fill_manual(values = fills, guide = "none") +
    scale_x_continuous(name = "LOT-R score (0 to 24)", limits = c(-0.5, 24.5),
                       breaks = seq(0, 24, by = 4), expand = expansion(mult = 0.01)) +
    scale_y_continuous(name = "Share of participants", expand = expansion(mult = c(0, 0.25)),
                       labels = function(x) paste0(round(100 * x), "%")) +
    facet_grid(sample ~ ., switch = "y",
               labeller = as_labeller(c("FN, full" = "FN Full",
                                        "FN, confirmatory" = "FN Confirmatory",
                                        "FP, full" = "FP Full"))) +
    theme_classic(base_size = 11) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 90, hjust = 0.5, size = 10),
          axis.title = element_text(size = 10),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

  f_png <- file.path(path_fig, "optimism-distribution.png")
  ggsave(f_png, fig, width = 7, height = 7.5, dpi = 300)
  msg("Saved: ", f_png)

  invisible(summary_tbl)
}
