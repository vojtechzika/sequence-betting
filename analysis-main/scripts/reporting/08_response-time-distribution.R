# ============================================================
# 08_response-time-distribution.R
#
# PURPOSE
#   Describes the raw distribution of the response time in the sequence
#   betting task (seconds from the onset of the decision screen to the
#   confirmation of the choice). The response time is shown at two levels:
#   the response time of a trial and the geometric mean response time of a
#   participant, that is, the exponential of the mean log response time over
#   all trials of the participant, which is the response-time predictor of the
#   exploratory analysis EX2.
#   Both are shown as histograms of the raw data for three samples: all
#   participants of the FN treatment, the confirmatory sample of the FN
#   treatment, and all participants of the FP treatment. One figure arranges
#   the two levels in columns and the samples in rows. A summary table gives
#   the distribution of the response time of a trial for all trials, for
#   trials with a bet and for trials without a bet, and of the participant
#   geometric means. The treatments are compared with rank-sum tests on participants.
#
# INPUT
#   path_der/master-sequences.csv         trials of the sequence task
#   model "mpl": scored-<tr>.csv          Holt-Laury consistency flag, per treatment
#
# OUTPUT
#   data/output/figures/response-time-distribution.png
#   data/output/csv/response-time-distribution-summary.csv
#                               one row per level and sample: number of
#                               observations, minimum, quartiles, maximum, mean
#                               and standard deviation, and the share of
#                               observations below fast_s seconds
#   data/output/csv/response-time-distribution-tests.csv
#                               one row per test of the FN and the FP treatment
#                               (all participants): Wilcoxon rank-sum test of
#                               the participant geometric mean response time and
#                               of the participant share of trials below fast_s seconds
#
# NOTES
#   - A bet is a trial with a positive stake; a missing stake counts as no bet
#   - The trial-level histogram uses a logarithmic axis because response times
#     are strongly right-skewed; the bars show the share of trials in bins of
#     0.1 on the base-10 logarithm, and the mean marks the arithmetic mean
#   - The participant-level histogram shows the share of participants in bins
#     of 1 second; the participant geometric mean is computed over all trials of
#     the participant, with and without a bet, as in EX2
#   - The standard deviation is across trials for the trial-level rows and
#     across participants for the participant geometric means; the mean of the
#     participant level is the arithmetic mean of the geometric means
#   - Trials without a positive recorded response time are excluded
#   - The confirmatory sample consists of the FN participants with a consistent
#     Holt-Laury choice pattern. The FP treatment has no confirmatory sample
#   - The row labels of the figure are written vertically (FN Full, FN
#     Confirmatory, FP Full); treatment colours: _reporting-style.R
# ============================================================

response_time_distribution <- function(cfg) {

  style  <- treatment_style(cfg)
  fast_s <- 2                                    # threshold of a fast response (s)

  tr_fn <- style$treat[style$type == "FN"]
  tr_fp <- style$treat[style$type == "FP"]
  stopifnot(length(tr_fn) == 1L, length(tr_fp) == 1L)

  # ---- Trials, with the consistency flag of the participant ----
  master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
  master[, `:=`(pid = as.character(pid), treat = as.character(treat),
                stake = as.numeric(stake), rt = as.numeric(screen_ms) / 1000)]
  master[is.na(stake), stake := 0]
  master[, bet := stake > 0]
  master <- master[is.finite(rt) & rt > 0]

  cons <- rbindlist(lapply(style$treat, function(tr) {
    sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
    sc[, .(pid = as.character(pid), consistent = as.integer(inconsistent) == 0L)]
  }))
  master <- merge(master, cons, by = "pid", all.x = TRUE)

  # ---- Samples: label, treatment, restriction to consistent, colour ----
  samples <- data.table(
    sample     = c("FN, full", "FN, confirmatory", "FP, full"),
    treat      = c(tr_fn, tr_fn, tr_fp),
    consistent = c(FALSE, TRUE, FALSE),
    fill       = c(style[treat == tr_fn, colour], style[treat == tr_fn, colour_light],
                   style[treat == tr_fp, colour])
  )

  trials <- rbindlist(lapply(seq_len(nrow(samples)), function(k) {
    d <- master[treat == samples$treat[k]]
    if (samples$consistent[k]) d <- d[consistent == TRUE]
    d[, sample := samples$sample[k]]
    d
  }))
  trials[, sample := factor(sample, levels = samples$sample)]

  part <- trials[, .(geo_rt = exp(mean(log(rt))), share_fast = mean(rt < fast_s)),
                 by = .(sample, pid)]

  # ---- Summary table ----
  describe <- function(x) {
    data.table(n = length(x), min = min(x), q25 = unname(quantile(x, 0.25)),
               median = median(x), q75 = unname(quantile(x, 0.75)), max = max(x),
               mean = mean(x), sd = if (length(x) >= 2L) sd(x) else NA_real_,
               share_below_fast = mean(x < fast_s))
  }
  summ <- function(d, level) {
    rbindlist(lapply(levels(d$sample), function(s) {
      cbind(data.table(level = level, sample = s), describe(d[sample == s, value]))
    }))
  }
  summary_tbl <- rbindlist(list(
    summ(trials[, .(sample, value = rt)],                 "trial, all"),
    summ(trials[bet == TRUE,  .(sample, value = rt)],     "trial, bet"),
    summ(trials[bet == FALSE, .(sample, value = rt)],     "trial, no bet"),
    summ(part[, .(sample, value = geo_rt)],              "participant geometric mean")
  ))
  summary_tbl[level == "participant geometric mean", share_below_fast := NA_real_]
  summary_tbl[, sample := factor(sample, levels = samples$sample)]
  summary_tbl[, level := factor(level, levels = c("trial, all", "trial, bet",
                                                  "trial, no bet", "participant geometric mean"))]
  setorder(summary_tbl, level, sample)

  f_csv <- file.path(path_csv, "response-time-distribution-summary.csv")
  fwrite(summary_tbl, f_csv)
  msg("Saved: ", f_csv)

  # ---- Tests: FN versus FP, all participants ----
  fn <- part[sample == "FN, full"]; fp <- part[sample == "FP, full"]
  wt <- function(a, b, what) {
    w <- suppressWarnings(wilcox.test(a, b, exact = FALSE))
    data.table(comparison = "FN full vs FP full", variable = what,
               n_fn = length(a), n_fp = length(b),
               mean_fn = mean(a), mean_fp = mean(b),
               statistic = unname(w$statistic), p_value = w$p.value)
  }
  tests <- rbindlist(list(
    wt(fn$geo_rt,    fp$geo_rt,    "participant geometric mean response time (s)"),
    wt(fn$share_fast, fp$share_fast, paste0("participant share of trials below ", fast_s, " s"))
  ))
  f_tests <- file.path(path_csv, "response-time-distribution-tests.csv")
  fwrite(tests, f_tests)
  msg("Saved: ", f_tests)

  # ---- Histogram of one level: one row per sample ----
  fills <- setNames(samples$fill, samples$sample)

  hist_col <- function(d, tf, bw, x_lim, x_breaks, x_labels, x_name, y_name,
                       mean_label, first = FALSE) {
    d <- copy(d)
    d[, x   := tf(value)]
    d[, bin := floor(x / bw + 1e-9) * bw + bw / 2]
    h <- d[, .N, by = .(sample, bin)]
    h[, share := N / sum(N), by = sample]

    gm <- d[, .(m = mean(value), n = .N), by = sample]
    gm[, pos := tf(m)]
    gm[, lab := mean_label(m)]
    gm[, hj  := fifelse(pos > mean(x_lim), 1.08, -0.08)]

    ggplot(h, aes(x = bin, y = share, fill = sample)) +
      geom_col(width = bw, colour = "white", linewidth = 0.25, alpha = 0.9) +
      geom_vline(data = gm, aes(xintercept = pos), colour = "grey20", linewidth = 0.7) +
      geom_text(data = gm, aes(x = pos, y = Inf, label = lab, hjust = hj), inherit.aes = FALSE,
                size = 3, colour = "grey20", vjust = 1.6) +
      geom_text(data = gm, aes(x = Inf, y = Inf, label = paste0("n = ", n)), inherit.aes = FALSE,
                size = 3, colour = "grey35", hjust = 1.15, vjust = 1.6) +
      scale_fill_manual(values = fills, guide = "none") +
      scale_x_continuous(name = x_name, limits = x_lim, breaks = x_breaks, labels = x_labels,
                         expand = expansion(mult = 0.01)) +
      scale_y_continuous(name = y_name, expand = expansion(mult = c(0, 0.25)),
                         labels = function(x) paste0(round(100 * x), "%")) +
      facet_grid(sample ~ ., switch = "y",
                 labeller = as_labeller(c("FN, full" = "FN Full",
                                          "FN, confirmatory" = "FN Confirmatory",
                                          "FP, full" = "FP Full"))) +
      theme_classic(base_size = 11) +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = if (first) element_text(angle = 90, hjust = 0.5, size = 10)
                                else element_blank(),
            axis.title = element_text(size = 10),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
  }

  # trial level on a base-10 logarithmic axis
  bw_t   <- 0.1
  lim_t  <- c(floor(min(log10(trials$rt)) / bw_t) * bw_t,
              ceiling(max(log10(trials$rt)) / bw_t) * bw_t)
  secs   <- c(0.5, 1, 2, 5, 10, 20, 50, 100, 200)
  secs   <- secs[log10(secs) >= lim_t[1] & log10(secs) <= lim_t[2]]
  p_trial <- hist_col(trials[, .(sample, value = rt)], log10, bw_t, lim_t, log10(secs), secs,
                      "Response time of a trial (seconds, logarithmic axis)", "Share of trials",
                      function(m) sprintf("mean %.1f s", m), first = TRUE)

  # participant level on a linear axis
  bw_p   <- 1
  lim_p  <- c(0, ceiling(max(part$geo_rt) / 5) * 5)
  p_part <- hist_col(part[, .(sample, value = geo_rt)], identity, bw_p, lim_p,
                     seq(lim_p[1], lim_p[2], by = 5), waiver(),
                     "Geometric mean response time of a participant (seconds)", "Share of participants",
                     function(m) sprintf("mean %.1f s", m))

  fig <- p_trial + p_part + plot_layout(ncol = 2)
  f_png <- file.path(path_fig, "response-time-distribution.png")
  ggsave(f_png, fig, width = 10, height = 7.5, dpi = 300)
  msg("Saved: ", f_png)

  invisible(summary_tbl)
}
