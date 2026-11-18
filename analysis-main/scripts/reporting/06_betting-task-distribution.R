# ============================================================
# 06_betting-task-distribution.R
#
# PURPOSE
#   Describes the raw distributions of three quantities of the sequence
#   betting task: the betting rate of a participant (share of the 64
#   sequences in which the participant bets), the side choice of a participant
#   (share of Heads among the bets placed), and the stake of a bet.
#   The quantities are shown as histograms of the raw data, for three
#   samples: all participants of the FN treatment, the confirmatory sample of
#   the FN treatment, and all participants of the FP treatment.
#   One figure arranges the quantities in columns and the samples in rows, so
#   that the samples can be compared within a quantity. A summary table gives
#   the distribution of each quantity in each sample.
#
# INPUT
#   path_der/master-sequences.csv         trials of the sequence task
#   model "mpl": scored-<tr>.csv          Holt-Laury consistency flag, per treatment
#
# OUTPUT
#   data/output/figures/betting-task-distribution.png
#   data/output/csv/betting-task-distribution-summary.csv
#                               one row per quantity and sample: number of
#                               observations, minimum, quartiles, maximum, mean
#                               and standard deviation, and the shares of
#                               observations at the boundaries of the quantity
#                               (share_at_0 and share_at_1 for the betting
#                               rate, share_at_half for the side choice,
#                               share_at_50 and share_at_100 for the stake)
#
# NOTES
#   - A bet is a trial with a positive stake; a missing stake counts as no bet
#   - Betting rate and side choice are participant-level, so the histograms
#     show the share of participants and the standard deviations are across
#     participants; the stake is trial-level, so its histogram shows the share
#     of bets and its standard deviation is across bets
#   - Stakes are shown in bins of 5 ECU (right-closed), so that stakes of 50
#     and 100 ECU fall in the bins ending at these values
#   - The side choice is shown for participants with at least
#     min_bets_side bets; the number of participants dropped is reported in the
#     summary (n_dropped) because shares based on few bets are mostly noise
#   - The confirmatory sample consists of the FN participants with a consistent
#     Holt-Laury choice pattern; all FN participants are normative betters
#     (see scripts/modelling/07_a-star-draws.R), so no further restriction applies. The FP
#     treatment has no confirmatory sample
#   - The row labels of the figure are written vertically (FN Full, FN
#     Confirmatory, FP Full) to leave more room for the panels
#   - Treatment labels and colours: _reporting-style.R; the confirmatory sample
#     uses the lighter shade of the FN colour
#   - Standard deviations of the summary differ from those of the descriptive
#     table of the main text, which are across trials for the sequence rows
# ============================================================

betting_task_distribution <- function(cfg) {
  
  style <- treatment_style(cfg)
  e     <- as.numeric(cfg$design$seq$endowment)
  heads_label   <- as.character(cfg$design$seq$side_labels$heads)
  min_bets_side <- 10L
  
  tr_fn <- style$treat[style$type == "FN"]
  tr_fp <- style$treat[style$type == "FP"]
  stopifnot(length(tr_fn) == 1L, length(tr_fp) == 1L)
  
  # ---- Trials, with the consistency flag of the participant ----
  master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
  master[, `:=`(pid = as.character(pid), treat = as.character(treat),
                stake = as.numeric(stake))]
  master[is.na(stake), stake := 0]
  master[, bet := stake > 0]
  
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
  
  # ---- Participant-level betting rate and side choice ----
  part <- trials[, .(
    n_trials = .N,
    n_bets   = sum(bet),
    rate     = mean(bet),
    heads    = if (sum(bet) > 0) mean(side[bet] == heads_label) else NA_real_
  ), by = .(sample, pid)]
  
  rate_d  <- part[, .(sample, value = rate)]
  side_d  <- part[n_bets >= min_bets_side, .(sample, value = heads)]
  stake_d <- trials[bet == TRUE, .(sample, value = stake)]
  
  # ---- Summary table ----
  describe <- function(x) {
    data.table(n = length(x), min = min(x), q25 = unname(quantile(x, 0.25)),
               median = median(x), q75 = unname(quantile(x, 0.75)), max = max(x),
               mean = mean(x), sd = if (length(x) >= 2L) sd(x) else NA_real_)
  }
  summ <- function(d, quantity, extra) {
    rbindlist(lapply(levels(d$sample), function(s) {
      x <- d[sample == s, value]
      cbind(data.table(quantity = quantity, sample = s), describe(x), extra(x))
    }), fill = TRUE)
  }
  
  s_rate  <- summ(rate_d, "betting rate", function(x)
    data.table(share_at_0 = mean(x == 0), share_at_1 = mean(x == 1)))
  s_stake <- summ(stake_d, "stake", function(x)
    data.table(share_at_50 = mean(x == 50), share_at_100 = mean(x == e)))
  s_side  <- summ(side_d, "share of Heads", function(x)
    data.table(share_at_half = mean(abs(x - 0.5) < 1e-9)))
  summary_tbl <- rbindlist(list(s_rate, s_stake, s_side), fill = TRUE)
  
  # participants dropped from the side choice for having too few bets
  n_part <- part[, .(n_participants = .N, n_dropped = sum(n_bets < min_bets_side)), by = sample]
  n_part[, sample := as.character(sample)]
  summary_tbl <- merge(summary_tbl, n_part, by = "sample", all.x = TRUE)
  summary_tbl[quantity != "share of Heads", n_dropped := NA_integer_]
  summary_tbl[, sample := factor(sample, levels = samples$sample)]
  summary_tbl[, quantity := factor(quantity, levels = c("betting rate", "share of Heads", "stake"))]
  setorder(summary_tbl, quantity, sample)
  
  f_csv <- file.path(path_csv, "betting-task-distribution-summary.csv")
  fwrite(summary_tbl, f_csv)
  msg("Saved: ", f_csv)
  
  # ---- Histogram of one quantity: one row per sample ----
  fills <- setNames(samples$fill, samples$sample)
  
  hist_col <- function(d, bw, x_lim, x_breaks, x_labels, x_name, y_name,
                       centre, mean_label, ref = NULL, first = FALSE) {
    d <- copy(d)
    d[, bin := centre(value)]
    h <- d[, .N, by = .(sample, bin)]
    h[, share := N / sum(N), by = sample]
    
    gm <- d[, .(m = mean(value), n = .N), by = sample]
    mid <- mean(x_lim)
    gm[, lab := mean_label(m)]
    gm[, hj := fifelse(m > mid, 1.08, -0.08)]
    
    g <- ggplot(h, aes(x = bin, y = share, fill = sample)) +
      geom_col(width = bw, colour = "white", linewidth = 0.25, alpha = 0.9) +
      geom_vline(data = gm, aes(xintercept = m), colour = "grey20", linewidth = 0.7) +
      geom_text(data = gm, aes(x = m, y = Inf, label = lab, hjust = hj), inherit.aes = FALSE,
                size = 3, colour = "grey20", vjust = 1.6) +
      geom_text(data = gm, aes(x = -Inf, y = Inf, label = paste0("n = ", n)), inherit.aes = FALSE,
                size = 3, colour = "grey35", hjust = -0.15, vjust = 1.6) +
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
    if (!is.null(ref)) {
      g <- g + geom_vline(xintercept = ref, linetype = "dashed", colour = "grey60", linewidth = 0.4)
    }
    g
  }
  
  bw_p <- 0.05                                   # bins of the shares
  centre_share <- function(x) pmin(floor(x / bw_p + 1e-9), 1 / bw_p - 1) * bw_p + bw_p / 2
  bw_s <- 5                                      # bins of the stakes (ECU)
  centre_stake <- function(x) ceiling(x / bw_s) * bw_s - bw_s / 2
  pct <- function(x) sprintf("%.0f%%", 100 * x)
  
  p_rate <- hist_col(rate_d, bw_p, c(0, 1), seq(0, 1, 0.25), scales::percent_format(accuracy = 1),
                     "Betting rate (share of sequences with a bet)", "Share of participants",
                     centre_share, function(m) paste0("mean ", pct(m)), first = TRUE)
  p_stake <- hist_col(stake_d, bw_s, c(0, e), seq(0, e, by = e / 5), waiver(),
                      "Stake of a bet (ECU)", "Share of bets",
                      centre_stake, function(m) sprintf("mean %.1f", m))
  p_side <- hist_col(side_d, bw_p, c(0, 1), seq(0, 1, 0.25), scales::percent_format(accuracy = 1),
                     "Share of Heads among bets", "Share of participants",
                     centre_share, function(m) paste0("mean ", pct(m)), ref = 0.5)
  
  fig <- p_rate + p_side + p_stake + plot_layout(ncol = 3)
  f_png <- file.path(path_fig, "betting-task-distribution.png")
  ggsave(f_png, fig, width = 13, height = 7.5, dpi = 300)
  msg("Saved: ", f_png)
  
  invisible(summary_tbl)
}
