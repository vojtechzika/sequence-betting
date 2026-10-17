# ============================================================
# 06_behavior-stakes-by-sequence.R
#
# PURPOSE
#   Figure with two panels showing, for every sequence, the betting rate
#   (left) and the mean stake conditional on betting (right) in each
#   treatment. Points of the same sequence are connected across treatments;
#   vertical lines mark the means of each treatment over all trials.
#   Sequences are ordered by the values of the first treatment of
#   cfg$run$treatment.
#
# INPUT
#   path_der/master-sequences.csv
#
# OUTPUT
#   data/output/figures/behavior-stakes-by-sequence.png
#
# NOTES
#   - Treatment labels and colours: _reporting-style.R
#   - Each panel is sorted by its own variable (first treatment), so the
#     order of sequences differs between the panels; labels are shown on both
#   - The vertical lines are trial-level means (every trial has the same
#     weight), not means of participant means
#   - The stake axis runs from 0 to the endowment (cfg$design$seq$endowment)
# ============================================================

behavior_stakes_by_sequence <- function(cfg) {
  
  style <- treatment_style(cfg)
  stopifnot(nrow(style) == 2L)
  e <- as.numeric(cfg$design$seq$endowment)
  
  master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
  master[, seq   := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, treat := as.character(treat)]
  master[, treatment := factor(style$label[match(treat, style$treat)], levels = style$label)]
  
  tr_colours <- setNames(style$colour, style$label)
  lab_a <- style$label[1]; lab_b <- style$label[2]
  
  # ---- Sequence-level aggregates ----
  seq_tbl <- master[, .(
    n_trials       = .N,
    bet_rate       = mean(stake > 0),
    mean_stake_pos = if (any(stake > 0)) mean(stake[stake > 0]) else NA_real_
  ), by = .(seq, treatment)]
  
  first <- seq_tbl[treatment == lab_a]
  seq_order_bet   <- first[order(bet_rate), seq]
  seq_order_stake <- first[order(mean_stake_pos, na.last = FALSE), seq]
  
  wide <- function(var) {
    w <- dcast(seq_tbl[!is.na(get(var))], seq ~ treatment, value.var = var)
    setnames(w, c(lab_a, lab_b), c("a", "b"))
    w
  }
  wide_bet <- wide("bet_rate")
  wide_stk <- wide("mean_stake_pos")
  
  gm <- master[, .(gm_bet_rate = mean(stake > 0),
                   gm_mean_stake_pos = mean(stake[stake > 0])), by = treatment]
  
  panel <- function(wide_tbl, points, xvar, gmvar, x_scale, seq_order, right_margin) {
    ggplot() +
      geom_vline(data = gm, aes(xintercept = .data[[gmvar]], colour = treatment), linewidth = 0.6) +
      geom_segment(data = wide_tbl, aes(x = a, xend = b, y = seq, yend = seq),
                   colour = "grey75", linewidth = 0.3) +
      geom_point(data = points, aes(x = .data[[xvar]], y = seq, colour = treatment),
                 size = 1.5, alpha = 0.85) +
      scale_colour_manual(values = tr_colours, name = NULL) +
      x_scale +
      scale_y_discrete(limits = seq_order) +
      theme_classic(base_size = 10) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(size = 6.5, family = "mono"),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
            plot.margin = right_margin)
  }
  
  p_bet <- panel(wide_bet, seq_tbl, "bet_rate", "gm_bet_rate",
                 scale_x_continuous(name = "Betting rate", labels = percent_format(accuracy = 1),
                                    limits = c(0, 1)),
                 seq_order_bet, margin(4, 8, 4, 8))
  p_stk <- panel(wide_stk, seq_tbl[!is.na(mean_stake_pos)], "mean_stake_pos", "gm_mean_stake_pos",
                 scale_x_continuous(name = "Mean stake | bet (ECU)", limits = c(0, e),
                                    breaks = seq(0, e, by = e / 5)),
                 seq_order_stake, margin(4, 8, 4, 4))
  
  final <- p_bet + p_stk + plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.text = element_text(size = 9))
  
  f_out <- file.path(path_fig, "behavior-stakes-by-sequence.png")
  ggsave(f_out, final, width = 12, height = 14, dpi = 300)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}
