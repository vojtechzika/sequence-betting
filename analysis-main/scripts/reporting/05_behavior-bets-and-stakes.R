# ============================================================
# 05_behavior-bets-and-stakes.R
#
# PURPOSE
#   Figure with two panels comparing the treatments.
#   Left:  distribution of participant-level betting rates (all decisions,
#          including no-bets).
#   Right: distribution of raw stakes (betting trials only).
#   Each panel marks the treatment means and reports a two-sample
#   Kolmogorov-Smirnov test. The tests use participant-level values (the
#   betting rate and the mean stake of each participant who bet at least
#   once), so that repeated trials of a participant are not treated as
#   independent observations. The stake density itself shows all betting
#   trials.
#
# INPUT
#   path_der/master-sequences.csv
#
# OUTPUT
#   data/output/figures/behavior-bets-and-stakes.png
#
# NOTES
#   - Treatment labels and colours: _reporting-style.R
#   - The test compares the first two treatments of cfg$run$treatment
#   - The stake axis runs from 1 to the endowment (cfg$design$seq$endowment)
# ============================================================

behavior_bets_and_stakes <- function(cfg) {
  
  style <- treatment_style(cfg)
  stopifnot(nrow(style) == 2L)
  e <- as.numeric(cfg$design$seq$endowment)
  
  master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, treat := as.character(treat)]
  master[, treatment := factor(style$label[match(treat, style$treat)],
                               levels = rev(style$label))]
  master[, bet := as.integer(stake > 0)]
  
  tr_colours <- setNames(style$colour, style$label)
  lab_a <- style$label[1]; lab_b <- style$label[2]
  
  ks_label <- function(x_a, x_b) {
    ks <- suppressWarnings(ks.test(x_a, x_b, exact = FALSE))
    paste0("KS: D = ", round(ks$statistic, 2), ", p = ", round(ks$p.value, 2))
  }
  
  density_panel <- function(d, var, summary_tbl, x_scale) {
    ggplot(d, aes(x = .data[[var]], fill = treatment, colour = treatment)) +
      geom_density(alpha = 0.4, adjust = 0.8, linewidth = 0, trim = TRUE) +
      geom_vline(data = summary_tbl, aes(xintercept = mean_val, colour = treatment),
                 linewidth = 0.7) +
      geom_text(data = summary_tbl[treatment == lab_a],
                aes(x = mean_val, y = Inf, label = label, colour = treatment),
                hjust = -0.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
      geom_text(data = summary_tbl[treatment == lab_b],
                aes(x = mean_val, y = Inf, label = label, colour = treatment),
                hjust = 1.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
      scale_fill_manual(values = tr_colours, name = NULL) +
      scale_colour_manual(values = tr_colours, name = NULL) +
      x_scale +
      scale_y_continuous(name = "Density") +
      theme_classic(base_size = 10) +
      theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "none", plot.margin = margin(4, 12, 8, 8))
  }
  
  # ---- Left panel: participant betting rates ----
  pid_bet <- master[, .(mean_bet = mean(bet), treatment = treatment[1], treat = treat[1]), by = pid]
  sum_bet <- pid_bet[, .(mean_val = mean(mean_bet)), by = treatment]
  sum_bet[, label := paste0(treatment, ": ", scales::percent(mean_val, accuracy = 0.1))]
  ks_bet <- ks_label(pid_bet[treat == style$treat[1], mean_bet],
                     pid_bet[treat == style$treat[2], mean_bet])
  
  p_bet <- density_panel(
    pid_bet, "mean_bet", sum_bet,
    scale_x_continuous(name = paste0("Average Betting Rate\n(", ks_bet, ", participant rates)"),
                       labels = scales::percent_format(accuracy = 1), limits = c(0, 1)))
  
  # ---- Right panel: stakes of betting trials ----
  pos <- master[stake > 0]
  sum_stk <- pos[, .(mean_val = mean(stake)), by = treatment]
  sum_stk[, label := paste0(treatment, ": ", round(mean_val, 1), " ECU")]
  pid_stk <- pos[, .(mean_stake = mean(stake), treat = treat[1]), by = pid]
  ks_stk <- ks_label(pid_stk[treat == style$treat[1], mean_stake],
                     pid_stk[treat == style$treat[2], mean_stake])
  
  p_stk <- density_panel(
    pos, "stake", sum_stk,
    scale_x_continuous(name = paste0("Stake (ECU, conditional on betting)\n(", ks_stk,
                                     ", participant means)"),
                       breaks = seq(0, e, by = e / 5), limits = c(1, e)))
  
  combined <- p_bet + p_stk + plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
  
  f_out <- file.path(path_fig, "behavior-bets-and-stakes.png")
  ggsave(f_out, combined, width = 12, height = 4, dpi = 300)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}
