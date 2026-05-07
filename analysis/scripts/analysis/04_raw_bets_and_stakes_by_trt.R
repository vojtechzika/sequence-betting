# ============================================================
# 04_raw_stakes_participants.R
#
# PURPOSE
#   Combined figure with two panels:
#   Left:  Distribution of participant-level average betting rates
#          by treatment (all decisions, including no-bets), with KS test.
#   Right: Overlapping distribution of raw stake amounts (1-100 ECU)
#          by treatment (positive stakes only, no-bets excluded), with KS test.
#
# INPUT
#   path_src/master_sequences.csv
#
# OUTPUT
#   path_fig/desc_bets_and_stakes_by_trt.png
# ============================================================
raw_bets_and_stakes_by_trt <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  
  master <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  master[, stake     := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, treat     := as.character(treat)]
  master[, treatment := factor(
    fifelse(treat == "m25", "FN (m = 2.5)", "FP (m = 1.9)"),
    levels = c("FP (m = 1.9)", "FN (m = 2.5)")
  )]
  master[, bet := as.integer(stake > 0)]
  
  tr_colours <- c("FN (m = 2.5)" = "#2166AC", "FP (m = 1.9)" = "#C0392B")
  
  # ---- Left panel: participant-level average betting rates ----------------
  pid_tbl <- master[, .(
    mean_bet  = mean(bet),
    treatment = treatment[1]
  ), by = pid]
  
  sum_bet <- pid_tbl[, .(mean_rate = mean(mean_bet)), by = treatment]
  sum_bet[, label := paste0(treatment, ": ", scales::percent(mean_rate, accuracy = 0.1))]
  
  ks_bet <- ks.test(
    pid_tbl[treatment == "FN (m = 2.5)", mean_bet],
    pid_tbl[treatment == "FP (m = 1.9)", mean_bet]
  )
  ks_bet_label <- paste0("KS: D = ", round(ks_bet$statistic, 2),
                         ", p = ", round(ks_bet$p.value, 2))
  
  p_bet <- ggplot(pid_tbl, aes(x = mean_bet, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.4, adjust = 0.8, linewidth = 0, trim = TRUE) +
    geom_vline(data = sum_bet,
               aes(xintercept = mean_rate, colour = treatment),
               linewidth = 0.7, linetype = "solid") +
    geom_text(data = sum_bet[treatment == "FN (m = 2.5)"],
              aes(x = mean_rate, y = Inf, label = label, colour = treatment),
              hjust = -0.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    geom_text(data = sum_bet[treatment == "FP (m = 1.9)"],
              aes(x = mean_rate, y = Inf, label = label, colour = treatment),
              hjust = 1.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    scale_fill_manual(values = tr_colours, name = NULL) +
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(
      name   = paste0("Average Betting Rate\n(", ks_bet_label, ")"),
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_y_continuous(name = "Density") +
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position    = "none",
      plot.margin        = margin(4, 12, 8, 8)
    )
  
  # ---- Right panel: raw stake distribution (positive stakes only) ---------
  pos <- master[stake > 0]
  
  sum_tbl <- pos[, .(mean_stake = mean(stake)), by = treatment]
  sum_tbl[, label := paste0(treatment, ": ", round(mean_stake, 1), " ECU")]
  
  ks_stake <- ks.test(
    pos[treatment == "FN (m = 2.5)", stake],
    pos[treatment == "FP (m = 1.9)", stake]
  )

  ks_stake_label <- paste0("KS: D = ", round(ks_stake$statistic, 2),
                           ", p = ", round(ks_stake$p.value, 2))
  

  
  p_pos <- ggplot(pos, aes(x = stake, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.4, adjust = 0.8, linewidth = 0, trim = TRUE) +
    geom_vline(data = sum_tbl,
               aes(xintercept = mean_stake, colour = treatment),
               linewidth = 0.7, linetype = "solid") +
    geom_text(data = sum_tbl[treatment == "FN (m = 2.5)"],
              aes(x = mean_stake, y = Inf, label = label, colour = treatment),
              hjust = -0.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    geom_text(data = sum_tbl[treatment == "FP (m = 1.9)"],
              aes(x = mean_stake, y = Inf, label = label, colour = treatment),
              hjust = 1.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    scale_fill_manual(values = tr_colours, name = NULL) +
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(
      name   = paste0("Stake (ECU, conditional on betting)\n(", ks_stake_label, ")"),
      breaks = seq(0, 100, by = 20),
      limits = c(1, 100)
    ) +
    scale_y_continuous(name = "Density") +
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position    = "none",
      plot.margin        = margin(4, 12, 8, 8)
    )
  
  # ---- Combine ------------------------------------------------------------
  combined <- p_bet + p_pos +
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
      theme = theme(plot.margin = margin(4, 4, 4, 4))
    )
  
  f_out <- file.path(path_fig, "desc_bets_and_stakes_by_trt.png")
  ggsave(f_out, combined, width = 12, height = 4, dpi = 300)
  msg("Saved: ", f_out)
  
 
  
  invisible(TRUE)
}