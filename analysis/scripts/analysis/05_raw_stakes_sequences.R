# ============================================================
# 76_sequence_descriptives.R
# ============================================================

raw_stakes_sequences <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(gridExtra)
  
  master <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  master[, seq := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, treat := as.character(treat)]
  master[, treatment := factor(
    fifelse(treat == "m25", "FN (m = 2.5)", "FP (m = 1.9)"),
    levels = c("FN (m = 2.5)", "FP (m = 1.9)")
  )]
  
  tr_colours <- c("FN (m = 2.5)" = "#2166AC", "FP (m = 1.9)" = "#C0392B")
  
  # ---- Sequence-level aggregates -----------------------------------------
  seq_tbl <- master[, .(
    n_trials       = .N,
    bet_rate       = mean(stake > 0),
    mean_stake_pos = ifelse(any(stake > 0), mean(stake[stake > 0]), NA_real_)
  ), by = .(seq, treatment)]
  
  fn_data <- seq_tbl[treatment == "FN (m = 2.5)"]
  
  # Ascending order in levels vector = highest value at TOP of ggplot
  # (ggplot draws factor level 1 at bottom, last level at top)
  seq_limits_bet   <- fn_data[order(bet_rate),         seq]
  seq_limits_stake <- fn_data[order(mean_stake_pos, na.last = FALSE), seq]
  
  # Wide tables for grey connectors — plain character, ordering via scale
  seq_wide_bet <- dcast(seq_tbl, seq ~ treatment, value.var = "bet_rate")
  setnames(seq_wide_bet, c("FN (m = 2.5)", "FP (m = 1.9)"), c("fn", "fp"))
  
  seq_wide_stk <- dcast(seq_tbl[!is.na(mean_stake_pos)],
                        seq ~ treatment, value.var = "mean_stake_pos")
  setnames(seq_wide_stk, c("FN (m = 2.5)", "FP (m = 1.9)"), c("fn", "fp"))
  
  # Grand means
  gm <- master[, .(
    gm_bet_rate       = mean(stake > 0),
    gm_mean_stake_pos = mean(stake[stake > 0])
  ), by = treatment]
  
  # ---- Panel 1: Betting rate ----------------------------------------------
  p_bet <- ggplot() +
    
    geom_vline(data = gm,
               aes(xintercept = gm_bet_rate, colour = treatment),
               linewidth = 0.6) +
    
    geom_segment(data = seq_wide_bet,
                 aes(x = fn, xend = fp, y = seq, yend = seq),
                 colour = "grey75", linewidth = 0.3) +
    
    geom_point(data = seq_tbl,
               aes(x = bet_rate, y = seq, colour = treatment),
               size = 1.5, alpha = 0.85) +
    
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(name   = "Betting rate",
                       labels = percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_y_discrete(limits = seq_limits_bet) +
    
    theme_classic(base_size = 10) +
    theme(
      axis.title.y       = element_blank(),
      axis.text.y        = element_text(size = 6.5, family = "mono"),
      axis.ticks.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      legend.position    = "none",
      plot.margin        = margin(4, 8, 4, 8)
    )
  
  # ---- Panel 2: Mean stake conditional on betting -------------------------
  p_stake <- ggplot() +
    
    geom_vline(data = gm,
               aes(xintercept = gm_mean_stake_pos, colour = treatment),
               linewidth = 0.6) +
    
    geom_segment(data = seq_wide_stk,
                 aes(x = fn, xend = fp, y = seq, yend = seq),
                 colour = "grey75", linewidth = 0.3) +
    
    geom_point(data = seq_tbl[!is.na(mean_stake_pos)],
               aes(x = mean_stake_pos, y = seq, colour = treatment),
               size = 1.5, alpha = 0.85) +
    
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(name   = "Mean stake | bet (ECU)",
                       limits = c(0, 100),
                       breaks = seq(0, 100, by = 20)) +
    scale_y_discrete(limits = seq_limits_stake) +
    
    theme_classic(base_size = 10) +
    theme(
      axis.title.y       = element_blank(),
      axis.text.y        = element_text(size = 6.5, family = "mono"),
      axis.ticks.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      legend.position    = "none",
      plot.margin        = margin(4, 8, 4, 4)
    )
  
  # ---- Shared legend ------------------------------------------------------
  legend_plot <- ggplot(
    data.table(x = 1:2, y = 0,
               treatment = factor(c("FN (m = 2.5)", "FP (m = 1.9)"),
                                  levels = c("FN (m = 2.5)", "FP (m = 1.9)"))),
    aes(x = x, y = y, colour = treatment)
  ) +
    geom_point(size = 3) +
    scale_colour_manual(values = tr_colours, name = NULL) +
    coord_cartesian(xlim = c(10, 20)) +
    theme_void() +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.text      = element_text(size = 9))
  
  legend_grob <- ggplotGrob(legend_plot)
  leg_idx     <- which(legend_grob$layout$name %in%
                         c("guide-box", "guide-box-bottom"))
  leg         <- legend_grob$grobs[[leg_idx[1]]]
  
  # ---- Assemble -----------------------------------------------------------
  panels <- p_bet + p_stake + plot_layout(ncol = 2, widths = c(1, 1))
  
  final <- arrangeGrob(
    patchworkGrob(panels),
    leg,
    ncol    = 1,
    heights = unit(c(1, 0.04), c("null", "npc"))
  )
  
  f_out <- file.path(path_fig, "sequence_descriptives_by_treatment.png")
  ggsave(f_out, plot = final, width = 12, height = 14, dpi = 300)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}