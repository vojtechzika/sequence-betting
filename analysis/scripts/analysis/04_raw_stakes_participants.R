# ============================================================
# 75_stake_distribution.R
#
# PURPOSE
#   Overlapping distribution of raw stake amounts (0-100 ECU)
#   by treatment (FN vs FP), observation level.
#   Zero stakes (no-bet) shown separately as a bar to avoid
#   compressing the positive stake distribution.
#
# INPUT
#   path_src/master_sequences.csv
#
# OUTPUT
#   path_fig/stake_distribution_by_treatment.png
# ============================================================

raw_stakes_participants <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(scales)
  
  master <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  master[, stake  := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, treat  := as.character(treat)]
  master[, treatment := factor(
    fifelse(treat == "m25", "FN (m = 2.5)", "FP (m = 1.9)"),
    levels = c("FN (m = 2.5)", "FP (m = 1.9)")
  )]
  
  tr_colours <- c("FN (m = 2.5)" = "#2166AC", "FP (m = 1.9)" = "#C0392B")
  
  # ---- Summary stats for annotation --------------------------------------
  sum_tbl <- master[, .(
    pct_zero  = mean(stake == 0),
    mean_stake = mean(stake),
    med_stake  = median(stake)
  ), by = treatment]
  
  # ---- Plot: two panels ---------------------------------------------------
  # Top panel: proportion of zero stakes (bar)
  # Bottom panel: density of positive stakes only
  
  p_zero <- ggplot(sum_tbl, aes(x = treatment, y = pct_zero, fill = treatment)) +
    geom_col(width = 0.5, alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * pct_zero)),
              vjust = -0.4, size = 3.2) +
    scale_fill_manual(values = tr_colours, guide = "none") +
    scale_y_continuous(name   = "Proportion no-bet (stake = 0)",
                       labels = percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_x_discrete(name = NULL) +
    theme_classic(base_size = 10) +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          plot.margin = margin(8, 12, 4, 8))
  
  # Positive stakes only
  pos <- master[stake > 0]
  
  p_pos <- ggplot(pos, aes(x = stake, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.35, adjust = 0.8, linewidth = 0.6) +
    geom_vline(data = sum_tbl,
               aes(xintercept = mean_stake, colour = treatment),
               linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = tr_colours, name = NULL) +
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(name   = "Stake (ECU, conditional on betting)",
                       breaks = seq(0, 100, by = 20),
                       limits = c(0, 100)) +
    scale_y_continuous(name = "Density") +
    theme_classic(base_size = 10) +
    theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          legend.position    = "bottom",
          legend.text        = element_text(size = 9),
          plot.margin        = margin(4, 12, 8, 8))
  
  # ---- Combine and save --------------------------------------------------
  library(patchwork)
  combined <- p_zero / p_pos + plot_layout(heights = c(1, 2))
  
  f_out <- file.path(path_fig, "stake_distribution_by_treatment.png")
  ggsave(f_out, combined, width = 7, height = 7, dpi = 300)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}