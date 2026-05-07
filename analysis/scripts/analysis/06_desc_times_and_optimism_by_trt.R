# ============================================================
# 04_desc_times_and_optimism_by_trt.R
#
# PURPOSE
#   Combined figure with two panels:
#   Left:  Distribution of sequence-level mean response times
#          by treatment, with KS test.
#   Right: Distribution of participant-level LOT-R optimism scores
#          by treatment, with KS test.
#
# INPUT
#   path_src/master_sequences.csv  -- for response times (screen_ms)
#   path_src/participants.csv      -- for treatment assignment
#   path_out/lotr_scored.csv       -- for LOT-R scores
#
# OUTPUT
#   path_fig/desc_times_and_optimism_by_trt.png
# ============================================================
desc_times_and_optimism_by_trt <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  
  tr_colours <- c("FN (m = 2.5)" = "#2166AC", "FP (m = 1.9)" = "#C0392B")
  
  treat_label <- function(tr) factor(
    fifelse(tr == "m25", "FN (m = 2.5)", "FP (m = 1.9)"),
    levels = c("FP (m = 1.9)", "FN (m = 2.5)")
  )
  
  # ---- Left panel: sequence-level response times --------------------------
  master <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  master[, treat     := as.character(treat)]
  master[, treatment := treat_label(treat)]
  master[, rt_ms := as.numeric(screen_ms)]
  master[, rt    := rt_ms / 1000]
  
  rt_tbl <- master[!is.na(rt) & is.finite(rt) & rt > 0]
  
  sum_rt <- rt_tbl[, .(mean_val = mean(rt)), by = treatment]
  sum_rt[, label := paste0(treatment, ": ", round(mean_val, 1), " s")]
  
  ks_rt <- ks.test(
    rt_tbl[treatment == "FN (m = 2.5)", rt_ms],
    rt_tbl[treatment == "FP (m = 1.9)", rt_ms]
  )
  ks_rt_label <- paste0("KS: D = ", round(ks_rt$statistic, 2),
                        ifelse(ks_rt$p.value < 0.001, ", p < 0.001",
                               paste0(", p = ", round(ks_rt$p.value, 2))))
  
  p_rt <- ggplot(rt_tbl, aes(x = rt, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.4, adjust = 0.8, linewidth = 0, trim = TRUE) +
    geom_vline(data = sum_rt,
               aes(xintercept = mean_val, colour = treatment),
               linewidth = 0.7, linetype = "solid") +
    geom_text(data = sum_rt[treatment == "FN (m = 2.5)"],
              aes(x = mean_val, y = Inf, label = label, colour = treatment),
              hjust = -0.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    geom_text(data = sum_rt[treatment == "FP (m = 1.9)"],
              aes(x = mean_val, y = Inf, label = label, colour = treatment),
              hjust = 1.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    scale_fill_manual(values = tr_colours, name = NULL) +
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(
      name   = paste0("Response Time in Seconds\n(", ks_rt_label, ")"),
      breaks = seq(0, 30, by = 2)
    ) +
    coord_cartesian(xlim = c(0, 30)) +
    scale_y_continuous(name = "Density") +
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position    = "none",
      plot.margin        = margin(4, 12, 8, 8)
    )
  
  # ---- Right panel: LOT-R optimism scores ---------------------------------
  lotr <- fread(file.path(path_out, "lotr_scored.csv"), encoding = "UTF-8")
  lotr[, pid := as.character(pid)]
  
  participants <- fread(file.path(path_src, "participants.csv"), encoding = "UTF-8")
  participants[, pid   := as.character(pid)]
  participants[, treat := as.character(treat)]
  
  lotr <- merge(lotr, participants[, .(pid, treat)], by = "pid", all.x = TRUE)
  lotr[, treatment := treat_label(treat)]
  lotr <- lotr[!is.na(lotr_score) & !is.na(treatment)]
  
  sum_lotr <- lotr[, .(mean_val = mean(lotr_score)), by = treatment]
  sum_lotr[, label := paste0(treatment, ": ", round(mean_val, 1))]
  
  ks_lotr <- ks.test(
    lotr[treatment == "FN (m = 2.5)", lotr_score],
    lotr[treatment == "FP (m = 1.9)", lotr_score]
  )
  ks_lotr_label <- paste0("KS: D = ", round(ks_lotr$statistic, 2),
                          ", p = ", round(ks_lotr$p.value, 2))
  

  
  p_lotr <- ggplot(lotr, aes(x = lotr_score, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.4, adjust = 0.8, linewidth = 0) +
    geom_vline(data = sum_lotr,
               aes(xintercept = mean_val, colour = treatment),
               linewidth = 0.7, linetype = "solid") +
    geom_text(data = sum_lotr[treatment == "FN (m = 2.5)"],
              aes(x = mean_val, y = Inf, label = label, colour = treatment),
              hjust = 1.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    geom_text(data = sum_lotr[treatment == "FP (m = 1.9)"],
              aes(x = mean_val, y = Inf, label = label, colour = treatment),
              hjust = -0.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
    # annotate("text", x = Inf, y = Inf,
    #          label = ks_lotr_label,
    #          hjust = 1.05, vjust = 2.5,
    #          size = 2.8, colour = "grey30") +
    scale_fill_manual(values = tr_colours, name = NULL) +
    scale_colour_manual(values = tr_colours, name = NULL) +
    scale_x_continuous(
                       paste0("Mean LOT-R Scores\n(", ks_lotr_label,")"),
                       breaks = seq(0, 24, by = 4),
                       limits = c(0, 24)) +
    scale_y_continuous(name = "Density") +
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position    = "none",
      plot.margin        = margin(4, 12, 8, 8)
    )
  
  # ---- Combine ------------------------------------------------------------
  combined <- p_rt + p_lotr +
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
      theme = theme(plot.margin = margin(4, 4, 4, 4))
    )
  
  f_out <- file.path(path_fig, "desc_times_and_optimism_by_trt.png")
  ggsave(f_out, combined, width = 12, height = 4, dpi = 300)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}