# ============================================================
# 07_response-times-and-optimism.R
#
# PURPOSE
#   Figure with two panels comparing the treatments.
#   Left:  distribution of response times of the sequence task (seconds,
#          all trials).
#   Right: distribution of participant-level LOT-R optimism scores.
#   Each panel marks the treatment means and reports a two-sample
#   Kolmogorov-Smirnov test. The test of response times uses participant
#   mean response times, so that repeated trials of a participant are not
#   treated as independent observations; the density shows all trials.
#
# INPUT
#   path_der/master-sequences.csv   response times (screen_ms)
#   path_src/participants.csv       treatment assignment
#   path_der/lotr-scored.csv        LOT-R scores
#
# OUTPUT
#   data/output/figures/response-times-and-optimism.png
#
# NOTES
#   - Treatment labels and colours: _reporting-style.R
#   - The test compares the first two treatments of cfg$run$treatment
#   - The response-time axis is cut at 30 seconds; longer times are not
#     shown but enter the means and the test
# ============================================================

response_times_and_optimism <- function(cfg) {
  
  style <- treatment_style(cfg)
  stopifnot(nrow(style) == 2L)
  
  tr_colours <- setNames(style$colour, style$label)
  lab_a <- style$label[1]; lab_b <- style$label[2]
  treat_label <- function(tr) factor(style$label[match(tr, style$treat)], levels = rev(style$label))
  
  ks_label <- function(x_a, x_b, note = "") {
    ks <- suppressWarnings(ks.test(x_a, x_b, exact = FALSE))
    p  <- if (ks$p.value < 0.001) ", p < 0.001" else paste0(", p = ", round(ks$p.value, 2))
    paste0("KS: D = ", round(ks$statistic, 2), p, note)
  }
  
  density_panel <- function(d, var, summary_tbl, x_scale, trim, lab_left) {
    # the treatment whose mean label is placed left of its line
    ggplot(d, aes(x = .data[[var]], fill = treatment, colour = treatment)) +
      geom_density(alpha = 0.4, adjust = 0.8, linewidth = 0, trim = trim) +
      geom_vline(data = summary_tbl, aes(xintercept = mean_val, colour = treatment),
                 linewidth = 0.7) +
      geom_text(data = summary_tbl[treatment == lab_left],
                aes(x = mean_val, y = Inf, label = label, colour = treatment),
                hjust = 1.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
      geom_text(data = summary_tbl[treatment != lab_left],
                aes(x = mean_val, y = Inf, label = label, colour = treatment),
                hjust = -0.08, vjust = 1.5, size = 2.8, show.legend = FALSE) +
      scale_fill_manual(values = tr_colours, name = NULL) +
      scale_colour_manual(values = tr_colours, name = NULL) +
      x_scale +
      scale_y_continuous(name = "Density") +
      theme_classic(base_size = 10) +
      theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "none", plot.margin = margin(4, 12, 8, 8))
  }
  
  # ---- Left panel: response times ----
  master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
  master[, pid       := as.character(pid)]
  master[, treat     := as.character(treat)]
  master[, treatment := treat_label(treat)]
  master[, rt        := as.numeric(screen_ms) / 1000]
  rt_tbl <- master[!is.na(rt) & is.finite(rt) & rt > 0]
  
  sum_rt <- rt_tbl[, .(mean_val = mean(rt)), by = treatment]
  sum_rt[, label := paste0(treatment, ": ", round(mean_val, 1), " s")]
  
  pid_rt <- rt_tbl[, .(mean_rt = mean(rt), treat = treat[1]), by = pid]
  ks_rt  <- ks_label(pid_rt[treat == style$treat[1], mean_rt],
                     pid_rt[treat == style$treat[2], mean_rt], ", participant means")
  
  p_rt <- density_panel(
    rt_tbl, "rt", sum_rt,
    list(scale_x_continuous(name = paste0("Response Time in Seconds\n(", ks_rt, ")"),
                            breaks = seq(0, 30, by = 2)),
         coord_cartesian(xlim = c(0, 30))),
    trim = TRUE, lab_left = lab_b)
  
  # ---- Right panel: LOT-R scores ----
  lotr <- fread(file.path(path_der, "lotr-scored.csv"), encoding = "UTF-8")
  lotr[, pid := as.character(pid)]
  participants <- fread(file.path(path_src, "participants.csv"), encoding = "UTF-8")
  participants[, pid   := as.character(pid)]
  participants[, treat := as.character(treat)]
  lotr <- merge(lotr, participants[, .(pid, treat)], by = "pid", all.x = TRUE)
  lotr[, treatment := treat_label(treat)]
  lotr <- lotr[!is.na(lotr_score) & !is.na(treatment)]
  
  sum_lotr <- lotr[, .(mean_val = mean(lotr_score)), by = treatment]
  sum_lotr[, label := paste0(treatment, ": ", round(mean_val, 1))]
  ks_lotr <- ks_label(lotr[treat == style$treat[1], lotr_score],
                      lotr[treat == style$treat[2], lotr_score])
  
  p_lotr <- density_panel(
    lotr, "lotr_score", sum_lotr,
    scale_x_continuous(name = paste0("LOT-R score (0 to 24)\n(", ks_lotr, ")"),
                       breaks = seq(0, 24, by = 4), limits = c(0, 24)),
    trim = FALSE, lab_left = lab_a)
  
  combined <- p_rt + p_lotr + plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
  
  f_out <- file.path(path_fig, "response-times-and-optimism.png")
  ggsave(f_out, combined, width = 12, height = 4, dpi = 300)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}
