# ============================================================
# 42_rq4-figures.R
#
# PURPOSE
#   Figures for RQ4 (side choice conditional on betting): a sequence-level
#   forest plot, split into Head-biased, neutral and Tail-biased sequences,
#   and a participant-level raincloud plot. The figures are drawn from the
#   tables of 41_rq4-tables.R and nothing is recomputed from the fit, so
#   labels agree with the tables. The reference line is the population
#   baseline hbar of the model (and parity, 0.5, for comparison), not a
#   normative benchmark.
#
# SAMPLE
#   full           all participants of each treatment in cfg$run$treatment
#
# INPUT (data/output/csv/)
#   rq4-<tr>-full-sequences.csv
#   rq4-<tr>-full-participants.csv
#   rq4-<tr>-full-model-summary.csv          (baseline hbar)
#
# OUTPUT
#   data/output/figures/rq4-<tr>-full-sequences-forest.png
#   data/output/figures/rq4-<tr>-full-participants-raincloud.png
#   data/output/csv/rq4-<tr>-full-forest-data.csv
#
# NOTES
#   - Sequence colours show the direction class of the sequence table
#     (Head-bias in red, Tail-bias in blue, neutral in grey)
#   - The horizontal axis of the raincloud covers the 2nd to 98th percentile
#     of the participant means with a padding of 5 percentage points; points
#     beyond it are drawn at the boundary
#   - Output files are overwritten at every run
# ============================================================

rq4_figures <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  sample <- "full"

  direction_colors <- c(strong_head = "#C0392B", moderate_head = "#E8A090", weak_head = "#F5C9C2",
                        neutral = "#AAAAAA",
                        weak_tail = "#C2D8F0", moderate_tail = "#74ADD1", strong_tail = "#2166AC")
  direction_labels <- c(strong_head = "Strong Head bias", moderate_head = "Moderate Head bias",
                        weak_head = "Weak Head bias", neutral = "Neutral",
                        weak_tail = "Weak Tail bias", moderate_tail = "Moderate Tail bias",
                        strong_tail = "Strong Tail bias")
  direction_levels <- names(direction_colors)

  for (tr in tr_vec) {

    stem  <- paste0("rq4-", tr, "-", sample)
    f_seq <- file.path(path_csv, paste0(stem, "-sequences.csv"))
    f_pid <- file.path(path_csv, paste0(stem, "-participants.csv"))
    f_mod <- file.path(path_csv, paste0(stem, "-model-summary.csv"))
    stopifnot(file.exists(f_seq), file.exists(f_pid), file.exists(f_mod))

    seq_tbl <- fread(f_seq, encoding = "UTF-8", colClasses = c(sequence = "character"))
    pid_tbl <- fread(f_pid, encoding = "UTF-8", colClasses = c(pid = "character"))
    mod_tbl <- fread(f_mod, encoding = "UTF-8")

    hb <- mod_tbl[parameter == "hbar"]
    stopifnot(nrow(hb) == 1L)
    hbar_val <- hb$median; hbar_lo <- hb$q025; hbar_hi <- hb$q975

    seq_tbl[, direction_label := factor(direction_label, levels = direction_levels)]

    # ---- Forest data (for use in combined plots) ----
    forest_data <- seq_tbl[, .(
      sequence = sequence, mu_mean = mu_h_mean, mu_median = mu_h_median,
      mu_q025 = mu_h_q025, mu_q975 = mu_h_q975, direction_label = direction_label,
      hbar = hbar_val, hbar_lo = hbar_lo, hbar_hi = hbar_hi,
      treatment = tr, sample = sample, rq = "rq4"
    )]
    f_forest_data <- file.path(path_csv, paste0(stem, "-forest-data.csv"))
    fwrite(forest_data, f_forest_data)
    msg("Saved: ", f_forest_data)

    # ---- Forest plot: Head-biased | neutral | Tail-biased ----
    dt <- copy(seq_tbl)
    dt[, sequence := factor(sequence, levels = dt[order(mu_h_mean), sequence])]
    dt[, panel := fcase(
      direction_label %in% c("strong_head", "moderate_head", "weak_head"), "Head-biased",
      direction_label == "neutral", "Neutral",
      direction_label %in% c("weak_tail", "moderate_tail", "strong_tail"), "Tail-biased")]
    dt[, panel := factor(panel, levels = c("Head-biased", "Neutral", "Tail-biased"))]

    x_lo <- floor(min(dt$mu_h_q025) * 10) / 10 - 0.03
    x_hi <- ceiling(max(dt$mu_h_q975) * 10) / 10 + 0.03

    p_forest <- ggplot(dt, aes(y = sequence)) +
      annotate("rect", xmin = hbar_lo, xmax = hbar_hi, ymin = -Inf, ymax = Inf,
               fill = "steelblue", alpha = 0.08) +
      geom_vline(xintercept = hbar_val, colour = "steelblue", linewidth = 0.65) +
      geom_vline(xintercept = 0.5, colour = "grey40", linewidth = 0.45, linetype = "dashed") +
      geom_segment(aes(x = mu_h_q025, xend = mu_h_q975, yend = sequence, colour = direction_label),
                   linewidth = 0.65, alpha = 0.90) +
      geom_point(aes(x = mu_h_mean, colour = direction_label), size = 1.8) +
      scale_colour_manual(values = direction_colors, labels = direction_labels, name = NULL, drop = FALSE) +
      scale_x_continuous(name = expression(hat(mu)[s]^h ~ "(posterior mean Heads probability)"),
                         limits = c(x_lo, x_hi), breaks = seq(0, 1, by = 0.1),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_y_discrete(name = NULL) +
      annotate("text", x = hbar_val, y = Inf, label = sprintf("bar(h) == %.2f", hbar_val),
               parse = TRUE, hjust = -0.10, vjust = 1.5, size = 2.8, colour = "steelblue") +
      annotate("text", x = 0.5, y = Inf, label = "Parity (0.50)",
               hjust = 1.08, vjust = 1.5, size = 2.8, colour = "grey40") +
      facet_wrap(~ panel, ncol = 3, scales = "free_y") +
      theme_classic(base_size = 10) +
      theme(axis.text.y = element_text(size = 7, family = "mono"),
            axis.title.x = element_text(size = 9),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "bottom", legend.text = element_text(size = 8),
            strip.text = element_text(size = 9, face = "bold"), strip.background = element_blank(),
            plot.margin = margin(10, 20, 8, 8)) +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                   override.aes = list(size = 2.5, linewidth = 0.8)))

    max_n <- max(dt[, .(n = uniqueN(sequence)), by = panel]$n)
    f_forest <- file.path(path_fig, paste0(stem, "-sequences-forest.png"))
    ggsave(f_forest, p_forest, width = 14, height = max(4, max_n * 0.28 + 2), dpi = 300)
    msg("Saved: ", f_forest)

    # ---- Participant raincloud ----
    x_lo_p <- max(0, quantile(pid_tbl$mu_h_mean, 0.02) - 0.05)
    x_hi_p <- min(1, quantile(pid_tbl$mu_h_mean, 0.98) + 0.05)

    p_rain <- ggplot(pid_tbl, aes(x = mu_h_mean)) +
      annotate("rect", xmin = hbar_lo, xmax = hbar_hi, ymin = -Inf, ymax = Inf,
               fill = "steelblue", alpha = 0.10) +
      geom_vline(xintercept = hbar_val, colour = "steelblue", linewidth = 0.70) +
      geom_vline(xintercept = 0.5, colour = "grey55", linewidth = 0.45, linetype = "longdash") +
      stat_slab(aes(y = 0.35), fill = "#74ADD1", colour = NA, adjust = 2.0, height = 0.40,
                alpha = 0.82, normalize = "panels", trim = FALSE) +
      geom_jitter(aes(y = 0.20), fill = "#2166AC", shape = 21, size = 1.7, alpha = 0.72,
                  stroke = 0.25, colour = "white", height = 0.055, width = 0) +
      annotate("text", x = hbar_val, y = Inf, label = sprintf("bar(h) == %.2f", hbar_val),
               parse = TRUE, hjust = -0.08, vjust = 1.5, size = 2.9, colour = "steelblue") +
      scale_x_continuous(name = expression(hat(mu)[i]^h ~ "(posterior mean Heads probability)"),
                         limits = c(x_lo_p, x_hi_p), oob = scales::squish,
                         breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(name = NULL, breaks = NULL, expand = expansion(mult = c(0.02, 0.12))) +
      theme_classic(base_size = 11) +
      theme(axis.title.x = element_text(size = 10), axis.text.y = element_blank(),
            axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "none", plot.margin = margin(14, 16, 6, 8))

    f_rain <- file.path(path_fig, paste0(stem, "-participants-raincloud.png"))
    ggsave(f_rain, p_rain, width = 8, height = 3.2, dpi = 300)
    msg("Saved: ", f_rain)
  }

  invisible(TRUE)
}
