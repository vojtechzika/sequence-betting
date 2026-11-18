# ============================================================
# 42_rq4-figures.R
#
# PURPOSE
#   Figures for RQ4 (side choice conditional on betting): a sequence-level
#   forest plot, split into Head-biased, neutral and Tail-biased sequences,
#   and a participant-level density plot of both treatments. The figures are drawn from the
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
#   data/output/figures/rq4-full-participants-distribution.png   (both treatments, side by side)
#   data/output/csv/rq4-<tr>-full-forest-data.csv
#
# NOTES
#   - Sequence colours show the direction class of the sequence table
#     (Head-bias in red, Tail-bias in blue, neutral in grey)
#   - Treatment labels and colours: _reporting-style.R
#   - The participant plot shows the distribution of the posterior mean Heads
#     probability of the participants as a density, one panel per treatment (side by side),
#     with parity (0.5, dashed) and the population baseline hbar of the
#     treatment (solid) with its 95% credible interval (shaded band)
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

  pid_all <- list()
  hbar_all <- list()

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

    # participants and baseline of this treatment, for the combined density plot
    pid_all[[tr]]  <- data.table(treatment = tr, mu_h_mean = pid_tbl$mu_h_mean)
    hbar_all[[tr]] <- data.table(treatment = tr, hbar = hbar_val, lo = hbar_lo, hi = hbar_hi)

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
    # number of sequences in each class is added to the panel titles
    n_panel <- dt[, .N, by = panel]
    dt[, panel := factor(panel, levels = levels(panel),
                         labels = sprintf("%s (n = %d)", levels(panel),
                                          n_panel$N[match(levels(panel), n_panel$panel)]))]

    # horizontal axis limits of the forest plot (Heads probability)
    x_lo <- 0.20
    x_hi <- 0.85

    p_forest <- ggplot(dt, aes(y = sequence)) +
      annotate("rect", xmin = hbar_lo, xmax = hbar_hi, ymin = -Inf, ymax = Inf,
               fill = "steelblue", alpha = 0.08) +
      geom_vline(xintercept = hbar_val, colour = "steelblue", linewidth = 0.65) +
      geom_vline(xintercept = 0.5, colour = "grey40", linewidth = 0.45, linetype = "dashed") +
      geom_segment(aes(x = mu_h_q025, xend = mu_h_q975, yend = sequence, colour = direction_label),
                   linewidth = 0.65, alpha = 0.90) +
      geom_point(aes(x = mu_h_mean, colour = direction_label), size = 2.2) +
      scale_colour_manual(values = direction_colors, labels = direction_labels, name = NULL, drop = FALSE) +
      scale_x_continuous(name = expression(hat(mu)[s]^h ~ "(posterior mean Heads probability)"),
                         limits = c(x_lo, x_hi), breaks = seq(0.2, 0.8, by = 0.2),
                         labels = scales::percent_format(accuracy = 1)) +
      # the two pure sequences are set in bold on the vertical axis
      scale_y_discrete(name = NULL, expand = expansion(add = c(0.6, 1.8)),
                       labels = function(x) as.expression(lapply(x, function(s)
                         if (s %in% c("HHHHHH", "TTTTTT")) bquote(bold(.(s))) else bquote(plain(.(s)))))) +
      annotate("text", x = hbar_val, y = Inf, label = sprintf("bar(h) == %.2f", hbar_val),
               parse = TRUE, hjust = -0.10, vjust = 1.5, size = 3.6, colour = "steelblue") +
      facet_wrap(~ panel, ncol = 3, scales = "free_y") +
      theme_classic(base_size = 10) +
      theme(axis.text.y = element_text(size = 9, family = "mono"),
            axis.text.x = element_text(size = 9),
            axis.title.x = element_text(size = 11),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "bottom", legend.text = element_text(size = 8),
            strip.text = element_text(size = 11, face = "bold"), strip.background = element_blank(),
            plot.margin = margin(10, 20, 8, 8)) +
      guides(colour = guide_legend(nrow = 1, byrow = TRUE,
                                   override.aes = list(size = 2.5, linewidth = 0.8)))

    max_n <- max(dt[, .(n = uniqueN(sequence)), by = panel]$n)
    f_forest <- file.path(path_fig, paste0(stem, "-sequences-forest.png"))
    ggsave(f_forest, p_forest, width = 9, height = max(4, max_n * 0.25 + 1.6), dpi = 300)
    msg("Saved: ", f_forest)
  }

  # ---- Participant distribution, both treatments in one plot ----
  # treatment labels and colours come from treatment_style() (_reporting-style.R)
  style    <- treatment_style(cfg)
  pid_all  <- rbindlist(pid_all)
  hbar_all <- rbindlist(hbar_all)
  pid_all[,  treatment := factor(treatment, levels = style$treat, labels = style$type)]
  hbar_all[, treatment := factor(treatment, levels = style$treat, labels = style$type)]
  n_tr    <- pid_all[, .(n = .N), by = treatment]
  fills   <- setNames(style$colour, style$type)
  hbar_all[, lab := sprintf("bar(h) == %.2f", hbar)]

  # density of the participant means, one panel per treatment (as in the other distribution figures)
  p_dist <- ggplot(pid_all, aes(x = mu_h_mean, fill = treatment)) +
    geom_rect(data = hbar_all, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey20", alpha = 0.15) +
    geom_density(adjust = 1.2, colour = "grey20", linewidth = 0.5, alpha = 0.85) +
    geom_vline(xintercept = 0.5, colour = "grey55", linewidth = 0.45, linetype = "dashed") +
    geom_vline(data = hbar_all, aes(xintercept = hbar), colour = "grey20", linewidth = 0.7) +
    geom_text(data = hbar_all, aes(x = hbar, y = Inf, label = lab), inherit.aes = FALSE, parse = TRUE,
              hjust = -0.08, vjust = 1.6, size = 3, colour = "grey20") +
    geom_text(data = n_tr, aes(x = -Inf, y = Inf, label = paste0("n = ", n)), inherit.aes = FALSE,
              size = 3, colour = "grey35", hjust = -0.15, vjust = 1.6) +
    scale_fill_manual(values = fills, guide = "none") +
    scale_x_continuous(name = expression(hat(mu)[i]^h ~ "(posterior mean Heads probability)"),
                       limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(name = "Density", expand = expansion(mult = c(0, 0.25))) +
    facet_grid(. ~ treatment) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 10),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

  f_dist <- file.path(path_fig, paste0("rq4-", sample, "-participants-distribution.png"))
  ggsave(f_dist, p_dist, width = 9, height = 2.8, dpi = 300)
  msg("Saved: ", f_dist)

  invisible(TRUE)
}
