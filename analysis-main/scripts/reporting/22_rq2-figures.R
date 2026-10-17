# ============================================================
# 22_rq2-figures.R
#
# PURPOSE
#   Figures for RQ2 (stake calibration): a sequence-level forest plot and
#   a participant-level raincloud plot for each treatment and sample.
#   The figures are drawn from the tables of 21_rq2-tables.R and nothing is
#   recomputed from the fit, so labels and classes agree with the tables.
#   The data of the forest plot are also exported for the combined forest
#   plot of several research questions.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#
# INPUT (data/output/csv/)
#   rq2-<tr>-<sample>-sequences.csv
#   rq2-<tr>-<sample>-participants.csv
#   rq2-<tr>-<sample>-model-summary.csv     (grand mean stake deviation)
#
# OUTPUT
#   data/output/figures/rq2-<tr>-<sample>-sequences-forest.png
#   data/output/figures/rq2-<tr>-<sample>-participants-raincloud.png
#   data/output/csv/rq2-<tr>-<sample>-forest-data.csv
#
# NOTES
#   - Sequence colours show the grand-mean class (above / likely above /
#     neutral / likely below / below) of the sequence table
#   - The horizontal axis is the stake deviation from the EU-optimal stake
#     as a proportion of the endowment (0 = optimal)
#   - Output files are overwritten at every run
# ============================================================

rq2_figures <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  figure_set <- function(tr, sample) {

    stem  <- paste0("rq2-", tr, "-", sample)
    f_seq <- file.path(path_csv, paste0(stem, "-sequences.csv"))
    f_pid <- file.path(path_csv, paste0(stem, "-participants.csv"))
    f_mod <- file.path(path_csv, paste0(stem, "-model-summary.csv"))
    stopifnot(file.exists(f_seq), file.exists(f_pid), file.exists(f_mod))

    seq_tbl <- fread(f_seq, encoding = "UTF-8", colClasses = c(sequence = "character"))
    pid_tbl <- fread(f_pid, encoding = "UTF-8", colClasses = c(pid = "character"))
    mod_tbl <- fread(f_mod, encoding = "UTF-8")

    gm <- mod_tbl[parameter == "grand_mean_a"]
    stopifnot(nrow(gm) == 1L)
    grand_mean <- gm$median
    grand_lo   <- gm$q025
    grand_hi   <- gm$q975

    # ---- Forest data, for the combined plot of several research questions ----
    forest_data <- seq_tbl[, .(
      sequence    = sequence,
      mu_mean     = mu_a_mean,
      mu_median   = mu_a_median,
      mu_q025     = mu_a_q025,
      mu_q975     = mu_a_q975,
      calib_label = calib_label,
      grand_label = grand_label,
      grand_mean  = grand_mean,
      grand_lo    = grand_lo,
      grand_hi    = grand_hi,
      treatment   = tr,
      sample      = sample,
      rq          = "rq2"
    )]
    f_forest_data <- file.path(path_csv, paste0(stem, "-forest-data.csv"))
    fwrite(forest_data, f_forest_data)
    msg("Saved: ", f_forest_data)

    # ---- Forest plot ----
    dt <- copy(seq_tbl)
    dt[, sequence    := factor(sequence, levels = dt[order(mu_a_mean), sequence])]
    dt[, grand_label := factor(grand_label,
                               levels = c("above", "likely_above", "neutral",
                                          "likely_below", "below"))]

    label_colors <- c(above = "#C0392B", likely_above = "#E8A090", neutral = "#AAAAAA",
                      likely_below = "#90B8D8", below = "#2166AC")
    label_names  <- c(above = "Above grand mean", likely_above = "Likely above grand mean",
                      neutral = "Neutral", likely_below = "Likely below grand mean",
                      below = "Below grand mean")

    # Label sides: "EU-optimal" goes left of the zero line and the grand-mean label
    # right of its line, unless there is little room on that side of the plot
    x_rng   <- range(c(dt$mu_a_q025, dt$mu_a_q975, 0))
    room    <- function(x, side) (if (side == "left") x - x_rng[1] else x_rng[2] - x) / diff(x_rng)
    h_eu    <- if (room(0, "left") < 0.15) -0.05 else 1.05
    h_grand <- if (room(grand_mean, "right") < 0.15) 1.05 else -0.05
    
    p <- ggplot(dt, aes(y = sequence)) +
      annotate("rect", xmin = grand_lo, xmax = grand_hi, ymin = -Inf, ymax = Inf,
               fill = "steelblue", alpha = 0.10) +
      geom_vline(xintercept = grand_mean, colour = "steelblue", linewidth = 0.6) +
      geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5, linetype = "dashed") +
      geom_segment(aes(x = mu_a_q025, xend = mu_a_q975, y = sequence, yend = sequence,
                       colour = grand_label), linewidth = 0.5, alpha = 0.7) +
      geom_point(aes(x = mu_a_mean, colour = grand_label), size = 1.8) +
      scale_colour_manual(values = label_colors, name = NULL, labels = label_names) +
      scale_x_continuous(
        name   = expression(hat(mu)[s]^a ~ "(posterior mean proportional stake deviation)"),
        labels = scales::percent_format(accuracy = 1)) +
      scale_y_discrete(name = NULL, expand = expansion(add = c(0.6, 2.4))) +
      annotate("text", x = grand_mean, y = nrow(dt) + 1.6,
               label = sprintf("Grand mean = %.1f%%", 100 * grand_mean),
               hjust = h_grand, size = 2.8, colour = "steelblue") +
      annotate("text", x = 0, y = nrow(dt) + 1.6, label = "EU-optimal",
               hjust = h_eu, size = 2.8, colour = "grey30") +
      theme_classic(base_size = 10) +
      theme(axis.text.y        = element_text(size = 7, family = "mono"),
            axis.title.x       = element_text(size = 9),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position    = "bottom",
            legend.text        = element_text(size = 8),
            plot.margin        = margin(8, 12, 8, 8))

    f_forest <- file.path(path_fig, paste0(stem, "-sequences-forest.png"))
    ggsave(f_forest, p, width = 7, height = 14, dpi = 300)
    msg("Saved: ", f_forest)

    # ---- Participant raincloud ----
    p <- ggplot(pid_tbl, aes(x = mu_a_mean, y = 1)) +
      geom_vline(xintercept = grand_mean, colour = "steelblue", linewidth = 0.6) +
      geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5, linetype = "dashed") +
      annotate("text", x = grand_mean, y = 1.75,
               label = sprintf("Grand mean = %.1f%%", 100 * grand_mean),
               hjust = -0.05, size = 2.8, colour = "steelblue") +
      annotate("text", x = 0, y = 1.75, label = "EU-optimal",
               hjust = 1.05, size = 2.8, colour = "grey30") +
      expand_limits(y = 1.85) +
      stat_halfeye(fill = "#74ADD1", adjust = 0.8, width = 0.5, .width = 0,
                   point_colour = NA, alpha = 0.80, position = position_nudge(y = 0.12)) +
      geom_jitter(fill = "#2166AC", shape = 21, size = 1.8, alpha = 0.75, stroke = 0.25,
                  colour = "white", height = 0.07, width = 0) +
      scale_x_continuous(
        name   = expression(hat(mu)[i]^a ~ "(posterior mean proportional stake deviation)"),
        labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(name = NULL, breaks = NULL) +
      theme_classic(base_size = 11) +
      theme(axis.title.x       = element_text(size = 10),
            axis.text.y        = element_blank(),
            axis.line.y        = element_blank(),
            axis.ticks.y       = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position    = "none",
            plot.margin        = margin(8, 12, 8, 8))

    f_rain <- file.path(path_fig, paste0(stem, "-participants-raincloud.png"))
    ggsave(f_rain, p, width = 8, height = 3, dpi = 300)
    msg("Saved: ", f_rain)

    invisible(NULL)
  }

  for (tr in tr_vec) figure_set(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    figure_set(tr, "confirmatory")
  }

  invisible(TRUE)
}
