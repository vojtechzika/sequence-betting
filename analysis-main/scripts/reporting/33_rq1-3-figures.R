# ============================================================
# 33_rq1-3-figures.R
#
# PURPOSE
#   Combined figures for RQ1 (betting), RQ2 (stake deviation) and RQ3
#   (welfare loss), drawn from the tables of 11_rq1-tables.R,
#   21_rq2-tables.R and 31_rq3-tables.R; nothing is recomputed from the fits.
#     1. Sequence forest plot with three panels (RQ3 | RQ1 | RQ2), sequences
#        ordered by the posterior mean welfare loss. Sequence labels are
#        coloured by the number of the RQ1 and RQ2 grand-mean classes that
#        are flagged (above, likely above, below, likely below).
#     2. Scatter companion: the deviations of the RQ1 and RQ2 sequence means
#        from their grand means against the rank of the sequence in welfare
#        loss, with an OLS trend fitted to the non-pure sequences. The OLS
#        slopes are exported.
#     3. Participant raincloud with three panels (RQ1, RQ2, RQ3) for the
#        confirmatory sample, with the full sample as a grey background
#        distribution.
#
# SAMPLES
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE; the full sample enters the raincloud as background
#
# INPUT (data/output/csv/)
#   rq<k>-<tr>-<sample>-forest-data.csv      k = 1, 2, 3
#   rq<k>-<tr>-<sample>-participants.csv     k = 1, 2, 3
#   rq<k>-<tr>-<sample>-model-summary.csv    grand means (grand_mean_bet,
#                                            grand_mean_a, grand_mean_c)
#
# OUTPUT
#   data/output/figures/rq1-3-<tr>-confirmatory-sequences-forest.png
#   data/output/figures/rq1-3-<tr>-confirmatory-sequences-scatter.png
#   data/output/figures/rq1-3-<tr>-confirmatory-participants-raincloud.png
#   data/output/csv/rq1-3-<tr>-confirmatory-scatter-trend.csv
#
# NOTES
#   - The sequences of the three tables must coincide; the script stops
#     otherwise
#   - Pure sequences (one symbol repeated) are labelled in the scatter plot
#     and excluded from the trend lines
#   - The density of the raincloud is bounded at 0 for the welfare loss
#     (requires ggdist >= 3.3.0)
#   - Output files are overwritten at every run
# ============================================================

rq1_3_figures <- function(cfg) {

  suppressPackageStartupMessages({
    library(gridExtra); library(grid); library(scales); library(ggdist); library(patchwork)
  })

  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))
  tr_vec <- tr_vec[vapply(tr_vec, function(tr) isTRUE(design$a_flags$betting_normative[[tr]]), logical(1))]

  read_csv <- function(rq, tr, sample, what, cls = NULL) {
    f <- file.path(path_csv, paste0(rq, "-", tr, "-", sample, "-", what, ".csv"))
    stopifnot(file.exists(f))
    fread(f, encoding = "UTF-8", colClasses = cls)
  }

  grand_levels  <- c("above", "likely_above", "neutral", "likely_below", "below")
  grand_colours <- c(above = "#C0392B", likely_above = "#E8A090", neutral = "#999999",
                     likely_below = "#90B8D8", below = "#2166AC")
  grand_names   <- c(above = "Above grand mean", likely_above = "Likely above",
                     neutral = "Equal", likely_below = "Likely below", below = "Below grand mean")

  # ============================================================
  # Forest plot and scatter companion
  # ============================================================
  forest_set <- function(tr, sample) {

    cls <- c(sequence = "character")
    rq1 <- read_csv("rq1", tr, sample, "forest-data", cls)
    rq2 <- read_csv("rq2", tr, sample, "forest-data", cls)
    rq3 <- read_csv("rq3", tr, sample, "forest-data", cls)
    stopifnot(setequal(rq1$sequence, rq3$sequence), setequal(rq2$sequence, rq3$sequence))

    flagged <- function(x) x %in% c("above", "likely_above", "below", "likely_below")
    flag_tbl <- data.table(
      sequence = rq3$sequence, mu_mean = rq3$mu_mean,
      n_flagged = flagged(rq1$grand_label[match(rq3$sequence, rq1$sequence)]) +
                  flagged(rq2$grand_label[match(rq3$sequence, rq2$sequence)])
    )
    setorder(flag_tbl, -mu_mean)
    seq_order  <- flag_tbl$sequence
    seq_levels <- rev(seq_order)                       # highest loss at the top
    label_col  <- c("0" = "grey60", "1" = "navyblue", "2" = "red")
    y_colours  <- label_col[as.character(flag_tbl$n_flagged[match(seq_levels, flag_tbl$sequence)])]

    prep <- function(dt) {
      dt <- copy(dt)
      dt[, sequence    := factor(sequence, levels = seq_levels)]
      dt[, grand_label := factor(grand_label, levels = grand_levels)]
      dt
    }
    rq1 <- prep(rq1); rq2 <- prep(rq2); rq3 <- prep(rq3)

    make_panel <- function(dt, x_label, show_y = FALSE) {
      gm <- dt$grand_mean[1]; lo <- dt$grand_lo[1]; hi <- dt$grand_hi[1]
      p <- ggplot(dt, aes(y = sequence)) +
        annotate("rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
                 fill = "steelblue", alpha = 0.08) +
        geom_vline(xintercept = gm, colour = "steelblue", linewidth = 0.55) +
        geom_segment(aes(x = mu_q025, xend = mu_q975, yend = sequence, colour = grand_label),
                     linewidth = 0.45, alpha = 0.7) +
        geom_point(aes(x = mu_mean, colour = grand_label), size = 1.8) +
        scale_colour_manual(values = grand_colours, labels = grand_names, drop = FALSE) +
        scale_x_continuous(name = x_label, labels = percent_format(accuracy = 1)) +
        annotate("text", x = gm, y = nlevels(dt$sequence) + 0.5, label = sprintf("%.1f%%", 100 * gm),
                 hjust = -0.1, size = 2.5, colour = "steelblue") +
        theme_classic(base_size = 10) +
        theme(axis.title.x = element_text(size = 8.5), axis.title.y = element_blank(),
              axis.ticks.y = element_blank(), legend.position = "none",
              panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
              plot.margin = margin(4, 8, 4, 4))
      if (show_y) {
        p + scale_y_discrete(limits = seq_levels) +
          theme(axis.text.y = element_text(family = "mono", face = "bold", size = 7, hjust = 1,
                                           colour = y_colours, margin = margin(r = 2)))
      } else p + theme(axis.text.y = element_blank())
    }

    p3 <- make_panel(rq3, expression(hat(mu)[s]^c ~ "(Welfare loss, share of endowment)"), TRUE)
    p1 <- make_panel(rq1, expression(hat(mu)[s]^b ~ "(Betting rate)"))
    p2 <- make_panel(rq2, expression(hat(mu)[s]^a ~ "(Stake deviation, share of endowment)"))

    legend_plot <- ggplot(data.table(x = 1:5, y = 0, label = factor(grand_levels, levels = grand_levels)),
                          aes(x = x, y = y, colour = label)) +
      geom_point(size = 4) +
      scale_colour_manual(values = grand_colours, labels = grand_names, name = "Grand-mean classification") +
      guides(colour = guide_legend(nrow = 1, title.hjust = 0.5, override.aes = list(size = 4))) +
      coord_cartesian(xlim = c(10, 20)) + theme_void() +
      theme(legend.position = "bottom", legend.direction = "horizontal",
            legend.text = element_text(size = 10), legend.title = element_text(size = 10.5, face = "bold"),
            plot.margin = margin(0, 0, 0, 0))

    aligned <- cbind(ggplotGrob(p3), ggplotGrob(p1), ggplotGrob(p2), size = "max")   # cbind.gtable
    final <- arrangeGrob(aligned, ggplotGrob(legend_plot), ncol = 1,
                         heights = unit(c(1, 0.06), c("null", "npc")))
    stem <- paste0("rq1-3-", tr, "-", sample)
    f_forest <- file.path(path_fig, paste0(stem, "-sequences-forest.png"))
    ggsave(f_forest, plot = final, width = 17, height = 18, dpi = 300)
    msg("Saved: ", f_forest)

    # ---- Scatter companion: deviations from the grand means against loss rank ----
    dt <- data.table(
      sequence  = seq_order, rank_loss = seq_along(seq_order),
      bet_dev   = rq1$mu_mean[match(seq_order, as.character(rq1$sequence))] - rq1$grand_mean[1],
      stake_dev = rq2$mu_mean[match(seq_order, as.character(rq2$sequence))] - rq2$grand_mean[1]
    )
    dl <- melt(dt, id.vars = c("sequence", "rank_loss"), measure.vars = c("bet_dev", "stake_dev"),
               variable.name = "margin", value.name = "deviation")
    dl[, margin := factor(margin, levels = c("bet_dev", "stake_dev"),
                          labels = c("Betting-rate deviation", "Stake deviation"))]
    dl[, is_pure := grepl("^(.)\\1*$", sequence)]
    margin_colours <- c("Betting-rate deviation" = "#2166AC", "Stake deviation" = "#C0392B")

    dtrend <- dl[is_pure == FALSE]
    trend <- rbindlist(lapply(levels(dtrend$margin), function(mg) {
      fit <- lm(deviation ~ rank_loss, data = dtrend[margin == mg])
      co  <- summary(fit)$coefficients
      data.table(margin = mg, n = nobs(fit), slope = co["rank_loss", 1], se = co["rank_loss", 2],
                 t_value = co["rank_loss", 3], p_value = co["rank_loss", 4],
                 r_squared = summary(fit)$r.squared)
    }))
    f_trend <- file.path(path_csv, paste0(stem, "-scatter-trend.csv"))
    fwrite(trend, f_trend); msg("Saved: ", f_trend)

    lab_data <- dl[is_pure == TRUE]
    lab_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = lab_data, aes(label = sequence), colour = "black", size = 2.8,
                               family = "mono", fontface = "bold", show.legend = FALSE,
                               min.segment.length = 0, seed = 1)
    } else geom_text(data = lab_data, aes(label = sequence), colour = "black", size = 2.8,
                     family = "mono", fontface = "bold", show.legend = FALSE, vjust = -0.8)

    ps <- ggplot(dl, aes(x = rank_loss, y = deviation, colour = margin)) +
      geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
      geom_point(aes(shape = is_pure, size = is_pure), alpha = 0.75) +
      geom_smooth(data = dtrend, aes(fill = margin), method = "lm", formula = y ~ x, se = TRUE,
                  linewidth = 0.9, alpha = 0.12) +
      lab_layer +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
      scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 3.2), guide = "none") +
      scale_colour_manual(values = margin_colours, name = NULL) +
      scale_fill_manual(values = margin_colours, guide = "none") +
      scale_x_continuous(name = "Sequence rank (1 = highest welfare loss)") +
      scale_y_continuous(name = "Deviation from grand mean", labels = percent_format(accuracy = 1)) +
      theme_classic(base_size = 10) +
      theme(legend.position = "bottom", axis.title = element_text(size = 9),
            plot.margin = margin(4, 8, 4, 4))
    f_sc <- file.path(path_fig, paste0(stem, "-sequences-scatter.png"))
    ggsave(f_sc, ps, width = 20, height = 10, units = "cm", dpi = 300)
    msg("Saved: ", f_sc)
  }

  # ============================================================
  # Participant raincloud
  # ============================================================
  raincloud_set <- function(tr) {

    cls <- c(pid = "character")
    get_pid <- function(rq, sample) read_csv(rq, tr, sample, "participants", cls)
    get_gm  <- function(rq, par, sample) {
      m <- read_csv(rq, tr, sample, "model-summary")[parameter == par]
      stopifnot(nrow(m) == 1L)
      list(median = m$median, lo = m$q025, hi = m$q975)
    }

    spec <- list(
      list(rq = "rq1", x = "mu_i_median", par = "grand_mean_bet", ref = 1,
           lab = expression(hat(mu)[i]^b ~ "(posterior median betting probability)"),
           lim = c(-0.05, 1.05), brk = seq(0, 1, by = 0.1), bounded = FALSE),
      list(rq = "rq2", x = "mu_a_median", par = "grand_mean_a", ref = 0,
           lab = expression(hat(mu)[i]^a ~ "(posterior median proportional stake deviation)"),
           lim = NULL, brk = NULL, bounded = FALSE),
      list(rq = "rq3", x = "mu_c_median", par = "grand_mean_c", ref = 0,
           lab = expression(hat(mu)[i]^c ~ "(posterior median welfare loss, share of endowment)"),
           lim = NULL, brk = NULL, bounded = TRUE)
    )

    if (utils::packageVersion("ggdist") < "3.3.0")
      stop("33_rq1-3-combined.R requires ggdist >= 3.3.0 (bounded densities).", call. = FALSE)
    col_conf <- "#74ADD1"; col_point <- "#2166AC"; col_full <- "grey60"; col_gm <- "steelblue"

    bet <- get_pid("rq1", "confirmatory")[, .(pid, mu_b_median = mu_i_median)]

    panels <- lapply(spec, function(s) {
      conf <- get_pid(s$rq, "confirmatory"); full <- get_pid(s$rq, "full")
      if (s$rq == "rq3") conf <- merge(conf[, setdiff(names(conf), "mu_b_median"), with = FALSE], bet, by = "pid", all.x = TRUE)
      gc <- get_gm(s$rq, s$par, "confirmatory"); gf <- get_gm(s$rq, s$par, "full")
      dens <- if (s$bounded) ggdist::density_bounded(bounds = c(0, Inf)) else ggdist::density_unbounded()

      p <- ggplot(conf, aes(x = .data[[s$x]], y = 1)) +
        stat_halfeye(data = full, aes(x = .data[[s$x]], y = 1), fill = col_full, adjust = 0.8,
                     width = 0.5, .width = 0, point_colour = NA, alpha = 0.40,
                     position = position_nudge(y = 0.12), inherit.aes = FALSE, density = dens) +
        geom_vline(xintercept = gf$median, colour = col_full, linewidth = 0.5) +
        geom_vline(xintercept = gc$median, colour = col_gm, linewidth = 0.6) +
        annotate("text", x = gf$median, y = Inf, label = sprintf("%.1f%%", 100 * gf$median),
                 hjust = -0.1, vjust = 1, size = 2.8, colour = col_full) +
        annotate("text", x = gc$median, y = Inf, label = sprintf("%.1f%%", 100 * gc$median),
                 hjust = 1.1, vjust = 1, size = 2.8, colour = col_gm) +
        stat_halfeye(fill = col_conf, adjust = 0.8, width = 0.5, .width = 0, point_colour = NA,
                     alpha = 0.80, position = position_nudge(y = 0.12), density = dens)
      p <- p + if (s$rq == "rq3") {
        geom_jitter(aes(fill = mu_b_median), shape = 21, size = 1.8, alpha = 0.75, stroke = 0.25,
                    colour = "white", height = 0.07, width = 0)
      } else {
        geom_jitter(fill = col_point, shape = 21, size = 1.8, alpha = 0.75, stroke = 0.25,
                    colour = "white", height = 0.07, width = 0)
      }
      p <- p + scale_x_continuous(name = s$lab, limits = s$lim, breaks = if (is.null(s$brk)) waiver() else s$brk,
                                  labels = percent_format(accuracy = 1), expand = c(0, 0)) +
        scale_y_continuous(name = NULL, breaks = NULL) +
        theme_classic(base_size = 11) +
        theme(axis.title.x = element_text(size = 10), axis.text.y = element_blank(),
              axis.line.y = element_blank(), axis.ticks.y = element_blank(),
              panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
              legend.position = "none")
      if (s$rq == "rq3")
        p <- p + scale_fill_gradient(low = "#E31A1C", high = "#2166AC", name = "Mean betting probability",
                                     labels = percent_format(accuracy = 1)) +
          theme(legend.position = "right")
      p
    })

    combined <- (panels[[1]] + theme(plot.margin = margin(8, 8, 8, 0))) +
      (panels[[2]] + theme(plot.margin = margin(8, 0, 8, 8))) +
      (panels[[3]] + theme(plot.margin = margin(8, 8, 8, 8))) +
      plot_layout(design = "AB\nCC") +
      plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
    f <- file.path(path_fig, paste0("rq1-3-", tr, "-confirmatory-participants-raincloud.png"))
    ggsave(f, combined, width = 8, height = 5, dpi = 600)
    msg("Saved: ", f)
  }

  for (tr in tr_vec) {
    forest_set(tr, "confirmatory")
    raincloud_set(tr)
  }

  invisible(NULL)
}
