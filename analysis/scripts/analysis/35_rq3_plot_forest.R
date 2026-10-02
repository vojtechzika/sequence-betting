# ============================================================
# 35_rq3_plot_forest.R
#
# Panel order:  RQ3 (with y-axis labels) | RQ1 | RQ2
# Sequences:    ordered by RQ3 posterior mean welfare loss (high→low)
# Label colours: 0 RQs flagged (RQ1+RQ2) = grey, 1 = #1B7837, 2 = #CC0000
# Legend:       standalone ggplot row below the panels
# No title or caption.
#
# Also produces a compact scatter companion: RQ1 (betting probability)
# and RQ2 (stake deviation) sequence-level deviations from their own
# grand means, plotted against each sequence's RANK on RQ3 welfare loss
# (1 = highest loss). Built to show the pattern described in the text --
# that betting-probability deviations track the welfare-loss ordering
# while stake deviations scatter around zero throughout -- in one
# compact panel rather than three tall columns of individual CIs.
# The two pure sequences (HHHHHH, TTTTTT) are highlighted and labelled,
# since they are the exception singled out in the text (both margins
# deviate); they are EXCLUDED from the fitted trend lines (OLS, not
# loess) so that their extreme values at the top of the welfare-loss
# ranking don't distort the trend that describes the other 62
# sequences. Intended for the main text, with the full forest plot
# above moved to the appendix as the detailed per-sequence reference.
#
# The OLS slope, SE, t, and p-value behind each trend line (fit on the
# same 62 non-pure sequences) are written to a small CSV alongside the
# figures, for reporting the numbers in the figure caption.
# ============================================================

rq1_3_forest <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(scales)
  
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  if (has_ggrepel) library(ggrepel)
  
  f_rq1 <- file.path(path_out, "rq1_m25_confirmatory_forest_data.csv")
  f_rq2 <- file.path(path_out, "rq2_m25_confirmatory_forest_data.csv")
  f_rq3 <- file.path(path_out, "rq3_m25_confirmatory_forest_data.csv")
  
  rq1 <- fread(f_rq1, encoding = "UTF-8")
  rq2 <- fread(f_rq2, encoding = "UTF-8")
  rq3 <- fread(f_rq3, encoding = "UTF-8")
  
  setnames(rq1, "underbet_label", "prereg_label")
  setnames(rq2, "calib_label",    "prereg_label")
  setnames(rq3, "loss_label",     "prereg_label")
  
  # ---- Flag table ---------------------------------------------------------
  is_flagged <- function(gl) gl %in% c("above", "likely_above", "below", "likely_below")
  
  all_seqs <- rq3[, sequence]
  
  flag_tbl <- data.table(
    sequence  = all_seqs,
    mu_mean   = rq3[match(all_seqs, sequence), mu_mean],
    flag_rq1  = is_flagged(rq1[match(all_seqs, sequence), grand_label]),
    flag_rq2  = is_flagged(rq2[match(all_seqs, sequence), grand_label])
  )
  flag_tbl[, n_flagged := flag_rq1 + flag_rq2]
  
  # ---- Ordering: RQ3 mu_mean desc -----------------------------------------
  flag_tbl  <- flag_tbl[order(-mu_mean)]
  seq_order  <- flag_tbl[, sequence]
  seq_levels <- rev(seq_order)   # rev so highest loss plots at top
  
  # ---- Label colours (RQ1+RQ2 flags only) ---------------------------------
  flag_palette <- c("0" = "grey60", "1" = "navyblue", "2" = "red")
  flag_tbl[, label_colour := flag_palette[as.character(n_flagged)]]
  
  # Named vector of colours in seq_levels order (bottom→top) for axis.text.y
  label_colours <- flag_tbl[match(seq_levels, sequence), label_colour]
  
  # ---- Factorise all panels -----------------------------------------------
  grand_levels <- c("above", "likely_above", "neutral", "likely_below", "below")
  
  factorize <- function(dt) {
    dt[, sequence    := factor(sequence,    levels = seq_levels)]
    dt[, grand_label := factor(grand_label, levels = grand_levels)]
    dt
  }
  rq1 <- factorize(rq1); rq2 <- factorize(rq2); rq3 <- factorize(rq3)
  
  grand_colours <- c(
    above        = "#C0392B",
    likely_above = "#E8A090",
    neutral      = "#999999",
    likely_below = "#90B8D8",
    below        = "#2166AC"
  )
  grand_labels_named <- c(
    above        = "Above grand mean",
    likely_above = "Likely above",
    neutral      = "Equal",
    likely_below = "Likely below",
    below        = "Below grand mean"
  )
  
  # ---- Data panel builder -------------------------------------------------
  make_panel <- function(dt, x_label,
                         show_y_labels  = FALSE,
                         y_label_colours = NULL,
                         x_label_fmt    = percent_format(accuracy = 1)) {
    
    gm    <- dt[1, grand_mean]
    gm_lo <- dt[1, grand_lo]
    gm_hi <- dt[1, grand_hi]
    n_seq <- nlevels(dt$sequence)
    
    p <- ggplot(dt, aes(y = sequence)) +
      annotate("rect",
               xmin = gm_lo, xmax = gm_hi, ymin = -Inf, ymax = Inf,
               fill = "steelblue", alpha = 0.08) +
      geom_vline(xintercept = gm, colour = "steelblue", linewidth = 0.55) +
      geom_segment(aes(x = mu_q025, xend = mu_q975,
                       y = sequence,  yend = sequence,
                       colour = grand_label),
                   linewidth = 0.45, alpha = 0.7) +
      geom_point(aes(x = mu_mean, colour = grand_label), size = 1.8) +
      scale_colour_manual(values = grand_colours, labels = grand_labels_named,
                          name = "Grand-mean classification", drop = FALSE) +
      scale_x_continuous(name = x_label, labels = x_label_fmt) +
      annotate("text",
               x = gm, y = n_seq + 0.5,
               label = sprintf("%.1f%%", 100 * gm),
               hjust = -0.1, vjust = 0.5, size = 2.5, colour = "steelblue") +
      theme_classic(base_size = 10) +
      theme(axis.title.x       = element_text(size = 8.5),
            axis.title.y       = element_blank(),
            axis.ticks.y       = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
            legend.position    = "none",
            plot.margin        = margin(4, 8, 4, 4))
    
    if (show_y_labels) {
      p <- p +
        scale_y_discrete(limits = seq_levels) +
        theme(axis.text.y = element_text(
          #colour   = y_label_colours,
          family   = "mono",
          face     = "bold",
          size     = 7,
          hjust    = 1,
          margin   = margin(r = 2)
        ))
    } else {
      p <- p + theme(axis.text.y = element_blank())
    }
    
    p
  }
  
  p3 <- make_panel(rq3,
                   x_label         = expression(hat(mu)[s]^c ~ "(Welfare loss, share of endowment)"),
                   show_y_labels   = TRUE,
                   y_label_colours = label_colours)
  p1 <- make_panel(rq1, expression(hat(mu)[s]^b ~ "(Betting rate)"))
  p2 <- make_panel(rq2, expression(hat(mu)[s]^a ~ "(Stake deviation, share of endowment)"))
  
  # ---- Legend as its own ggplot -------------------------------------------
  legend_plot <- ggplot(
    data.table(x     = 1:5,
               y     = 0,
               label = factor(grand_levels, levels = grand_levels)),
    aes(x = x, y = y, colour = label)
  ) +
    geom_point(size = 4) +
    scale_colour_manual(values      = grand_colours,
                        labels      = grand_labels_named,
                        name        = "Grand-mean classification") +
    guides(colour = guide_legend(nrow         = 1,
                                 title.hjust  = 0.5,
                                 override.aes = list(size = 4))) +
    coord_cartesian(xlim = c(10, 20)) +
    theme_void() +
    theme(legend.position   = "bottom",
          legend.direction  = "horizontal",
          legend.text       = element_text(size = 10),
          legend.title      = element_text(size = 10.5, face = "bold"),
          legend.key.size   = unit(1.1, "lines"),
          legend.box.margin = margin(0, 0, 0, 0),
          plot.margin       = margin(0, 0, 0, 0))
  
  # ---- Assemble forest plot -------------------------------------------------
  grob_p3     <- ggplotGrob(p3)
  grob_p1     <- ggplotGrob(p1)
  grob_p2     <- ggplotGrob(p2)
  grob_legend <- ggplotGrob(legend_plot)
  
  aligned <- gtable_cbind(grob_p3, grob_p1, grob_p2)
  
  final <- arrangeGrob(
    aligned,
    grob_legend,
    ncol    = 1,
    heights = unit(c(1, 0.06), c("null", "npc"))
  )
  
  f_out_forest <- file.path(path_fig, "combined_rq1_3_sequences_forest.png")
  ggsave(f_out_forest, plot = final, width = 17, height = 18, dpi = 300)
  message("Saved: ", f_out_forest)
  
  # ---- Compact scatter companion -------------------------------------------
  # Same three data.tables, still in seq_order (RQ3 desc). Rank is just
  # position in that ordering; grand means come straight off row 1 of
  # rq1/rq2 (constant across rows, as used in make_panel() above).
  
  rank_tbl <- data.table(sequence = seq_order, rank_loss = seq_along(seq_order))
  
  gm_rq1 <- rq1[1, grand_mean]
  gm_rq2 <- rq2[1, grand_mean]
  
  dt <- data.table(
    sequence  = seq_order,
    rank_loss = rank_tbl[, rank_loss],
    bet_dev   = rq1[match(seq_order, sequence), mu_mean] - gm_rq1,
    stake_dev = rq2[match(seq_order, sequence), mu_mean] - gm_rq2
  )
  
  dt_long <- melt(dt, id.vars = c("sequence", "rank_loss"),
                  measure.vars  = c("bet_dev", "stake_dev"),
                  variable.name = "margin", value.name = "deviation")
  
  dt_long[, margin := factor(margin,
                             levels = c("bet_dev", "stake_dev"),
                             labels = c("Betting-rate Δ",
                                        "Stake-deviation Δ"))]
  
  dt_long[, is_pure := sequence %in% c("HHHHHH", "TTTTTT")]
  
  # Project's standard binary-contrast pair (used for FN/FP, above/below
  # grand mean, and Head/Tail bias throughout the other figure scripts).
  margin_colours <- c("Betting-rate Δ"    = "#2166AC",
                      "Stake-deviation Δ" = "#C0392B")
  
  # Trend lines are OLS, fit only on the 62 non-pure sequences: the two
  # pure sequences sit at the top of the welfare-loss ranking with an
  # extreme, opposite-signed betting deviation (they're the acknowledged
  # exception, not part of the gradient the line is meant to summarise),
  # and including them -- especially in a loess -- pulls the fitted curve
  # into a non-monotonic shape at the boundary that overstates structure
  # the text doesn't claim. A straight line fit on the remainder states
  # the "one margin trends, one doesn't" contrast directly.
  dt_trend <- dt_long[is_pure == FALSE]
  
  # ---- OLS slope table, written to CSV for the figure caption -------------
  trend_stats <- rbindlist(lapply(levels(dt_trend$margin), function(m) {
    fit <- lm(deviation ~ rank_loss, data = dt_trend[margin == m])
    co  <- summary(fit)$coefficients
    data.table(
      margin    = m,
      n         = nobs(fit),
      slope     = co["rank_loss", "Estimate"],
      se        = co["rank_loss", "Std. Error"],
      t_value   = co["rank_loss", "t value"],
      p_value   = co["rank_loss", "Pr(>|t|)"],
      r_squared = summary(fit)$r.squared
    )
  }))
  
  f_out_stats <- file.path(path_out, "combined_rq1_3_sequences_scatter_trend_stats.csv")
  fwrite(trend_stats, f_out_stats)
  message("Saved: ", f_out_stats)
  
  p_scatter <- ggplot(dt_long, aes(x = rank_loss, y = deviation, colour = margin)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
    geom_point(aes(shape = is_pure, size = is_pure), alpha = 0.75) +
    geom_smooth(data = dt_trend, method = "lm", se = TRUE,
                linewidth = 0.9, alpha = 0.12) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 3.2), guide = "none") +
    scale_colour_manual(values = margin_colours, name = NULL) +
    scale_fill_manual(values = margin_colours, guide = "none") +
    scale_x_continuous(name = "Sequence rank (1 = highest welfare loss)") +
    scale_y_continuous(name   = "Deviation from grand mean",
                       labels = percent_format(accuracy = 1)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",
          axis.title      = element_text(size = 9),
          plot.margin     = margin(4, 8, 4, 4))
  
  if (has_ggrepel) {
    p_scatter <- p_scatter +
      ggrepel::geom_text_repel(
        data                = dt_long[is_pure == TRUE],
        aes(label = sequence), colour = "black", size = 2.8,
        family = "mono", fontface = "bold", show.legend = FALSE,
        min.segment.length = 0, seed = 1
      )
  } else {
    p_scatter <- p_scatter +
      geom_text(
        data                = dt_long[is_pure == TRUE],
        aes(label = sequence), colour = "black", size = 2.8,
        family = "mono", fontface = "bold", show.legend = FALSE,
        vjust = -0.8
      )
  }
  
  f_out_scatter <- file.path(path_fig, "combined_rq1_3_sequences_scatter.png")
  ggsave(f_out_scatter, plot = p_scatter, width = 20, height = 10, units = "cm", dpi = 300)
  message("Saved: ", f_out_scatter)
  
  invisible(TRUE)
}