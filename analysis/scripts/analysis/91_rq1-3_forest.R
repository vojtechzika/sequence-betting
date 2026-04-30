# ============================================================
# 91_rq1_3_forest.R
#
# Panel order:  RQ3 (with y-axis labels) | RQ1 | RQ2
# Sequences:    ordered by RQ3 posterior mean welfare loss (high→low)
# Label colours: 0 RQs flagged (RQ1+RQ2) = grey, 1 = #1B7837, 2 = #CC0000
# Legend:       standalone ggplot row below the panels
# No title or caption.
# ============================================================

rq1_3_forest <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(scales)
  
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
                   x_label         = expression(hat(mu)[s]^c ~ "(welfare loss, share of endowment)"),
                   show_y_labels   = TRUE,
                   y_label_colours = label_colours)
  p1 <- make_panel(rq1, expression(hat(mu)[s]^b ~ "(betting probability)"))
  p2 <- make_panel(rq2, expression(hat(mu)[s]^a ~ "(stake deviation, share of endowment)"))
  
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
  
  # ---- Assemble -----------------------------------------------------------
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
  
  f_out <- file.path(path_fig, "combined_rq1_3_sequences_forest.png")
  ggsave(f_out, plot = final, width = 17, height = 18, dpi = 300)
  message("Saved: ", f_out)
  
  invisible(TRUE)
}