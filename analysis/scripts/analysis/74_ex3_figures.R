# ============================================================
# 74_ex3_figures.R
#
# PURPOSE
#   EX3 combined figures comparing FN (m25) and FP (m19)
#   treatments using the full datasets. Produces only the
#   combined multi-panel figures (individual RQ figures are
#   produced by the standard RQ pipelines).
#
#   (A) ex3_combined_raincloud.png
#       Three panels side by side (RQ1, RQ2, RQ3).
#       Each panel: overlapping transparent distributions for
#       both treatments, grand median lines, jitter strips.
#
#   (B) ex3_combined_forest.png
#       Three panels side by side (RQ1, RQ2, RQ3).
#       Each panel: x-wing forest — sequences on y-axis,
#       two dodged horizontal CIs per row (one per treatment),
#       grand median lines.
#
# INPUT
#   path_out/rq{1,2,3}_{m25,m19}_full_{sequences,participants,model_summary}.csv
#
# OUTPUT
#   path_fig/ex3_combined_raincloud.png
#   path_fig/ex3_combined_forest.png
# ============================================================

ex3_figures <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(ggdist)
  library(scales)
  library(patchwork)
  
  tr_fn <- "m25"
  tr_fp <- "m19"
  tag   <- "full"
  
  tr_colours <- c(m25 = "#2166AC", m19 = "#C0392B")
  tr_labels  <- c(m25 = "FN (m = 2.5)", m19 = "FP (m = 1.9)")
  
  # ---- RQ specs -----------------------------------------------------------
  # All point-estimate columns refer to *_median; CI columns are unchanged.
  rq_specs <- list(
    rq1 = list(
      seq_mu    = "mu_b_median",
      seq_q025  = "mu_b_q025",
      seq_q975  = "mu_b_q975",
      pid_mu    = "mu_i_median",
      gm_param  = "grand mean bet probability",
      x_label_r = expression(hat(mu)[i]^b ~ "(betting probability)"),
      x_label_f = expression(hat(mu)[s]^b ~ "(betting probability)"),
      x_fmt     = percent_format(accuracy = 1)
    ),
    rq2 = list(
      seq_mu    = "mu_a_median",
      seq_q025  = "mu_a_q025",
      seq_q975  = "mu_a_q975",
      pid_mu    = "mu_a_median",
      gm_param  = "grand mean mu_a (absolute scale, proportion of endowment)",
      x_label_r = expression(hat(mu)[i]^a ~ "(stake deviation)"),
      x_label_f = expression(hat(mu)[s]^a ~ "(stake deviation)"),
      x_fmt     = percent_format(accuracy = 1)
    ),
    rq3 = list(
      seq_mu    = "mu_c_median",
      seq_q025  = "mu_c_q025",
      seq_q975  = "mu_c_q975",
      pid_mu    = "mu_c_median",
      gm_param  = "grand mean mu_c (proportion of endowment)",
      x_label_r = expression(hat(mu)[i]^c ~ "(welfare loss)"),
      x_label_f = expression(hat(mu)[s]^c ~ "(welfare loss)"),
      x_fmt     = percent_format(accuracy = 1)
    )
  )
  
  # ---- Helpers ------------------------------------------------------------
  load_rq <- function(rq, tr) {
    stem  <- paste0(rq, "_", tr, "_", tag)
    f_seq <- file.path(path_out, paste0(stem, "_sequences.csv"))
    f_pid <- file.path(path_out, paste0(stem, "_participants.csv"))
    f_mod <- file.path(path_out, paste0(stem, "_model_summary.csv"))
    if (!all(file.exists(c(f_seq, f_pid, f_mod)))) {
      warning("EX3: missing inputs for ", rq, "/", tr)
      return(NULL)
    }
    list(seq = fread(f_seq, encoding = "UTF-8"),
         pid = fread(f_pid, encoding = "UTF-8"),
         mod = fread(f_mod, encoding = "UTF-8"))
  }
  
  # Returns the posterior median plus 95% CI from the model summary.
  get_gm <- function(mod, param) {
    row <- mod[parameter == param]
    if (nrow(row) == 0L) return(list(median = NA_real_, q025 = NA_real_, q975 = NA_real_))
    list(median = row$median[[1]], q025 = row$q025[[1]], q975 = row$q975[[1]])
  }
  
  # =========================================================
  # (A) RAINCLOUD PLOTS
  # =========================================================
  
  # label_pos controls the text annotation attached to each grand-median vline.
  #
  # Structure: a named list with entries $fn and $fp, each containing:
  #   hjust  — horizontal justification of the text relative to its x position.
  #             Negative values push right of the line; values > 1 push left.
  #   vjust  — vertical justification. Combined with y = Inf this determines
  #             how far the label sits below the top edge of the panel.
  #             Larger values move the label further down.
  #   y      — y coordinate for the label. Use Inf for the top of the panel,
  #             -Inf for the bottom, or any numeric value on the ggdist y scale
  #             (densities sit around y = 0; jitter strips around y = -0.15).
  #   x_off  — scalar offset added to the vline x position before placing the
  #             label. Useful when two lines are close and their labels overlap.
  #
  # The defaults below reproduce the original appearance. Edit rain_label_pos
  # (further down) to customise each panel independently.
  default_label_pos <- list(
    fn = list(hjust = -0.1, vjust = 1.5, y = Inf, x_off = 0),
    fp = list(hjust = -0.1, vjust = 3.2, y = Inf, x_off = 0)
  )
  
  make_raincloud <- function(rq, spec, show_legend = FALSE,
                             label_pos = default_label_pos) {
    
    # Merge caller overrides with defaults — partial lists are fine.
    lp <- modifyList(default_label_pos, label_pos)
    
    arts_fn <- load_rq(rq, tr_fn)
    arts_fp <- load_rq(rq, tr_fp)
    if (is.null(arts_fn) || is.null(arts_fp)) return(NULL)
    
    gm_fn <- get_gm(arts_fn$mod, spec$gm_param)
    gm_fp <- get_gm(arts_fp$mod, spec$gm_param)
    
    pid_fn <- copy(arts_fn$pid)
    pid_fp <- copy(arts_fp$pid)
    setnames(pid_fn, spec$pid_mu, "mu_val")
    setnames(pid_fp, spec$pid_mu, "mu_val")
    pid_fn[, treatment := factor(tr_fn, levels = c(tr_fn, tr_fp))]
    pid_fp[, treatment := factor(tr_fp, levels = c(tr_fn, tr_fp))]
    combined <- rbindlist(list(pid_fn[, .(mu_val, treatment)],
                               pid_fp[, .(mu_val, treatment)]))
    
    p <- ggplot(combined, aes(x = mu_val, fill = treatment, colour = treatment)) +
      
      # Grand median bands
      annotate("rect",
               xmin = gm_fn$q025, xmax = gm_fn$q975,
               ymin = -Inf, ymax = Inf,
               fill = tr_colours[tr_fn], alpha = 0.07) +
      annotate("rect",
               xmin = gm_fp$q025, xmax = gm_fp$q975,
               ymin = -Inf, ymax = Inf,
               fill = tr_colours[tr_fp], alpha = 0.07) +
      
      # Grand median lines — both solid
      geom_vline(xintercept = gm_fn$median,
                 colour = tr_colours[tr_fn], linewidth = 0.7) +
      geom_vline(xintercept = gm_fp$median,
                 colour = tr_colours[tr_fp], linewidth = 0.7) +
      
      # Grand median labels — fully driven by label_pos$fn / label_pos$fp
      annotate("text",
               x      = gm_fn$median + lp$fn$x_off,
               y      = lp$fn$y,
               label  = sprintf("%.1f%%", 100 * gm_fn$median),
               hjust  = lp$fn$hjust,
               vjust  = lp$fn$vjust,
               size   = 2.6,
               colour = tr_colours[tr_fn]) +
      annotate("text",
               x      = gm_fp$median + lp$fp$x_off,
               y      = lp$fp$y,
               label  = sprintf("%.1f%%", 100 * gm_fp$median),
               hjust  = lp$fp$hjust,
               vjust  = lp$fp$vjust,
               size   = 2.6,
               colour = tr_colours[tr_fp]) +
      
      # Overlapping halfeye densities — both at y = 0, transparent
      stat_halfeye(
        aes(y = 0),
        adjust       = 0.8,
        width        = 0.6,
        .width       = 0,
        point_colour = NA,
        alpha        = 0.40,
        position     = position_identity()
      ) +
      
      # Jitter strips just below the density baseline
      geom_jitter(
        aes(y = -0.15),
        shape  = 21,
        size   = 1.4,
        alpha  = 0.55,
        stroke = 0.2,
        colour = "white",
        height = 0.04,
        width  = 0
      ) +
      
      scale_fill_manual(values = tr_colours, labels = tr_labels, name = NULL) +
      scale_colour_manual(values = tr_colours, labels = tr_labels, name = NULL) +
      scale_x_continuous(name = spec$x_label_r, labels = spec$x_fmt ) +
      labs(y = NULL) +
      
      theme_classic(base_size = 10) +
      theme(
        axis.title.x       = element_text(size = 9),
        axis.text.y        = element_blank(),
        axis.line.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        legend.position    = if (show_legend) "bottom" else "none",
        legend.text        = element_text(size = 8.5),
        plot.margin        = margin(8, 12, 4, 8)
      )
    
    p
  }
  
  # ---- Per-panel label positions ------------------------------------------
  # Edit any field below to reposition a label in a specific panel.
  # Only the fields you change need to be listed; everything else falls back
  # to default_label_pos defined above.
  #
  # Quick reference:
  #   hjust  < 0   → label starts to the right of the line (most common)
  #   hjust  > 1   → label ends to the left  of the line
  #   vjust  > 0   → moves label down from y  (at y = Inf: away from top edge)
  #   y      = Inf → top of panel; -Inf → bottom; numeric → exact y position
  #                  (density peak ~0.55, density baseline ~0, jitter ~-0.15)
  #   x_off        → nudge label left (negative) or right (positive) of the line
  rain_label_pos <- list(
    rq1 = list(
      fn = list(hjust = -0.2, vjust = 1, y = Inf, x_off = 0),
      fp = list(hjust = 1.2, vjust = 1, y = Inf, x_off = 0)
    ),
    rq2 = list(
      fn = list(hjust = 1.2, vjust = 1, y = Inf, x_off = 0),
      fp = list(hjust = -0.2, vjust = 1, y = Inf, x_off = 0)
    ),
    rq3 = list(
      fn = list(hjust = 1.2, vjust = 1, y = Inf, x_off = 0),
      fp = list(hjust = -0.2, vjust = 1, y = Inf, x_off = 0)
    )
  )
  
  rain_plots <- list(
    rq1 = make_raincloud("rq1", rq_specs$rq1, show_legend = FALSE,
                         label_pos = rain_label_pos$rq1),
    rq2 = make_raincloud("rq2", rq_specs$rq2, show_legend = FALSE,
                         label_pos = rain_label_pos$rq2),
    rq3 = make_raincloud("rq3", rq_specs$rq3, show_legend = FALSE,
                         label_pos = rain_label_pos$rq3)
  )
  
  p1 <- rain_plots$rq1
  p2 <- rain_plots$rq2
  p3 <- rain_plots$rq3
  
  valid_rain <- Filter(Negate(is.null), rain_plots)
  if (length(valid_rain) == 3) {
    
    combined_rain <- (p1 + theme(plot.margin = margin(8, 8,  8, 0))) +
      (p2 + theme(plot.margin = margin(8, 0, 8, 8)))  +
      (p3 + theme(plot.margin = margin(8, 0, 8, 0))) +
      plot_layout(ncol = 2, nrow = 2, heights = c(1, 1)) +
      plot_layout(design = "AB
                        CC")
    
    f_out <- file.path(path_fig, "ex3_combined_raincloud.png")
    ggsave(f_out, combined_rain, width = 8, height = 5, dpi = 300)
    msg("Saved: ", f_out)
  }
  
  # =========================================================
  # (B) FOREST PLOTS
  # =========================================================
  
  # Shared sequence order: descending by RQ3/FN posterior median welfare loss.
  arts_rq3_fn      <- load_rq("rq3", tr_fn)
  seq_order_shared  <- arts_rq3_fn$seq[order(-mu_c_median), sequence]
  seq_levels_shared <- rev(seq_order_shared)
  
  make_forest <- function(rq, spec, seq_levels, show_y_labels = FALSE, show_legend = FALSE) {
    
    arts_fn <- load_rq(rq, tr_fn)
    arts_fp <- load_rq(rq, tr_fp)
    if (is.null(arts_fn) || is.null(arts_fp)) return(NULL)
    
    gm_fn <- get_gm(arts_fn$mod, spec$gm_param)
    gm_fp <- get_gm(arts_fp$mod, spec$gm_param)
    
    seq_fn <- copy(arts_fn$seq)[, treatment := tr_fn]
    seq_fp <- copy(arts_fp$seq)[, treatment := tr_fp]
    
    setnames(seq_fn, spec$seq_mu,   "mu_median")
    setnames(seq_fn, spec$seq_q025, "mu_q025")
    setnames(seq_fn, spec$seq_q975, "mu_q975")
    setnames(seq_fp, spec$seq_mu,   "mu_median")
    setnames(seq_fp, spec$seq_q025, "mu_q025")
    setnames(seq_fp, spec$seq_q975, "mu_q975")
    
    seq_fn[, sequence := factor(sequence, levels = seq_levels)]
    seq_fp[, sequence := factor(sequence, levels = seq_levels)]
    
    combined <- rbindlist(list(
      seq_fn[, .(sequence, mu_median, mu_q025, mu_q975, treatment)],
      seq_fp[, .(sequence, mu_median, mu_q025, mu_q975, treatment)]
    ))
    combined[, treatment := factor(treatment, levels = c(tr_fn, tr_fp))]
    
    n_seq <- length(seq_levels)
    dodge <- position_dodge(width = 0.6)
    
    p <- ggplot(combined,
                aes(x = mu_median, y = sequence,
                    colour = treatment, group = treatment)) +
      
      # Grand median bands
      annotate("rect",
               xmin = gm_fn$q025, xmax = gm_fn$q975,
               ymin = -Inf, ymax = Inf,
               fill = tr_colours[tr_fn], alpha = 0.07) +
      annotate("rect",
               xmin = gm_fp$q025, xmax = gm_fp$q975,
               ymin = -Inf, ymax = Inf,
               fill = tr_colours[tr_fp], alpha = 0.07) +
      
      # Grand median lines — both solid
      geom_vline(xintercept = gm_fn$median,
                 colour = tr_colours[tr_fn], linewidth = 0.55) +
      geom_vline(xintercept = gm_fp$median,
                 colour = tr_colours[tr_fp], linewidth = 0.55) +
      
      # Grand median labels at top
      annotate("text",
               x = gm_fn$median, y = n_seq + 0.6,
               label = sprintf("FN: %.1f%%", 100 * gm_fn$median),
               hjust = -0.1, vjust = 0.5, size = 2.4,
               colour = tr_colours[tr_fn]) +
      annotate("text",
               x = gm_fp$median, y = n_seq + 0.6,
               label = sprintf("FP: %.1f%%", 100 * gm_fp$median),
               hjust = -0.1, vjust = 1.8, size = 2.4,
               colour = tr_colours[tr_fp]) +
      
      # Dodged CI segments
      geom_errorbarh(aes(xmin = mu_q025, xmax = mu_q975),
                     height    = 0,
                     linewidth = 0.4,
                     alpha     = 0.65,
                     position  = dodge) +
      
      # Dodged point estimates (posterior medians)
      geom_point(size = 1.6, position = dodge) +
      
      scale_colour_manual(values = tr_colours, labels = tr_labels,
                          name = "Treatment") +
      scale_x_continuous(name = spec$x_label_f, labels = spec$x_fmt) +
      scale_y_discrete(name = NULL) +
      
      theme_classic(base_size = 10) +
      theme(
        axis.title.x       = element_text(size = 8.5),
        axis.title.y       = element_blank(),
        axis.text.y        = if (show_y_labels)
          element_text(size = 6.5, family = "mono")
        else
          element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
        legend.position    = if (show_legend) "bottom" else "none",
        legend.text        = element_text(size = 9),
        legend.title       = element_text(size = 9, face = "bold"),
        plot.margin        = margin(4, 12, 4, 8)
      )
    
    p
  }
  
  forest_plots <- list(
    rq1 = make_forest("rq1", rq_specs$rq1, seq_levels_shared, show_y_labels = FALSE, show_legend = FALSE),
    rq2 = make_forest("rq2", rq_specs$rq2, seq_levels_shared, show_y_labels = FALSE, show_legend = TRUE),
    rq3 = make_forest("rq3", rq_specs$rq3, seq_levels_shared, show_y_labels = TRUE,  show_legend = FALSE)
  )
  
  valid_forest <- Filter(Negate(is.null), forest_plots)
  if (length(valid_forest) == 3) {
    
    combined_forest <- forest_plots$rq3 + forest_plots$rq1 + forest_plots$rq2 +
      plot_layout(ncol = 3, guides = "collect") &
      theme(legend.position = "bottom")
    
    f_out <- file.path(path_fig, "ex3_combined_forest.png")
    ggsave(f_out, combined_forest, width = 18, height = 14, dpi = 300)
    msg("Saved: ", f_out)
  }
  
  invisible(TRUE)
}