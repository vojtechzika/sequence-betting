# ============================================================
# 44_rq4_figures.R
#   RQ4 figures: sequence-level forest plot and participant-level
#   raincloud plot for side choice (Heads vs Tails)
#
# PURPOSE
#   Mirrors 14_rq1_figures.R in structure and style.
#   Key difference from RQ1: the reference line is hbar (the
#   population baseline from the model), not a normative
#   benchmark. The colour scale is two-sided: Head-biased
#   sequences in red, Tail-biased in blue, neutral in grey.
#   Also exports forest data for the combined multi-RQ plot
#   produced by 91_rq1-3_forest.R (extended to include RQ4).
#
# INPUT
#   path_out/rq4_<tr>_sequences.csv
#   path_out/rq4_<tr>_participants.csv
#   path_out/rq4_<tr>_model_summary.csv
#
# OUTPUT
#   path_fig/rq4_<tr>_sequences_forest.png
#   path_fig/rq4_<tr>_participants_raincloud.png
#   path_out/rq4_<tr>_forest_data.csv
#
# CONFIG KEYS USED
#   cfg$run$treatment
#   cfg$run$overwrite_outputs
#
# CALL ORDER
#   rq4_tables(cfg)  -- 43_
#   rq4_figures(cfg) -- this script
# ============================================================

library(data.table)
library(ggplot2)
library(ggdist)
library(scales)

rq4_figures <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L)
  
  dir.create(path_fig, showWarnings = FALSE, recursive = TRUE)
  
  # Two-sided colour palette: Head-biased -> warm reds, Tail-biased -> cool blues
  direction_colors <- c(
    strong_head   = "#C0392B",
    moderate_head = "#E8A090",
    weak_head     = "#F5C9C2",
    neutral       = "#AAAAAA",
    weak_tail     = "#C2D8F0",
    moderate_tail = "#74ADD1",
    strong_tail   = "#2166AC"
  )
  
  direction_labels <- c(
    strong_head   = "Strong Head bias",
    moderate_head = "Moderate Head bias",
    weak_head     = "Weak Head bias",
    neutral       = "Neutral",
    weak_tail     = "Weak Tail bias",
    moderate_tail = "Moderate Tail bias",
    strong_tail   = "Strong Tail bias"
  )
  
  direction_levels <- names(direction_colors)
  
  for (tr in tr_vec) {
    
    f_seq_csv <- file.path(path_out, paste0("rq4_", tr, "_sequences.csv"))
    f_pid_csv <- file.path(path_out, paste0("rq4_", tr, "_participants.csv"))
    f_mod_csv <- file.path(path_out, paste0("rq4_", tr, "_model_summary.csv"))
    
    if (!file.exists(f_seq_csv) || !file.exists(f_pid_csv) || !file.exists(f_mod_csv)) {
      warning("RQ4 figures: missing table inputs for tr='", tr, "'. Skipping.")
      next
    }
    
    seq_tbl <- fread(f_seq_csv, encoding = "UTF-8")
    pid_tbl <- fread(f_pid_csv, encoding = "UTF-8")
    mod_tbl <- fread(f_mod_csv, encoding = "UTF-8")
    
    # Reference line: hbar (population baseline from model summary)
    hbar_row  <- mod_tbl[parameter == "hbar (population baseline Heads prob)"]
    hbar_val  <- hbar_row$median
    hbar_lo   <- hbar_row$q025
    hbar_hi   <- hbar_row$q975
    
    stopifnot(length(hbar_val) == 1L, is.finite(hbar_val))
    
    seq_tbl[, direction_label := factor(direction_label, levels = direction_levels)]
    pid_tbl[, side_label      := factor(side_label,      levels = direction_levels)]
    
    # ----------------------------------------
    # Forest plot data export (for 91_ combined plot)
    # ----------------------------------------
    f_forest_data <- file.path(path_out, paste0("rq4_", tr, "_forest_data.csv"))
    if (!should_skip(f_forest_data, cfg, "output",
                     paste0("RQ4 forest data (", tr, ")"))) {
      
      forest_data <- seq_tbl[, .(
        sequence        = sequence,
        mu_mean         = mu_h_mean,
        mu_median       = mu_h_median,
        mu_q025         = mu_h_q025,
        mu_q975         = mu_h_q975,
        direction_label = direction_label,
        hbar            = hbar_val,
        hbar_lo         = hbar_lo,
        hbar_hi         = hbar_hi,
        treatment       = tr,
        rq              = "rq4"
      )]
      fwrite(forest_data, f_forest_data)
      msg("Saved: ", f_forest_data)
    }
    
    # ----------------------------------------
    # Forest plot — three panels
    # ----------------------------------------
    f_forest <- file.path(path_fig, paste0("rq4_", tr, "_sequences_forest.png"))
    if (!should_skip(f_forest, cfg, "output",
                     paste0("RQ4 forest plot (", tr, ")"))) {
      
      dt <- copy(seq_tbl)
      dt[, sequence := factor(sequence,
                              levels = dt[order(mu_h_mean), sequence])]
      
      # Assign panel
      dt[, panel := fcase(
        direction_label %in% c("strong_head", "moderate_head", "weak_head"), "Head-biased",
        direction_label == "neutral",                                          "Neutral",
        direction_label %in% c("weak_tail", "moderate_tail", "strong_tail"),  "Tail-biased"
      )]
      dt[, panel := factor(panel, levels = c("Head-biased", "Neutral", "Tail-biased"))]
      
      p_forest <- ggplot(dt, aes(y = sequence)) +
        
        # Posterior uncertainty band around hbar
        annotate("rect",
                 xmin = hbar_lo, xmax = hbar_hi,
                 ymin = -Inf,    ymax = Inf,
                 fill = "steelblue", alpha = 0.08) +
        
        # hbar reference line
        geom_vline(xintercept = hbar_val,
                   colour    = "steelblue",
                   linewidth = 0.65,
                   linetype  = "solid") +
        
        # 0.5 chance line
        geom_vline(xintercept = 0.5,
                   colour    = "grey40",
                   linewidth = 0.45,
                   linetype  = "dashed") +
        
        # CIs
        geom_segment(aes(x     = mu_h_q025,
                         xend  = mu_h_q975,
                         y     = sequence,
                         yend  = sequence,
                         colour = direction_label),
                     linewidth = 0.65,
                     alpha     = 0.90) +
        
        # Point estimates
        geom_point(aes(x = mu_h_mean, colour = direction_label),
                   size = 1.8) +
        
        scale_colour_manual(
          values = direction_colors,
          labels = direction_labels,
          name   = NULL,
          drop   = FALSE
        ) +
        
        scale_x_continuous(
          name   = expression(hat(mu)[s]^h ~
                                "(posterior mean Heads probability)"),
          limits = c(
            floor(min(dt$mu_h_q025) * 10) / 10 - 0.03,
            ceiling(max(dt$mu_h_q975) * 10) / 10 + 0.03
          ),
          breaks = seq(0.2, 0.8, by = 0.1),
          labels = percent_format(accuracy = 1)
        ) +
        
        scale_y_discrete(name = NULL) +
        
        annotate("text",
                 x     = hbar_val,
                 y     = Inf,
                 label = sprintf("bar(h) == %.2f", hbar_val),
                 parse = TRUE,
                 hjust = -0.10,
                 vjust = 1.5,
                 size  = 2.8,
                 colour = "steelblue") +
        
        annotate("text",
                 x     = 0.5,
                 y     = Inf,
                 label = "Fair coin (0.50)",
                 hjust =  1.08,
                 vjust =  1.5,
                 size  = 2.8,
                 colour = "grey40") +
        
        facet_wrap(
          ~ panel,
          ncol   = 3,
          scales = "free_y"
        ) +
        
        theme_classic(base_size = 10) +
        theme(
          axis.text.y        = element_text(size = 7, family = "mono"),
          axis.title.x       = element_text(size = 9),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          legend.position    = "bottom",
          legend.text        = element_text(size = 8),
          strip.text         = element_text(size = 9, face = "bold"),
          strip.background   = element_blank(),
          plot.margin        = margin(10, 20, 8, 8)
        ) +
        
        guides(colour = guide_legend(
          nrow         = 2,
          byrow        = TRUE,
          override.aes = list(size = 2.5, linewidth = 0.8)
        ))
      
      # Height scales with the largest panel; width covers three columns
      n_panels  <- dt[, .(n = uniqueN(sequence)), by = panel]
      max_n     <- max(n_panels$n)
      fig_h     <- max(4, max_n * 0.28 + 2)   # ~0.28 in per sequence row + margins
      
      ggsave(f_forest, p_forest, width = 14, height = fig_h, dpi = 300)
      msg("Saved: ", f_forest)
    }
    
    # ----------------------------------------
    # Participant raincloud
    # ----------------------------------------
    f_rain <- file.path(path_fig, paste0("rq4_", tr, "_participants_raincloud.png"))
    if (!should_skip(f_rain, cfg, "output",
                     paste0("RQ4 raincloud (", tr, ")"))) {
      
      dt_pid <- copy(pid_tbl)
      
      # Raincloud layout:
      #   y = 0.20  jitter strip
      #   y = 0.35  density baseline (slab grows upward)
      # x-axis shows the 2nd-98th percentile range with 5pp padding;
      # outlier points beyond those limits are squished to the boundary.
      
      x_lo <- max(0,   quantile(dt_pid$mu_h_mean, 0.02)  - 0.05)
      x_hi <- min(1,   quantile(dt_pid$mu_h_mean, 0.98)  + 0.05)
      
      p_rain <- ggplot(dt_pid, aes(x = mu_h_mean)) +
        
        # hbar posterior band
        annotate("rect",
                 xmin = hbar_lo, xmax = hbar_hi,
                 ymin = -Inf,    ymax = Inf,
                 fill = "steelblue", alpha = 0.10) +
        
        # hbar line
        geom_vline(xintercept = hbar_val,
                   colour    = "steelblue",
                   linewidth = 0.70,
                   linetype  = "solid") +
        
        # 0.5 reference line: longdash clearly distinct from solid hbar
        geom_vline(xintercept = 0.5,
                   colour    = "grey55",
                   linewidth = 0.45,
                   linetype  = "longdash") +
        
        # Density slab: height relative to data range, trimmed to data support
        stat_slab(
          aes(y = 0.35),
          fill      = "#74ADD1",
          colour    = NA,
          adjust    = 2.0,
          height    = 0.40,
          alpha     = 0.82,
          normalize = "panels",
          trim      = FALSE
        ) +
        
        # Jittered points at y = 0.20, below the density baseline
        geom_jitter(
          aes(y = 0.20),
          fill   = "#2166AC",
          shape  = 21,
          size   = 1.7,
          alpha  = 0.72,
          stroke = 0.25,
          colour = "white",
          height = 0.055,
          width  = 0
        ) +
        
        # hbar label top of panel
        annotate("text",
                 x     = hbar_val,
                 y     = Inf,
                 label = sprintf("bar(h) == %.2f", hbar_val),
                 parse = TRUE,
                 hjust = -0.08,
                 vjust = 1.5,
                 size  = 2.9,
                 colour = "steelblue") +
        
        scale_x_continuous(
          name   = expression(hat(mu)[i]^h ~
                                "(posterior mean Heads probability)"),
          limits = c(x_lo, x_hi),
          oob    = scales::squish,
          breaks = seq(0, 1, by = 0.1),
          labels = percent_format(accuracy = 1)
        ) +
        
        scale_y_continuous(name = NULL, breaks = NULL,
                           expand = expansion(mult = c(0.02, 0.12))) +
        
        theme_classic(base_size = 11) +
        theme(
          axis.title.x       = element_text(size = 10),
          axis.text.y        = element_blank(),
          axis.line.y        = element_blank(),
          axis.ticks.y       = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          legend.position    = "none",
          plot.margin        = margin(14, 16, 6, 8)
        )
      
      ggsave(f_rain, p_rain, width = 8, height = 3.2, dpi = 300)
      msg("Saved: ", f_rain)
    }
  }
  
  invisible(TRUE)
}