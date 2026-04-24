# ============================================================
# 14_rq1_figures.R
#
# PURPOSE
#   Produces paper-ready figures for RQ1 (extensive margin of
#   betting): sequence-level forest plot and participant-level
#   raincloud plot. Also exports forest plot data for combined
#   multi-RQ sequence forest plot.
#
# INPUT
#   path_out/rq1_<tr>_<tag>[_bb]_sequences.csv
#   path_out/rq1_<tr>_<tag>[_bb]_participants.csv
#   path_out/rq1_<tr>_<tag>[_bb]_model_summary.csv
#
# OUTPUT
#   path_fig/rq1_<tr>_<tag>[_bb]_sequences_forest.png
#   path_fig/rq1_<tr>_<tag>[_bb]_participants_raincloud.png
#   path_out/rq1_<tr>_<tag>[_bb]_forest_data.csv
#
# TAGS
#   full          -- all participants
#   confirmatory  -- normative betters only
#
# ROBUSTNESS
#   Called with robustness = TRUE to produce _bb outputs.
# ============================================================

rq1_figures <- function(cfg, robustness = FALSE) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  suffix <- if (robustness) "_bb" else ""
  tags   <- c("full", "confirmatory")
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      file_stem <- paste0("rq1_", tr, "_", tag, suffix)
      
      f_seq_csv <- file.path(path_out, paste0(file_stem, "_sequences.csv"))
      f_pid_csv <- file.path(path_out, paste0(file_stem, "_participants.csv"))
      f_mod_csv <- file.path(path_out, paste0(file_stem, "_model_summary.csv"))
      
      if (!file.exists(f_seq_csv) || !file.exists(f_pid_csv) || !file.exists(f_mod_csv)) {
        warning("RQ1 figures: missing table inputs for tr='", tr, "', tag='", tag, "'. Skipping.")
        next
      }
      
      seq_tbl <- fread(f_seq_csv, encoding = "UTF-8")
      pid_tbl <- fread(f_pid_csv, encoding = "UTF-8")
      mod_tbl <- fread(f_mod_csv, encoding = "UTF-8")
      
      gm         <- mod_tbl[parameter == "grand mean bet probability"]
      grand_mean <- gm$median
      grand_lo   <- gm$q025
      grand_hi   <- gm$q975
      
      # ----------------------------------------
      # Forest plot data export
      # ----------------------------------------
      f_forest_data <- file.path(path_out, paste0(file_stem, "_forest_data.csv"))
      if (!should_skip(f_forest_data, cfg, "output",
                       paste0("RQ1 forest data (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        forest_data <- seq_tbl[, .(
          sequence    = sequence,
          mu_mean     = mu_b_mean,
          mu_median   = mu_b_median,
          mu_q025     = mu_b_q025,
          mu_q975     = mu_b_q975,
          underbet_label = underbet_label,
          grand_label    = grand_label,
          grand_mean     = grand_mean,
          grand_lo       = grand_lo,
          grand_hi       = grand_hi,
          treatment      = tr,
          tag            = tag,
          rq             = "rq1"
        )]
        fwrite(forest_data, f_forest_data)
        msg("Saved: ", f_forest_data)
      }
      
      # ----------------------------------------
      # Forest plot
      # ----------------------------------------
      f_forest <- file.path(path_fig, paste0(file_stem, "_sequences_forest.png"))
      if (!should_skip(f_forest, cfg, "output",
                       paste0("RQ1 forest plot (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        dt <- copy(seq_tbl)
        dt[, sequence    := factor(sequence, levels = dt[order(mu_b_mean), sequence])]
        dt[, grand_label := factor(grand_label, levels = c("above", "neutral", "below"))]
        
        label_colors <- c(above = "#2166AC", neutral = "#AAAAAA", below = "#C0392B")
        label_names  <- c(above = "Above grand mean", neutral = "Neutral",
                          below = "Below grand mean")
        
        p <- ggplot(dt, aes(y = sequence)) +
          
          annotate("rect",
                   xmin = grand_lo, xmax = grand_hi,
                   ymin = -Inf, ymax = Inf,
                   fill = "steelblue", alpha = 0.10) +
          
          geom_vline(xintercept = grand_mean,
                     colour = "steelblue", linewidth = 0.6, linetype = "solid") +
          
          geom_vline(xintercept = 1.0,
                     colour = "grey30", linewidth = 0.5, linetype = "dashed") +
          
          geom_segment(aes(x = mu_b_q025, xend = mu_b_q975,
                           y = sequence,  yend = sequence,
                           colour = grand_label),
                       linewidth = 0.5, alpha = 0.7) +
          
          geom_point(aes(x = mu_b_mean, colour = grand_label), size = 1.8) +
          
          scale_colour_manual(values = label_colors, name = NULL, labels = label_names) +
          
          scale_x_continuous(
            name   = expression(hat(mu)[s]^b ~ "(posterior mean betting probability)"),
            limits = c(0.55, 1.02),
            breaks = seq(0.6, 1.0, by = 0.1),
            labels = scales::percent_format(accuracy = 1)
          ) +
          
          scale_y_discrete(name = NULL) +
          
          annotate("text", x = grand_mean, y = 0.4,
                   label = sprintf("Grand mean = %.2f", grand_mean),
                   hjust = -0.05, size = 2.8, colour = "steelblue") +
          
          annotate("text", x = 1.0, y = 0.4,
                   label = "Normative\nbenchmark",
                   hjust = 1.05, size = 2.8, colour = "grey30") +
          
          theme_classic(base_size = 10) +
          theme(
            axis.text.y        = element_text(size = 7, family = "mono"),
            axis.title.x       = element_text(size = 9),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position    = "bottom",
            legend.text        = element_text(size = 8),
            plot.margin        = margin(8, 12, 8, 8)
          )
        
        ggsave(f_forest, p, width = 7, height = 14, dpi = 300)
        msg("Saved: ", f_forest)
      }
      
      # ----------------------------------------
      # Participant raincloud
      # ----------------------------------------
      f_rain <- file.path(path_fig, paste0(file_stem, "_participants_raincloud.png"))
      if (!should_skip(f_rain, cfg, "output",
                       paste0("RQ1 raincloud (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        dt_pid <- copy(pid_tbl)
        
        p <- ggplot(dt_pid, aes(x = mu_i_mean, y = 1)) +
          
          geom_vline(xintercept = grand_mean,
                     colour = "steelblue", linewidth = 0.6, linetype = "solid") +
          
          annotate("text", x = grand_mean, y = 1.55,
                   label = sprintf("Grand mean = %.2f", grand_mean),
                   hjust = -0.05, size = 2.8, colour = "steelblue") +
          
          stat_halfeye(
            fill         = "#74ADD1",
            adjust       = 0.8,
            width        = 0.5,
            .width       = 0,
            point_colour = NA,
            alpha        = 0.80,
            position     = position_nudge(y = 0.12)
          ) +
          
          geom_jitter(
            fill   = "#2166AC",
            shape  = 21,
            size   = 1.8,
            alpha  = 0.75,
            stroke = 0.25,
            colour = "white",
            height = 0.07,
            width  = 0
          ) +
          
          scale_x_continuous(
            name   = expression(hat(mu)[i]^b ~ "(posterior mean betting probability)"),
            limits = c(-0.05, 1.05),
            breaks = seq(0, 1, by = 0.1),
            labels = scales::percent_format(accuracy = 1)
          ) +
          
          scale_y_continuous(name = NULL, breaks = NULL) +
          
          theme_classic(base_size = 11) +
          theme(
            axis.title.x       = element_text(size = 10),
            axis.text.y        = element_blank(),
            axis.line.y        = element_blank(),
            axis.ticks.y       = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position    = "none",
            plot.margin        = margin(8, 12, 8, 8)
          )
        
        ggsave(f_rain, p, width = 8, height = 3, dpi = 300)
        msg("Saved: ", f_rain)
      }
    }
  }
  
  invisible(TRUE)
}