# ============================================================
# 34_rq3_figures.R
#
# PURPOSE
#   Produces paper-ready figures for RQ3 (certainty-equivalent
#   welfare loss): sequence-level forest plot, participant-level
#   raincloud plot, and forest data export for combined plots.
#   Reads selected model from rq3_diagnostics.csv.
#
# INPUT
#   path_out/rq3_diagnostics.csv
#   path_mod/rq3_fit_sequences_<tr>_<tag>[_gamma|_alt].rds
#   path_mod/rq3_pid_levels_<tr>_<tag>[_gamma|_alt].rds
#   path_mod/rq3_seq_levels_<tr>_<tag>[_gamma|_alt].rds
#   path_out/rq3_<tr>_<tag>_participants.csv
#   path_out/rq3_<tr>_<tag>_model_summary.csv
#
# OUTPUT
#   path_fig/rq3_<tr>_<tag>_sequences_forest.png
#   path_fig/rq3_<tr>_<tag>_participants_raincloud.png
#   path_out/rq3_<tr>_<tag>_forest_data.csv
#
# TAGS
#   full          -- all participants
#   confirmatory  -- normative betters only
# ============================================================

rq3_figures <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  tags   <- c("full", "confirmatory")
  
  # ---- Model suffix map ----
  suffix_map <- list(
    primary     = "",
    gamma_only  = "_gamma",
    alternative = "_alt"
  )
  
  # ---- Read diagnostics to determine selected model ----
  f_diag <- file.path(path_out, "rq3_diagnostics.csv")
  diag   <- if (file.exists(f_diag)) fread(f_diag) else NULL
  
  get_selected_suffix <- function(tr, tg) {
    if (is.null(diag)) return("")
    row <- diag[treatment == tr & tag == tg]
    if (nrow(row) == 0L) return("")
    sfx <- suffix_map[[row$selected_model]]
    if (is.null(sfx)) "" else sfx
  }
  
  rho_vec     <- as.numeric(cfg$design$rq3$rho)
  rho_main    <- rho_vec[1]
  rho_nm_main <- gsub("\\.", "", sprintf("%.2f", rho_main))
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      suffix    <- get_selected_suffix(tr, tag)
      file_stem <- paste0("rq3_", tr, "_", tag)
      
      f_fit     <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tag, suffix, ".rds"))
      f_pid     <- file.path(path_mod, paste0("rq3_pid_levels_",    tr, "_", tag, suffix, ".rds"))
      f_seq     <- file.path(path_mod, paste0("rq3_seq_levels_",    tr, "_", tag, suffix, ".rds"))
      f_pid_csv <- file.path(path_out, paste0(file_stem, "_participants.csv"))
      f_mod_csv <- file.path(path_out, paste0(file_stem, "_model_summary.csv"))
      
      if (!all(file.exists(c(f_fit, f_pid, f_seq, f_pid_csv, f_mod_csv)))) {
        warning("RQ3 figures: missing inputs for tr='", tr, "', tag='", tag, "'. Skipping.")
        next
      }
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      pid_tbl    <- fread(f_pid_csv, encoding = "UTF-8")
      mod_tbl    <- fread(f_mod_csv, encoding = "UTF-8")
      
      # ---- Extract posterior draws ----
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      stopifnot(!is.null(post$mu_c), !is.null(post$mu_c_i))
      
      mu_s_draws <- post$mu_c
      mu_i_draws <- post$mu_c_i
      
      # ---- Grand mean ----
      grand_draws    <- apply(mu_s_draws, 1, mean)
      grand_mean_val <- median(grand_draws)
      grand_lo       <- quantile(grand_draws, 0.025)
      grand_hi       <- quantile(grand_draws, 0.975)
      
      # ---- Sequence table ----
      infile <- file.path(path_src, "master_sequences.csv")
      master <- fread(infile, encoding = "UTF-8")
      master[, pid   := as.character(pid)]
      master[, treat := as.character(treat)]
      master[, seq   := as.character(seq)]
      
      d <- master[treat == tr & pid %in% pid_levels]
      n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
      setkey(n_trials_by_seq, seq)
      
      seq_tbl <- data.table(
        sequence    = seq_levels,
        n_trials    = as.integer(n_trials_by_seq[.(seq_levels), n_trials]),
        mu_c_median = apply(mu_s_draws, 2, median),
        mu_c_mean   = apply(mu_s_draws, 2, mean),
        mu_c_q025   = apply(mu_s_draws, 2, quantile, probs = 0.025),
        mu_c_q975   = apply(mu_s_draws, 2, quantile, probs = 0.975)
      )
      seq_tbl[is.na(n_trials), n_trials := 0L]
      
      for (rho in rho_vec) {
        nm <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        seq_tbl[, (nm) := apply(mu_s_draws, 2, function(x) mean(x > rho))]
      }
      
      colL <- paste0("L_rho_", rho_nm_main)
      seq_tbl[, loss_label :=
                fifelse(get(colL) >= 0.95, "strong",
                        fifelse(get(colL) >= 0.80, "moderate",
                                fifelse(get(colL) >= 0.50, "weak", "neutral")))]
      
      seq_tbl[, p_above_grand := apply(mu_s_draws, 2,
                                       function(x) mean(x > grand_draws))]
      seq_tbl[, p_below_grand := apply(mu_s_draws, 2,
                                       function(x) mean(x < grand_draws))]
      seq_tbl[, grand_label := fifelse(
        p_above_grand >= 0.95, "above",
        fifelse(p_above_grand >= 0.80, "likely_above",
                fifelse(p_below_grand >= 0.95, "below",
                        fifelse(p_below_grand >= 0.80, "likely_below",
                                "neutral")))
      )]
      
      setorder(seq_tbl, sequence)
      
      # ---- Forest data export ----
      f_forest_data <- file.path(path_out, paste0(file_stem, "_forest_data.csv"))
      if (!should_skip(f_forest_data, cfg, "output",
                       paste0("RQ3 forest data (", tr, "/", tag, ")"))) {
        forest_data <- seq_tbl[, .(
          sequence    = sequence,
          mu_mean     = mu_c_mean,
          mu_median   = mu_c_median,
          mu_q025     = mu_c_q025,
          mu_q975     = mu_c_q975,
          loss_label  = loss_label,
          grand_label = grand_label,
          grand_mean  = grand_mean_val,
          grand_lo    = as.numeric(grand_lo),
          grand_hi    = as.numeric(grand_hi),
          treatment   = tr,
          tag         = tag,
          rq          = "rq3"
        )]
        fwrite(forest_data, f_forest_data)
        msg("Saved: ", f_forest_data)
      }
      
      # ---- Forest plot ----
      f_forest <- file.path(path_fig, paste0(file_stem, "_sequences_forest.png"))
      if (!should_skip(f_forest, cfg, "output",
                       paste0("RQ3 forest plot (", tr, "/", tag, ")"))) {
        
        dt <- copy(seq_tbl)
        dt[, sequence    := factor(sequence, levels = dt[order(mu_c_mean), sequence])]
        
        dt[, grand_label := factor(grand_label,
                                   levels = c("above", "likely_above", "neutral",
                                              "likely_below", "below"))]
        
        label_colors <- c(above        = "#C0392B",
                          likely_above = "#E8A090",
                          neutral      = "#AAAAAA",
                          likely_below = "#90B8D8",
                          below        = "#2166AC")
        label_names  <- c(above        = "Above grand mean",
                          likely_above = "Likely above grand mean",
                          neutral      = "Neutral",
                          likely_below = "Likely below grand mean",
                          below        = "Below grand mean")
        
        p <- ggplot(dt, aes(y = sequence)) +
          
          annotate("rect",
                   xmin = as.numeric(grand_lo), xmax = as.numeric(grand_hi),
                   ymin = -Inf, ymax = Inf,
                   fill = "steelblue", alpha = 0.10) +
          
          geom_vline(xintercept = grand_mean_val,
                     colour = "steelblue", linewidth = 0.6, linetype = "solid") +
          
          geom_vline(xintercept = 0,
                     colour = "grey30", linewidth = 0.5, linetype = "dashed") +
          
          geom_segment(aes(x = mu_c_q025, xend = mu_c_q975,
                           y = sequence,  yend = sequence,
                           colour = grand_label),
                       linewidth = 0.5, alpha = 0.7) +
          
          geom_point(aes(x = mu_c_mean, colour = grand_label), size = 1.8) +
          
          scale_colour_manual(values = label_colors, name = NULL,
                              labels = label_names) +
          
          scale_x_continuous(
            name   = expression(hat(mu)[s]^c ~ "(posterior mean welfare loss, share of endowment)"),
            labels = scales::percent_format(accuracy = 1)
          ) +
          
          scale_y_discrete(name = NULL) +
          
          annotate("text", x = grand_mean_val, y = Inf,
                   label = sprintf("Grand mean = %.1f%%", 100 * grand_mean_val),
                   hjust = -0.05, vjust = 1.5,
                   size = 2.8, colour = "steelblue") +
          
          annotate("text", x = 0, y = Inf,
                   label = "Zero loss",
                   hjust = -0.05, vjust = 1.5,
                   size = 2.8, colour = "grey30") +
          
          theme_classic(base_size = 10) +
          theme(
            axis.text.y        = element_text(size = 7, family = "mono"),
            axis.title.x       = element_text(size = 9),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position    = "bottom",
            legend.text        = element_text(size = 8),
            plot.margin = margin(8, 12, 8, 20)  # extra left margin
          )
        
        ggsave(f_forest, p, width = 7, height = 14, dpi = 300)
        msg("Saved: ", f_forest)
      }
      
      # ---- Participant raincloud ----
      f_rain <- file.path(path_fig, paste0(file_stem, "_participants_raincloud.png"))
      if (!should_skip(f_rain, cfg, "output",
                       paste0("RQ3 raincloud (", tr, "/", tag, ")"))) {
        
        p <- ggplot(pid_tbl, aes(x = mu_c_mean, y = 1)) +
          
          geom_vline(xintercept = grand_mean_val,
                     colour = "steelblue", linewidth = 0.6, linetype = "solid") +
          
          geom_vline(xintercept = 0,
                     colour = "grey30", linewidth = 0.5, linetype = "dashed") +
          
          annotate("text", x = grand_mean_val, y = 1.55,
                   label = sprintf("Grand mean = %.1f%%", 100 * grand_mean_val),
                   hjust = -0.05, size = 2.8, colour = "steelblue") +
          
          annotate("text", x = 0, y = 1.55,
                   label = "Zero loss",
                   hjust = 1.05, size = 2.8, colour = "grey30") +
          
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
            aes(fill = mu_b_mean),
            shape  = 21,
            size   = 1.8,
            alpha  = 0.75,
            stroke = 0.25,
            colour = "white",
            height = 0.07,
            width  = 0
          ) +
          scale_fill_gradient(
            low  = "#E31A1C",   # red = low betting
            high = "#2166AC",   # blue = high betting
            name = "Mean betting\nprobability",
            labels = scales::percent_format(accuracy = 1)
          ) +
          
          scale_x_continuous(
            name   = expression(hat(mu)[i]^c ~ "(posterior mean welfare loss, share of endowment)"),
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
            plot.margin = margin(8, 12, 8, 0)  # extra left margin
          )
        
        ggsave(f_rain, p, width = 8, height = 3, dpi = 300)
        msg("Saved: ", f_rain)
      }
    }
  }
  
  invisible(TRUE)
}