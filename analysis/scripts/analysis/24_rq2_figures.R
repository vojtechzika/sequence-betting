# ============================================================
# 24_rq2_figures.R
#
# PURPOSE
#   Produces paper-ready figures for RQ2 (intensive margin of
#   betting): sequence-level forest plot and participant-level
#   raincloud plot. Computes sequence table directly from fit
#   objects and saves forest data for combined multi-RQ plot.
#
# INPUT
#   path_mod/rq2_fit_sequences_<tr>_<tag>[_bb].rds
#   path_mod/rq2_pid_levels_<tr>_<tag>[_bb].rds
#   path_mod/rq2_seq_levels_<tr>_<tag>[_bb].rds
#   path_out/rq2_prepared_<tr>_<tag>[_bb].csv
#   path_out/rq2_<tr>_<tag>[_bb]_participants.csv
#   path_out/rq2_<tr>_<tag>[_bb]_model_summary.csv
#
# OUTPUT
#   path_fig/rq2_<tr>_<tag>[_bb]_sequences_forest.png
#   path_fig/rq2_<tr>_<tag>[_bb]_participants_raincloud.png
#   path_out/rq2_<tr>_<tag>[_bb]_forest_data.csv
#
# TAGS
#   full          -- all participants
#   confirmatory  -- normative betters only
#
# ROBUSTNESS
#   Called with robustness = TRUE to produce _bb outputs.
# ============================================================

rq2_figures <- function(cfg, robustness = FALSE) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  suffix <- if (robustness) "_bb" else ""
  tags   <- c("full", "confirmatory")
  
  e <- as.numeric(cfg$design$seq$endowment)
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      file_stem <- paste0("rq2_", tr, "_", tag, suffix)
      
      f_fit  <- file.path(path_mod, paste0("rq2_fit_sequences_", tr, "_", tag, suffix, ".rds"))
      f_pid  <- file.path(path_mod, paste0("rq2_pid_levels_",    tr, "_", tag, suffix, ".rds"))
      f_seq  <- file.path(path_mod, paste0("rq2_seq_levels_",    tr, "_", tag, suffix, ".rds"))
      f_prep <- file.path(path_mod, paste0("rq2_prepared_", tr, "_", tag, suffix, ".rds"))
      f_pid_csv <- file.path(path_out, paste0(file_stem, "_participants.csv"))
      f_mod_csv <- file.path(path_out, paste0(file_stem, "_model_summary.csv"))
      
      if (!all(file.exists(c(f_fit, f_pid, f_seq, f_prep, f_pid_csv, f_mod_csv)))) {
        warning("RQ2 figures: missing inputs for tr='", tr, "', tag='", tag, "'. Skipping.")
        next
      }
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      prep <- readRDS(f_prep)
      pid_tbl    <- fread(f_pid_csv, encoding = "UTF-8")
      mod_tbl    <- fread(f_mod_csv, encoding = "UTF-8")
      
      # ---- Extract posterior draws ----
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      
      alpha <- as.numeric(post$alpha)
      u     <- post$u
      b     <- post$b
      K     <- length(alpha)
      N     <- length(pid_levels)
      S     <- length(seq_levels)
      
      # ---- pid constants from prep ----
      prep[, pid := as.character(pid)]
      pid_map <- prep[, .(delta_bar = delta_bar[1], sd_star = sd_star[1]), by = pid]
      setkey(pid_map, pid)
      pid_map       <- pid_map[.(pid_levels)]
      delta_bar_vec <- as.numeric(pid_map$delta_bar)
      sd_star_vec   <- as.numeric(pid_map$sd_star)
      
      # ---- Compute mu_s draws ----
      eta_base   <- sweep(u, 1, alpha, "+")  # K x N
      mu_s_draws <- matrix(NA_real_, nrow = K, ncol = S)
      for (s in seq_len(S)) {
        eta_mat        <- eta_base + b[, s]
        mapped         <- sweep(eta_mat, 2, sd_star_vec, "*")
        mapped         <- sweep(mapped,  2, delta_bar_vec, "+")
        mu_s_draws[, s] <- rowMeans(mapped) / e
      }
      
      mu_i_draws <- (sweep(eta_base, 2, sd_star_vec, "*") +
                       rep(delta_bar_vec, each = K)) / e
      
      # ---- Grand mean ----
      grand_draws    <- apply(mu_s_draws, 1, mean)
      grand_mean_val <- median(grand_draws)
      grand_lo       <- quantile(grand_draws, 0.025)
      grand_hi       <- quantile(grand_draws, 0.975)
      
      # ---- Sequence table ----
      rho_vec     <- as.numeric(cfg$design$rq2$rho)
      rho_main    <- rho_vec[1]
      rho_nm_main <- gsub("\\.", "", sprintf("%.2f", rho_main))
      
      n_trials_by_seq <- prep[pid %in% pid_levels, .(n_trials = .N), by = seq]
      setkey(n_trials_by_seq, seq)
      
      seq_tbl <- data.table(
        sequence    = seq_levels,
        n_trials    = as.integer(n_trials_by_seq[.(seq_levels), n_trials]),
        mu_a_median = apply(mu_s_draws, 2, median),
        mu_a_mean   = apply(mu_s_draws, 2, mean),
        mu_a_q025   = apply(mu_s_draws, 2, quantile, probs = 0.025),
        mu_a_q975   = apply(mu_s_draws, 2, quantile, probs = 0.975)
      )
      seq_tbl[is.na(n_trials), n_trials := 0L]
      
      for (rho in rho_vec) {
        nmU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        nmO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        seq_tbl[, (nmU) := apply(mu_s_draws, 2, function(x) mean(x < -rho))]
        seq_tbl[, (nmO) := apply(mu_s_draws, 2, function(x) mean(x >  rho))]
      }
      
      colU <- paste0("U_a_rho_", rho_nm_main)
      colO <- paste0("O_a_rho_", rho_nm_main)
      U    <- seq_tbl[[colU]]
      O    <- seq_tbl[[colO]]
      
      seq_tbl[, calib_label := {
        lab        <- rep("neutral", .N)
        pick_under <- (U >= 0.50 & O < 0.50) | (U >= 0.50 & O >= 0.50 & U >= O)
        pick_over  <- (O >= 0.50 & U < 0.50) | (U >= 0.50 & O >= 0.50 & O >  U)
        lab[pick_under & U >= 0.95] <- "under_strong"
        lab[pick_under & U >= 0.80 & U < 0.95] <- "under_moderate"
        lab[pick_under & U >= 0.50 & U < 0.80] <- "under_weak"
        lab[pick_over  & O >= 0.95] <- "over_strong"
        lab[pick_over  & O >= 0.80 & O < 0.95] <- "over_moderate"
        lab[pick_over  & O >= 0.50 & O < 0.80] <- "over_weak"
        lab
      }]
      
      seq_tbl[, grand_label := fifelse(
        mu_a_q025 > grand_mean_val, "above",
        fifelse(mu_a_q975 < grand_mean_val, "below", "neutral")
      )]
      
      setorder(seq_tbl, sequence)
      
      # ---- Forest data export ----
      f_forest_data <- file.path(path_out, paste0(file_stem, "_forest_data.csv"))
      if (!should_skip(f_forest_data, cfg, "output",
                       paste0("RQ2 forest data (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        forest_data <- seq_tbl[, .(
          sequence    = sequence,
          mu_mean     = mu_a_mean,
          mu_median   = mu_a_median,
          mu_q025     = mu_a_q025,
          mu_q975     = mu_a_q975,
          calib_label = calib_label,
          grand_label = grand_label,
          grand_mean  = grand_mean_val,
          grand_lo    = as.numeric(grand_lo),
          grand_hi    = as.numeric(grand_hi),
          treatment   = tr,
          tag         = tag,
          rq          = "rq2"
        )]
        fwrite(forest_data, f_forest_data)
        msg("Saved: ", f_forest_data)
      }
      
      # ---- Forest plot ----
      f_forest <- file.path(path_fig, paste0(file_stem, "_sequences_forest.png"))
      if (!should_skip(f_forest, cfg, "output",
                       paste0("RQ2 forest plot (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        dt <- copy(seq_tbl)
        dt[, sequence    := factor(sequence, levels = dt[order(mu_a_mean), sequence])]
        dt[, grand_label := factor(grand_label, levels = c("above", "neutral", "below"))]
        
        label_colors <- c(above = "#C0392B", neutral = "#AAAAAA", below = "#2166AC")
        label_names  <- c(above = "Above grand mean", neutral = "Neutral",
                          below = "Below grand mean")
        
        p <- ggplot(dt, aes(y = sequence)) +
          
          annotate("rect",
                   xmin = as.numeric(grand_lo), xmax = as.numeric(grand_hi),
                   ymin = -Inf, ymax = Inf,
                   fill = "steelblue", alpha = 0.10) +
          
          geom_vline(xintercept = grand_mean_val,
                     colour = "steelblue", linewidth = 0.6, linetype = "solid") +
          
          geom_vline(xintercept = 0,
                     colour = "grey30", linewidth = 0.5, linetype = "dashed") +
          
          geom_segment(aes(x = mu_a_q025, xend = mu_a_q975,
                           y = sequence,  yend = sequence,
                           colour = grand_label),
                       linewidth = 0.5, alpha = 0.7) +
          
          geom_point(aes(x = mu_a_mean, colour = grand_label), size = 1.8) +
          
          scale_colour_manual(values = label_colors, name = NULL, labels = label_names) +
          
          scale_x_continuous(
            name   = expression(hat(mu)[s]^a ~ "(posterior mean proportional stake deviation)"),
            labels = scales::percent_format(accuracy = 1)
          ) +
          
          scale_y_discrete(name = NULL) +
          
          annotate("text", x = grand_mean_val, y = 1,
                   label = sprintf("Grand mean = %.1f%%", 100 * grand_mean_val),
                   hjust = -0.05, size = 2.8, colour = "steelblue") +
          
          annotate("text", x = 0, y = 1,
                   label = "EU-optimal",
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
      
      # ---- Participant raincloud ----
      f_rain <- file.path(path_fig, paste0(file_stem, "_participants_raincloud.png"))
      if (!should_skip(f_rain, cfg, "output",
                       paste0("RQ2 raincloud (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        p <- ggplot(pid_tbl, aes(x = mu_a_mean, y = 1)) +
          
          geom_vline(xintercept = grand_mean_val,
                     colour = "steelblue", linewidth = 0.6, linetype = "solid") +
          
          geom_vline(xintercept = 0,
                     colour = "grey30", linewidth = 0.5, linetype = "dashed") +
          
          annotate("text", x = grand_mean_val, y = 1.55,
                   label = sprintf("Grand mean = %.1f%%", 100 * grand_mean_val),
                   hjust = -0.05, size = 2.8, colour = "steelblue") +
          
          annotate("text", x = 0, y = 1.55,
                   label = "EU-optimal",
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
            name   = expression(hat(mu)[i]^a ~ "(posterior mean proportional stake deviation)"),
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