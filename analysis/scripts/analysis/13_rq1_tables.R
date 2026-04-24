# ============================================================
# 13_rq1_tables.R
#
# PURPOSE
#   Produces sequence-level, participant-level, and model summary
#   tables for RQ1 (extensive margin of betting).
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/rq1_fit_sequences_<tr>_<tag>[_bb].rds
#   path_mod/rq1_pid_levels_<tr>_<tag>[_bb].rds
#   path_mod/rq1_seq_levels_<tr>_<tag>[_bb].rds
#
# OUTPUT
#   path_out/rq1_<tr>_<tag>[_bb]_sequences.csv
#   path_out/rq1_<tr>_<tag>[_bb]_sequences_summary.csv
#   path_out/rq1_<tr>_<tag>[_bb]_participants.csv
#   path_out/rq1_<tr>_<tag>[_bb]_participants_summary.csv
#   path_out/rq1_<tr>_<tag>[_bb]_model_summary.csv
#
# TAGS
#   full          -- all participants
#   confirmatory  -- normative betters only
#
# ROBUSTNESS
#   Called with robustness = TRUE to produce _bb outputs from
#   Beta-Binomial fits. Mirrors primary output structure exactly.
# ============================================================

rq1_tables <- function(cfg, robustness = FALSE) {
  
  tr_vec   <- unique(as.character(cfg$run$treatment))
  rho_vec  <- as.numeric(cfg$design$rq1$rho)
  rho_main <- rho_vec[1]
  suffix   <- if (robustness) "_bb" else ""
  
  f_master <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  master <- fread(f_master, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "seq") %in% names(master)))
  
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, y := as.integer(stake > 0)]
  
  tags    <- c("full", "confirmatory")
  outputs <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      f_fit <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_", tag, suffix, ".rds"))
      f_pid <- file.path(path_mod, paste0("rq1_pid_levels_",    tr, "_", tag, suffix, ".rds"))
      f_seq <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_", tag, suffix, ".rds"))
      
      if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) next
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      stopifnot(length(seq_levels) == 64L)
      
      d <- master[treat == tr & pid %in% pid_levels]
      if (nrow(d) == 0) next
      
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      stopifnot(!is.null(post$mu_b), !is.null(post$mu_b_i))
      
      mu_b_draws  <- post$mu_b
      mu_bi_draws <- post$mu_b_i
      stopifnot(is.matrix(mu_b_draws),  ncol(mu_b_draws)  == length(seq_levels))
      stopifnot(is.matrix(mu_bi_draws), ncol(mu_bi_draws) == length(pid_levels))
      
      file_stem <- paste0("rq1_", tr, "_", tag, suffix)
      
      # ----------------------------------------
      # Sequence table
      # ----------------------------------------
      f_seq_csv <- file.path(path_out, paste0(file_stem, "_sequences.csv"))
      if (!should_skip(f_seq_csv, cfg, "output",
                       paste0("RQ1 sequences (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
        setkey(n_trials_by_seq, seq)
        
        seq_tbl <- data.table(
          sequence    = seq_levels,
          mu_b_median = apply(mu_b_draws, 2, median),
          mu_b_mean   = apply(mu_b_draws, 2, mean),
          mu_b_q025   = apply(mu_b_draws, 2, quantile, probs = 0.025),
          mu_b_q975   = apply(mu_b_draws, 2, quantile, probs = 0.975)
        )
        
        seq_tbl[, n_trials := as.integer(n_trials_by_seq[.(sequence), n_trials])]
        
        for (rho in rho_vec) {
          nm  <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
          thr <- 1 - rho
          seq_tbl[, (nm) := apply(mu_b_draws, 2, function(x) mean(x < thr))]
        }
        
        col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
        seq_tbl[, underbet_label :=
                  fifelse(get(col_main) >= 0.95, "strong",
                          fifelse(get(col_main) >= 0.80, "moderate",
                                  fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
        
        stopifnot(!is.null(post$grand_mean_bet))
        grand_draws    <- post$grand_mean_bet
        grand_mean_val <- median(grand_draws)
        
        seq_tbl[, p_above_grand := apply(mu_b_draws, 2,
                                         function(x) mean(x > grand_draws))]
        seq_tbl[, p_below_grand := apply(mu_b_draws, 2,
                                         function(x) mean(x < grand_draws))]
        seq_tbl[, grand_label := fifelse(
          mu_b_q025 > grand_mean_val, "above",
          fifelse(mu_b_q975 < grand_mean_val, "below", "neutral")
        )]
        
        setorder(seq_tbl, sequence)
        fwrite(seq_tbl, f_seq_csv)
        msg("Saved: ", f_seq_csv)
        
        # Sequence summary
        f_seq_sum <- file.path(path_out, paste0(file_stem, "_sequences_summary.csv"))
        if (!should_skip(f_seq_sum, cfg, "output",
                         paste0("RQ1 sequences summary (", tr, "/", tag,
                                if (robustness) "/bb" else "", ")"))) {
          
          seq_prereg <- seq_tbl[, .N, by = underbet_label]
          seq_prereg[, underbet_label := factor(underbet_label,
                                                levels = c("strong", "moderate", "weak", "neutral"))]
          setorder(seq_prereg, underbet_label)
          seq_prereg[, pct := round(100 * N / nrow(seq_tbl), 1)]
          
          seq_prereg_stats <- seq_tbl[, .(
            mu_mean   = round(mean(mu_b_mean),   3),
            mu_median = round(median(mu_b_mean), 3),
            mu_min    = round(min(mu_b_mean),    3),
            mu_max    = round(max(mu_b_mean),    3)
          ), by = underbet_label]
          
          seq_prereg_ids <- seq_tbl[, .(
            sequences = paste(sort(sequence), collapse = ", ")
          ), by = underbet_label]
          
          seq_prereg <- merge(seq_prereg, seq_prereg_stats, by = "underbet_label")
          seq_prereg <- merge(seq_prereg, seq_prereg_ids,   by = "underbet_label")
          seq_prereg[, classification := "preregistered"]
          
          seq_grand <- seq_tbl[, .N, by = grand_label]
          seq_grand[, grand_label := factor(grand_label,
                                            levels = c("above", "neutral", "below"))]
          setorder(seq_grand, grand_label)
          seq_grand[, pct := round(100 * N / nrow(seq_tbl), 1)]
          
          seq_grand_stats <- seq_tbl[, .(
            mu_mean   = round(mean(mu_b_mean),   3),
            mu_median = round(median(mu_b_mean), 3),
            mu_min    = round(min(mu_b_mean),    3),
            mu_max    = round(max(mu_b_mean),    3)
          ), by = grand_label]
          
          seq_grand_ids <- seq_tbl[, .(
            sequences = paste(sort(sequence), collapse = ", ")
          ), by = grand_label]
          
          seq_grand <- merge(seq_grand, seq_grand_stats, by = "grand_label")
          seq_grand <- merge(seq_grand, seq_grand_ids,   by = "grand_label")
          setnames(seq_grand, "grand_label", "underbet_label")
          seq_grand[, classification := "grand_mean"]
          
          seq_summary <- rbindlist(list(seq_prereg, seq_grand), fill = TRUE)
          fwrite(seq_summary, f_seq_sum)
          msg("Saved: ", f_seq_sum)
        }
      }
      
      # ----------------------------------------
      # Participant table
      # ----------------------------------------
      f_pid_csv <- file.path(path_out, paste0(file_stem, "_participants.csv"))
      if (!should_skip(f_pid_csv, cfg, "output",
                       paste0("RQ1 participants (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        n_trials_by_pid <- d[, .(n_trials = .N), by = pid]
        setkey(n_trials_by_pid, pid)
        
        pid_tbl <- data.table(
          pid         = pid_levels,
          mu_i_median = apply(mu_bi_draws, 2, median),
          mu_i_mean   = apply(mu_bi_draws, 2, mean),
          mu_i_q025   = apply(mu_bi_draws, 2, quantile, probs = 0.025),
          mu_i_q975   = apply(mu_bi_draws, 2, quantile, probs = 0.975)
        )
        
        for (rho in rho_vec) {
          nm  <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
          thr <- 1 - rho
          pid_tbl[, (nm) := apply(mu_bi_draws, 2, function(x) mean(x < thr))]
        }
        
        col_main <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
        pid_tbl[, underbet_label :=
                  fifelse(get(col_main) >= 0.95, "solid",
                          fifelse(get(col_main) >= 0.90, "likely",
                                  fifelse(get(col_main) >= 0.75, "leaning", "neutral")))]
        
        pid_tbl[, n_trials := as.integer(n_trials_by_pid[.(pid), n_trials])]
        setorder(pid_tbl, pid)
        fwrite(pid_tbl, f_pid_csv)
        msg("Saved: ", f_pid_csv)
        
        # Participant summary
        f_pid_sum <- file.path(path_out, paste0(file_stem, "_participants_summary.csv"))
        if (!should_skip(f_pid_sum, cfg, "output",
                         paste0("RQ1 participants summary (", tr, "/", tag,
                                if (robustness) "/bb" else "", ")"))) {
          
          pid_summary <- pid_tbl[, .N, by = underbet_label]
          pid_summary[, underbet_label := factor(underbet_label,
                                                 levels = c("solid", "likely", "leaning", "neutral"))]
          setorder(pid_summary, underbet_label)
          pid_summary[, pct := round(100 * N / nrow(pid_tbl), 1)]
          
          pid_stats <- pid_tbl[, .(
            mu_mean   = round(mean(mu_i_mean),   3),
            mu_median = round(median(mu_i_mean), 3),
            mu_min    = round(min(mu_i_mean),    3),
            mu_max    = round(max(mu_i_mean),    3)
          ), by = underbet_label]
          
          pid_ids <- pid_tbl[, .(
            pids = paste(sort(pid), collapse = ", ")
          ), by = underbet_label]
          
          pid_summary <- merge(pid_summary, pid_stats, by = "underbet_label")
          pid_summary <- merge(pid_summary, pid_ids,   by = "underbet_label")
          fwrite(pid_summary, f_pid_sum)
          msg("Saved: ", f_pid_sum)
        }
      }
      
      # ----------------------------------------
      # Model summary
      # ----------------------------------------
      f_mod_csv <- file.path(path_out, paste0(file_stem, "_model_summary.csv"))
      if (!should_skip(f_mod_csv, cfg, "output",
                       paste0("RQ1 model summary (", tr, "/", tag,
                              if (robustness) "/bb" else "", ")"))) {
        
        stopifnot(!is.null(post$alpha), !is.null(post$sigma_u),
                  !is.null(post$sigma_s), !is.null(post$grand_mean_bet))
        
        sum_draw <- function(x, nm) {
          data.table(
            parameter = nm,
            median    = median(x),
            mean      = mean(x),
            q025      = quantile(x, 0.025),
            q975      = quantile(x, 0.975)
          )
        }
        
        mod_tbl <- rbindlist(list(
          sum_draw(post$alpha,          "alpha (log-odds intercept)"),
          sum_draw(post$grand_mean_bet, "grand mean bet probability"),
          sum_draw(post$sigma_u,        "sigma_u (between-participant SD)"),
          sum_draw(post$sigma_s,        "sigma_s (between-sequence SD)"),
          sum_draw(post$gamma_drift,    "gamma_drift (linear drift per block)")
        ))
        
        mod_tbl[, treatment := tr]
        mod_tbl[, tag       := tag]
        if (robustness) mod_tbl[, likelihood := "beta_binomial"]
        
        fwrite(mod_tbl, f_mod_csv)
        msg("Saved: ", f_mod_csv)
      }
      
      outputs[[paste(tr, tag, sep = "_")]] <- list(
        sequences_csv         = f_seq_csv,
        sequences_summary_csv = f_seq_sum,
        participants_csv      = f_pid_csv,
        participants_summary_csv = f_pid_sum,
        model_summary_csv     = f_mod_csv
      )
    }
  }
  
  invisible(outputs)
}