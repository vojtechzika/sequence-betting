# ============================================================
# 33_rq3_tables.R
#
# PURPOSE
#   Produces sequence-level, participant-level, and model summary
#   tables for RQ3 (certainty-equivalent welfare loss).
#   Reads rq3_diagnostics.csv to determine selected model.
#   Selected model tables -> path_out/
#   Other model tables    -> path_out/alternatives/
#
# INPUT
#   path_src/master_sequences.csv
#   path_out/rq3_diagnostics.csv
#   path_mod/rq3_fit_sequences_<tr>_<tag>[_gamma|_alt].rds
#   path_mod/rq3_pid_levels_<tr>_<tag>[_gamma|_alt].rds
#   path_mod/rq3_seq_levels_<tr>_<tag>[_gamma|_alt].rds
#
# OUTPUT
#   path_out/rq3_<tr>_<tag>_sequences.csv
#   path_out/rq3_<tr>_<tag>_sequences_summary.csv
#   path_out/rq3_<tr>_<tag>_participants.csv
#   path_out/rq3_<tr>_<tag>_participants_summary.csv
#   path_out/rq3_<tr>_<tag>_model_summary.csv
#   path_out/alternatives/alt_rq3_<tr>_<tag>_*.csv
#
# TAGS
#   full          -- all participants
#   confirmatory  -- normative betters only
# ============================================================

rq3_tables <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  stopifnot(!is.null(design$rq3), !is.null(design$rq3$rho))
  
  rho_vec     <- as.numeric(design$rq3$rho)
  rho_main    <- rho_vec[1]
  rho_nm_main <- gsub("\\.", "", sprintf("%.2f", rho_main))
  
  # ---- Model suffix map ----
  suffix_map <- list(
    primary    = "",
    gamma_only = "_gamma",
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
  
  get_other_suffixes <- function(selected_suffix) {
    all_suffixes <- unlist(suffix_map)
    all_suffixes[all_suffixes != selected_suffix]
  }
  
  # ---- Output dirs ----
  path_out_alt <- file.path(path_out, "alternatives")
  dir.create(path_out_alt, showWarnings = FALSE, recursive = TRUE)
  
  # ---- Load master once ----
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  master <- fread(infile, encoding = "UTF-8")
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  
  tags    <- c("full", "confirmatory")
  outputs <- list()
  
  # ---- Helper: produce tables for one fit ----
  produce_tables <- function(tr, tag, suffix, outdir, is_alt = FALSE) {
    
    f_fit <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tag, suffix, ".rds"))
    f_pid <- file.path(path_mod, paste0("rq3_pid_levels_",    tr, "_", tag, suffix, ".rds"))
    f_seq <- file.path(path_mod, paste0("rq3_seq_levels_",    tr, "_", tag, suffix, ".rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) return(invisible(NULL))
   
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    
    d <- master[treat == tr & pid %in% pid_levels]
    if (nrow(d) == 0) return(invisible(NULL))
    
    fit  <- readRDS(f_fit)
    post <- rstan::extract(fit)
    stopifnot(!is.null(post$mu_c), !is.null(post$mu_c_i))
    
    mu_s_draws <- post$mu_c
    mu_i_draws <- post$mu_c_i
    stopifnot(is.matrix(mu_s_draws), ncol(mu_s_draws) == length(seq_levels))
    stopifnot(is.matrix(mu_i_draws), ncol(mu_i_draws) == length(pid_levels))
    
    file_prefix <- if (is_alt) "alt_" else ""
    file_stem   <- paste0(file_prefix, "rq3_", tr, "_", tag)
    
    grand_draws    <- apply(mu_s_draws, 1, mean)
    grand_mean_val <- median(grand_draws)
    
    sum_draw <- function(x, nm) {
      data.table(parameter = nm, median = median(x), mean = mean(x),
                 q025 = quantile(x, 0.025), q975 = quantile(x, 0.975))
    }
    
    # Determine likelihood label from suffix
    likelihood_label <- if (suffix == "")        "hurdle_gamma" else
      if (suffix == "_gamma")   "gamma_only"   else
        if (suffix == "_alt")     "gaussian"     else
          "unknown"
    
    # ---- Sequence table ----
    f_seq_csv <- file.path(outdir, paste0(file_stem, "_sequences.csv"))
    if (!should_skip(f_seq_csv, cfg, "output",
                     paste0("RQ3 sequences (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      
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
      fwrite(seq_tbl, f_seq_csv)
      msg("Saved: ", f_seq_csv)
      
      # Sequence summary
      f_seq_sum <- file.path(outdir, paste0(file_stem, "_sequences_summary.csv"))
      if (!should_skip(f_seq_sum, cfg, "output",
                       paste0("RQ3 sequences summary (", tr, "/", tag,
                              if (is_alt) "/alt" else "", ")"))) {
        
        seq_prereg <- seq_tbl[, .N, by = loss_label]
        seq_prereg[, loss_label := factor(loss_label,
                                          levels = c("strong", "moderate", "weak", "neutral"))]
        setorder(seq_prereg, loss_label)
        seq_prereg[, pct := round(100 * N / nrow(seq_tbl), 1)]
        
        seq_prereg_stats <- seq_tbl[, .(
          mu_mean   = round(mean(mu_c_mean),   3),
          mu_median = round(median(mu_c_mean), 3),
          mu_min    = round(min(mu_c_mean),    3),
          mu_max    = round(max(mu_c_mean),    3)
        ), by = loss_label]
        
        seq_prereg_ids <- seq_tbl[, .(
          sequences = paste(sort(sequence), collapse = ", ")
        ), by = loss_label]
        
        seq_prereg <- merge(seq_prereg, seq_prereg_stats, by = "loss_label")
        seq_prereg <- merge(seq_prereg, seq_prereg_ids,   by = "loss_label")
        seq_prereg[, classification := "preregistered"]
        
        seq_grand <- seq_tbl[, .N, by = grand_label]
        seq_grand[, grand_label := factor(grand_label,
                                          levels = c("above", "likely_above", "neutral",
                                                     "likely_below", "below"))]
        setorder(seq_grand, grand_label)
        seq_grand[, pct := round(100 * N / nrow(seq_tbl), 1)]
        
        seq_grand_stats <- seq_tbl[, .(
          mu_mean   = round(mean(mu_c_mean),   3),
          mu_median = round(median(mu_c_mean), 3),
          mu_min    = round(min(mu_c_mean),    3),
          mu_max    = round(max(mu_c_mean),    3)
        ), by = grand_label]
        
        seq_grand_ids <- seq_tbl[, .(
          sequences = paste(sort(sequence), collapse = ", ")
        ), by = grand_label]
        
        seq_grand <- merge(seq_grand, seq_grand_stats, by = "grand_label")
        seq_grand <- merge(seq_grand, seq_grand_ids,   by = "grand_label")
        setnames(seq_grand, "grand_label", "loss_label")
        seq_grand[, classification := "grand_mean"]
        
        seq_summary <- rbindlist(list(seq_prereg, seq_grand), fill = TRUE)
        fwrite(seq_summary, f_seq_sum)
        msg("Saved: ", f_seq_sum)
      }
      
    } else {
      seq_tbl <- fread(f_seq_csv)
    }
    
    # ---- Participant table ----
    f_pid_csv <- file.path(outdir, paste0(file_stem, "_participants.csv"))
    if (!should_skip(f_pid_csv, cfg, "output",
                     paste0("RQ3 participants (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      
      n_trials_by_pid <- d[, .(n_trials = .N), by = pid]
      setkey(n_trials_by_pid, pid)
      
      pid_tbl <- data.table(
        pid         = pid_levels,
        mu_c_median = apply(mu_i_draws, 2, median),
        mu_c_mean   = apply(mu_i_draws, 2, mean),
        mu_c_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
        mu_c_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975)
      )
      
      for (rho in rho_vec) {
        nm <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        pid_tbl[, (nm) := apply(mu_i_draws, 2, function(x) mean(x > rho))]
      }
      
      colL <- paste0("L_rho_", rho_nm_main)
      pid_tbl[, loss_label :=
                fifelse(get(colL) >= 0.95, "solid",
                        fifelse(get(colL) >= 0.90, "likely",
                                fifelse(get(colL) >= 0.75, "leaning", "neutral")))]
      
      pid_tbl[, n_trials := as.integer(n_trials_by_pid[.(pid), n_trials])]
      setorder(pid_tbl, pid)
      
      
      # ---- Merge RQ1 betting probability ----
      f_pid_rq1 <- file.path(path_out, paste0("rq1_", tr, "_", tag, "_participants.csv"))
      if (file.exists(f_pid_rq1)) {
        pid_rq1 <- fread(f_pid_rq1, encoding = "UTF-8")
        pid_tbl <- merge(pid_tbl,
                         pid_rq1[, .(pid, mu_b_mean = mu_i_mean)],
                         by = "pid", all.x = TRUE)
      }
      
      #save
      fwrite(pid_tbl, f_pid_csv)
      msg("Saved: ", f_pid_csv)
      
      # Participant summary
      f_pid_sum <- file.path(outdir, paste0(file_stem, "_participants_summary.csv"))
      if (!should_skip(f_pid_sum, cfg, "output",
                       paste0("RQ3 participants summary (", tr, "/", tag,
                              if (is_alt) "/alt" else "", ")"))) {
        
        pid_summary <- pid_tbl[, .N, by = loss_label]
        pid_summary[, loss_label := factor(loss_label,
                                           levels = c("solid", "likely", "leaning", "neutral"))]
        setorder(pid_summary, loss_label)
        pid_summary[, pct := round(100 * N / nrow(pid_tbl), 1)]
        
        pid_stats <- pid_tbl[, .(
          mu_mean   = round(mean(mu_c_mean),   3),
          mu_median = round(median(mu_c_mean), 3),
          mu_min    = round(min(mu_c_mean),    3),
          mu_max    = round(max(mu_c_mean),    3)
        ), by = loss_label]
        
        pid_ids <- pid_tbl[, .(
          pids = paste(sort(pid), collapse = ", ")
        ), by = loss_label]
        
        pid_summary <- merge(pid_summary, pid_stats, by = "loss_label")
        pid_summary <- merge(pid_summary, pid_ids,   by = "loss_label")
        fwrite(pid_summary, f_pid_sum)
        msg("Saved: ", f_pid_sum)
      }
    }
    
    # ---- Model summary ----
    f_mod_csv <- file.path(outdir, paste0(file_stem, "_model_summary.csv"))
    if (!should_skip(f_mod_csv, cfg, "output",
                     paste0("RQ3 model summary (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      
      if (likelihood_label == "hurdle_gamma") {
        mod_tbl <- rbindlist(list(
          sum_draw(post$a0,    "a0 (hurdle intercept)"),
          sum_draw(post$ap,    "ap (positive mean intercept)"),
          sum_draw(post$su0,   "su0 (hurdle between-participant SD)"),
          sum_draw(post$sb0,   "sb0 (hurdle between-sequence SD)"),
          sum_draw(post$sup,   "sup (positive between-participant SD)"),
          sum_draw(post$sbp,   "sbp (positive between-sequence SD)"),
          sum_draw(post$shape, "shape (Gamma shape)")
        ))
      } else if (likelihood_label == "gaussian") {
        mod_tbl <- rbindlist(list(
          sum_draw(post$alpha,   "alpha (grand mean)"),
          sum_draw(post$sigma_u, "sigma_u (between-participant SD)"),
          sum_draw(post$sigma_s, "sigma_s (between-sequence SD)"),
          sum_draw(post$sigma,   "sigma (residual SD)")
        ))
      } else {
        # gamma_only
        mod_tbl <- rbindlist(list(
          sum_draw(post$ap,          "ap (positive mean intercept)"),
          sum_draw(post$sup,         "sup (between-participant SD)"),
          sum_draw(post$sbp,         "sbp (between-sequence SD)"),
          sum_draw(post$shape,       "shape (Gamma shape)"),
          sum_draw(post$gamma_drift, "gamma_drift (linear drift per block)")
        ))
      }
      
      mod_tbl[, treatment  := tr]
      mod_tbl[, tag        := tag]
      mod_tbl[, likelihood := likelihood_label]
      mod_tbl[, n_trials       := nrow(d)]
      mod_tbl[, n_participants := length(pid_levels)]
      mod_tbl[, n_sequences    := length(seq_levels)]

      
      grand_row <- data.table(
        parameter  = "grand mean mu_c (proportion of endowment)",
        median     = median(grand_draws),
        mean       = mean(grand_draws),
        q025       = quantile(grand_draws, 0.025),
        q975       = quantile(grand_draws, 0.975),
        treatment  = tr,
        tag        = tag,
        likelihood = likelihood_label
      )
      
      mod_tbl <- rbindlist(list(mod_tbl, grand_row), fill = TRUE)
      fwrite(mod_tbl, f_mod_csv)
      msg("Saved: ", f_mod_csv)
    }
    
    list(sequences_csv            = f_seq_csv,
         sequences_summary_csv    = f_seq_sum,
         participants_csv         = f_pid_csv,
         participants_summary_csv = f_pid_sum,
         model_summary_csv        = f_mod_csv)
  }
  
  # ============================================================
  # Main loop
  # ============================================================
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(design$a_flags$betting_normative[[tr]])) next
      
      selected_suffix <- get_selected_suffix(tr, tag)
      other_suffixes  <- get_other_suffixes(selected_suffix)
      
      # Selected model -> path_out
      res <- produce_tables(tr, tag,
                            suffix = selected_suffix,
                            outdir = path_out,
                            is_alt = FALSE)
      
      # Other models -> path_out/alternatives
      for (sfx in other_suffixes) {
        produce_tables(tr, tag,
                       suffix = sfx,
                       outdir = path_out_alt,
                       is_alt = TRUE)
      }
      
      if (!is.null(res)) outputs[[paste(tr, tag, sep = "_")]] <- res
    }
  }
  
  invisible(outputs)
}