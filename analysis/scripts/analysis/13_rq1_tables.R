rq1_tables <- function(cfg) {
  
  ds    <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  # thresholds
  rho_vec <- as.numeric(cfg$design$rq1$rho)
  rho_main <- rho_vec[1]
  
  # master
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  master <- data.table::fread(f_master, encoding = "UTF-8")
  stopifnot(all(c("pid","treat","stake","seq") %in% names(master)))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, y := as.integer(stake > 0)]
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  tags <- c("full", "conf")
  outputs <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      # if you only produce conf for betting_normative==TRUE, optionally skip early:
      if (tag == "conf" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, "_", tag, ".rds"))
      f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, "_", tag, ".rds"))
      f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, "_", tag, ".rds"))
      
      if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) next
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      stopifnot(length(seq_levels) == 64L)
      
      # data slice used for n_trials
      d <- master[treat == tr & pid %in% pid_levels]
      if (nrow(d) == 0) next
      
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      
      stopifnot(!is.null(post$mu_b), !is.null(post$mu_b_i))
      mu_b_draws  <- post$mu_b    # iters x S
      mu_bi_draws <- post$mu_b_i  # iters x N
      
      stopifnot(is.matrix(mu_b_draws),  ncol(mu_b_draws)  == length(seq_levels))
      stopifnot(is.matrix(mu_bi_draws), ncol(mu_bi_draws) == length(pid_levels))
      
      # map tag -> output type label
      type_nm <- if (tag == "full") "exploratory" else "confirmatory"
      
      # --------------------
      # sequences table
      # --------------------
      f_seq_csv <- file.path(out_dir, paste0("rq1_", tr, "_", type_nm, "_sequences.csv"))
      if (!should_skip(f_seq_csv, cfg, "output",
                       paste0("RQ1 sequences ", type_nm, " (", ds, "/", tr, ")"))) {
        
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
          nm <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
          thr <- 1 - rho
          seq_tbl[, (nm) := apply(mu_b_draws, 2, function(x) mean(x < thr))]
        }
        
        col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
        seq_tbl[, underbet_label :=
                  fifelse(get(col_main) >= 0.95, "strong",
                          fifelse(get(col_main) >= 0.80, "moderate",
                                  fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
        
        setorder(seq_tbl, sequence)
        fwrite(seq_tbl, f_seq_csv)
        msg("Saved: ", f_seq_csv)
      }
      
      # --------------------
      # participants table
      # --------------------
      f_pid_csv <- file.path(out_dir, paste0("rq1_", tr, "_", type_nm, "_participants.csv"))
      if (!should_skip(f_pid_csv, cfg, "output",
                       paste0("RQ1 participants ", type_nm, " (", ds, "/", tr, ")"))) {
        
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
          nm <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
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
      }
      
      outputs[[paste(tr, tag, sep = "_")]] <- list(
        fit = f_fit, pid = f_pid, seq = f_seq,
        sequences_csv = f_seq_csv,
        participants_csv = f_pid_csv
      )
    }
  }
  
  invisible(outputs)
}