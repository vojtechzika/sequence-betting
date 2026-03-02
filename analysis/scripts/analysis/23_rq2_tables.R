# ============================================================
# scripts/analysis/23_rq2_tables.R
# RQ2 tables (exploratory + confirmatory)
# Naming: rq2_<tr>_<type>_<unit>.csv
# ============================================================

library(data.table)
library(rstan)

rq2_tables <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  # Confirmatory treatment = first in run$treatment (same as RQ1)
  confirmatory_tr <- as.character(cfg$run$treatment)[1]
  
  e <- as.numeric(cfg$design$seq$endowment)
  rho_vec <- as.numeric(cfg$design$rq2$rho)
  rho_main <- rho_vec[1]
  
  tau_main <- as.numeric(cfg$design$a_flags$tau[1])
  tau_nm <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (tr in tr_vec) {
    
    f_fit  <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_pid  <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_seq  <- file.path(mod_dir, paste0("rq2_seq_levels_", tr, ".rds"))
    f_prep <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    
    stopifnot(file.exists(f_fit), file.exists(f_pid),
              file.exists(f_seq), file.exists(f_prep))
    
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    prep <- fread(f_prep)
    
    # STRICT schema
    required_cols <- c("pid","seq","delta_bar","sd_star","n_bets")
    if (!all(required_cols %in% names(prep))) {
      stop("Prepared file missing required columns. Rerun 21_rq2_stan after deleting artifacts.")
    }
    
    prep[, pid := as.character(pid)]
    
    pid_map <- prep[, .(
      delta_bar = delta_bar[1],
      sd_star   = sd_star[1],
      n_bets    = n_bets[1]
    ), by = pid]
    
    setkey(pid_map, pid)
    pid_map <- pid_map[.(pid_levels)]
    stopifnot(!anyNA(pid_map$pid))
    
    delta_bar_vec <- pid_map$delta_bar
    sd_star_vec   <- pid_map$sd_star
    n_bets_vec    <- pid_map$n_bets
    
    fit  <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    alpha <- post$alpha
    u     <- post$u
    b     <- post$b
    
    iters <- length(alpha)
    N <- length(pid_levels)
    S <- length(seq_levels)
    
    # =======================
    # SEQUENCES – EXPLORATORY
    # =======================
    
    eta_base <- sweep(u, 1, alpha, "+")
    mu_s_draws <- matrix(NA_real_, iters, S)
    
    for (s in 1:S) {
      eta_mat <- eta_base + b[, s]
      mapped  <- sweep(eta_mat, 2, sd_star_vec, "*")
      mapped  <- sweep(mapped, 2, delta_bar_vec, "+")
      mu_s_draws[, s] <- rowMeans(mapped) / e
    }
    
    seq_tbl <- data.table(
      sequence = seq_levels,
      mu_a_median = apply(mu_s_draws, 2, median),
      mu_a_mean   = apply(mu_s_draws, 2, mean),
      mu_a_q025   = apply(mu_s_draws, 2, quantile, 0.025),
      mu_a_q975   = apply(mu_s_draws, 2, quantile, 0.975)
    )
    
    for (rho in rho_vec) {
      nmU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      nmO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      seq_tbl[, (nmU) := apply(mu_s_draws, 2, function(x) mean(x < -rho))]
      seq_tbl[, (nmO) := apply(mu_s_draws, 2, function(x) mean(x >  rho))]
    }
    
    f_seq_expl <- file.path(out_dir, paste0("rq2_", tr, "_exploratory_sequences.csv"))
    
    if (!should_skip(
      paths = f_seq_expl,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ2 exploratory sequences (", tr, ")")
    )) {
      fwrite(seq_tbl, f_seq_expl)
      msg("Saved: ", f_seq_expl)
    }
    
    # =========================
    # PARTICIPANTS – EXPLORATORY
    # =========================
    
    mu_i_draws <- matrix(NA_real_, iters, N)
    for (i in 1:N) {
      mu_i_draws[, i] <- (delta_bar_vec[i] +
                            (alpha + u[, i]) * sd_star_vec[i]) / e
    }
    
    part_tbl <- data.table(
      pid = pid_levels,
      n_bets = n_bets_vec,
      mu_a_median = apply(mu_i_draws, 2, median),
      mu_a_mean   = apply(mu_i_draws, 2, mean),
      mu_a_q025   = apply(mu_i_draws, 2, quantile, 0.025),
      mu_a_q975   = apply(mu_i_draws, 2, quantile, 0.975)
    )
    
    f_part_expl <- file.path(out_dir, paste0("rq2_", tr, "_exploratory_participants.csv"))
    
    if (!should_skip(
      paths = f_part_expl,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ2 exploratory participants (", tr, ")")
    )) {
      fwrite(part_tbl, f_part_expl)
      msg("Saved: ", f_part_expl)
    }
    
    # =========================
    # CONFIRMATORY
    # =========================
    
    if (tr == confirmatory_tr) {
      
      f_flags <- file.path(mod_dir, paste0("a_star_pid_flags_", tr, ".rds"))
      stopifnot(file.exists(f_flags))
      
      flags <- readRDS(f_flags)
      keep_name <- paste0("pid_keep_tau", tau_nm)
      
      if (!(keep_name %in% names(flags$pid_sets))) {
        stop("Missing confirmatory pid set in flags RDS.")
      }
      
      keep_pid <- intersect(pid_levels, flags$pid_sets[[keep_name]])
      keep_idx <- match(keep_pid, pid_levels)
      
      # Sequences confirmatory
      mu_s_conf <- matrix(NA_real_, iters, S)
      for (s in 1:S) {
        eta_mat <- eta_base[, keep_idx, drop = FALSE] + b[, s]
        mapped  <- sweep(eta_mat, 2, sd_star_vec[keep_idx], "*")
        mapped  <- sweep(mapped, 2, delta_bar_vec[keep_idx], "+")
        mu_s_conf[, s] <- rowMeans(mapped) / e
      }
      
      seq_conf <- data.table(
        sequence = seq_levels,
        mu_a_median = apply(mu_s_conf, 2, median),
        mu_a_mean   = apply(mu_s_conf, 2, mean)
      )
      
      f_seq_conf <- file.path(out_dir, paste0("rq2_", tr, "_confirmatory_sequences.csv"))
      fwrite(seq_conf, f_seq_conf)
      msg("Saved: ", f_seq_conf)
      
      # Participants confirmatory
      part_conf <- part_tbl[pid %in% keep_pid]
      f_part_conf <- file.path(out_dir, paste0("rq2_", tr, "_confirmatory_participants.csv"))
      fwrite(part_conf, f_part_conf)
      msg("Saved: ", f_part_conf)
    }
  }
  
  invisible(TRUE)
}