# scripts/analysis/13_rq1_tables.R
# ============================================================
# RQ1 tables (exploratory + confirmatory) in ONE script
#
# Exploratory (FULL sample; per treatment in cfg$plan$by):
#   - Reads Stan generated quantities:
#       mu_b   (iters x S)   sequence-level mean betting probability
#       mu_b_i (iters x N)   participant-level mean betting probability
#   - Writes:
#       output/rq1_sequences_exploratory_<tr>.csv
#       output/rq1_participants_exploratory_<tr>.csv
#
# Confirmatory (SUBSET; confirmatory treatment only):
#   - Confirmatory treatment = first element of cfg$plan$by
#   - Keep set is derived from:
#       models/a_star_pid_flags_<tr>.rds  (pid_sets)
#     using tau_main = first element of cfg$design$a_flags$tau
#     keep_name = pid_keep_tauXX  (XX from sprintf("%.2f", tau_main) with '.' removed)
#   - Uses FULL RQ1 Stan fit; recomputes within-draw means averaging ONLY over keep set:
#       mu_b_sub[k,s] = mean_{i in keep} inv_logit(alpha[k] + u[k,i] + b[k,s])
#       mu_i_sub[k,j] = mean_s           inv_logit(alpha[k] + u[k,keep_j] + b[k,s])
#   - Writes:
#       output/rq1_sequences_confirmatory_<tr>.csv
#       output/rq1_participants_confirmatory_<tr>.csv
#
# Inputs (per treatment):
#   models/rq1_fit_sequences_<tr>.rds
#   models/rq1_pid_levels_<tr>.rds
#   models/rq1_seq_levels_<tr>.rds
#   clean/<ds>/master_sequences.csv   (for n_trials)
#
# Notes:
#   - Uses should_skip() for all outputs.
#   - No model fitting; post-processing only.
# ============================================================

library(data.table)
library(rstan)

rq1_tables <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  # ----------------------------
  # Rho thresholds
  # ----------------------------
  stopifnot(!is.null(cfg$design$rq1$rho), !is.null(cfg$design$rq1$rho))
  rho_vec <- as.numeric(cfg$design$rq1$rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  rho_main <- rho_vec[1]
  
  # ----------------------------
  # Master (for n_trials)
  # ----------------------------
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  master <- fread(f_master, encoding = "UTF-8")
  req <- c("pid", "treat", "stake", "seq")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0) stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, y := as.integer(stake > 0)]  # RQ1 definition
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list(exploratory = list(), confirmatory = NULL)
  
  # ============================================================
  # EXPLORATORY (full sample) for ALL treatments in cfg$plan$by
  # ============================================================
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, ".rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) {
      stop(
        "Missing RQ1 Stan artifacts for dataset='", ds, "', treatment='", tr, "'.\n",
        "Expected:\n  ", f_fit, "\n  ", f_pid, "\n  ", f_seq, "\n",
        "Run rq1_stan(cfg) first."
      )
    }
    
    d <- master[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ1 tables: No rows for treat='", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    stopifnot(length(pid_levels) >= 1L, length(seq_levels) >= 1L)
    
    # If your design always has 64 sequences per treatment, keep strict:
    stopifnot(length(seq_levels) == 64L)
    
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    if (is.null(post$mu_b))   stop("RQ1 tables: fit missing 'mu_b' (rerun RQ1 Stan).")
    if (is.null(post$mu_b_i)) stop("RQ1 tables: fit missing 'mu_b_i' (rerun RQ1 Stan).")
    
    mu_b_draws  <- post$mu_b    # iters x S
    mu_bi_draws <- post$mu_b_i  # iters x N
    
    stopifnot(is.matrix(mu_b_draws),  ncol(mu_b_draws)  == length(seq_levels))
    stopifnot(is.matrix(mu_bi_draws), ncol(mu_bi_draws) == length(pid_levels))
    
    # ----------------------------
    # SEQUENCES exploratory
    # ----------------------------
    f_seq_csv <- file.path(out_dir, paste0("rq1_", tr, "_exploratory_sequences.csv"))
    
    seq_tbl <- NULL
    if (!should_skip(
      paths = f_seq_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ1 sequences exploratory (", ds, "/", tr, ")")
    )) {
      
      n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
      setkey(n_trials_by_seq, seq)
      
      seq_tbl <- data.table(
        sequence    = seq_levels,
        mu_b_median = apply(mu_b_draws, 2, median),
        mu_b_mean   = apply(mu_b_draws, 2, mean),
        mu_b_q025   = apply(mu_b_draws, 2, quantile, probs = 0.025),
        mu_b_q975   = apply(mu_b_draws, 2, quantile, probs = 0.975)
      )
      
      seq_tbl[, n_trials := n_trials_by_seq[.(sequence), n_trials]]
      seq_tbl[, n_trials := as.integer(n_trials)]
      
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
    
    # ----------------------------
    # PARTICIPANTS exploratory
    # ----------------------------
    f_pid_csv <- file.path(out_dir, paste0("rq1_", tr, "_exploratory_participants.csv"))
    
    pid_tbl <- NULL
    if (!should_skip(
      paths = f_pid_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ1 participants exploratory (", ds, "/", tr, ")")
    )) {
      
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
      
      pid_tbl[, n_trials := n_trials_by_pid[.(pid), n_trials]]
      pid_tbl[, n_trials := as.integer(n_trials)]
      
      setorder(pid_tbl, pid)
      fwrite(pid_tbl, f_pid_csv)
      msg("Saved: ", f_pid_csv)
    }
    
    outputs$exploratory[[tr]] <- list(
      sequences_csv    = f_seq_csv,
      participants_csv = f_pid_csv,
      sequences        = seq_tbl,
      participants     = pid_tbl
    )
  }
  
  # ============================================================
  # CONFIRMATORY (subset) for CONFIRMATORY treatment only
  # - confirmatory treatment = first element of cfg$plan$by
  # ============================================================
  tr_conf <- tr_vec[1]
  
  # tau thresholds (first is main)
  stopifnot(!is.null(cfg$design$a_flags), !is.null(cfg$design$a_flags$tau))
  tau_vec <- as.numeric(cfg$design$a_flags$tau)
  stopifnot(length(tau_vec) >= 1L, all(is.finite(tau_vec)), all(tau_vec > 0), all(tau_vec < 1))
  tau_main <- tau_vec[1]
  tau_nm <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr_conf, ".rds"))
  f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr_conf, ".rds"))
  f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr_conf, ".rds"))
  stopifnot(file.exists(f_fit), file.exists(f_pid), file.exists(f_seq))
  
  f_flags <- file.path(mod_dir, paste0("a_star_pid_flags_", tr_conf, ".rds"))
  stopifnot(file.exists(f_flags))
  
  d_conf <- master[treat == tr_conf]
  if (nrow(d_conf) == 0) stop("RQ1 confirmatory: no rows for treat='", tr_conf, "' in master_sequences.csv.")
  
  f_seq_out <- file.path(out_dir, paste0("rq1_", tr_conf, "_confirmatory_sequences.csv"))
  f_pid_out <- file.path(out_dir, paste0("rq1_", tr_conf, "_confirmatory_participants.csv"))
  
  if (!should_skip(
    paths = c(f_seq_out, f_pid_out),
    cfg   = cfg,
    type  = "output",
    label = paste0("RQ1 confirmatory tables (", ds, "/", tr_conf, ")")
  )) {
    
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    stopifnot(length(seq_levels) == 64L)
    
    # --- derive keep set from flags rds ---
    flags <- readRDS(f_flags)
    stopifnot(is.list(flags), !is.null(flags$pid_sets))
    
    keep_name <- paste0("pid_keep_tau", tau_nm)
    if (!(keep_name %in% names(flags$pid_sets))) {
      stop("a_star_pid_flags_<tr>.rds missing pid set '", keep_name, "'. File=", f_flags)
    }
    
    keep_pid <- as.character(flags$pid_sets[[keep_name]])
    keep_pid <- intersect(pid_levels, keep_pid)  # enforce fit-level alignment
    
    if (length(keep_pid) == 0L) {
      stop("RQ1 confirmatory: keep set is empty after intersecting with pid_levels.\n",
           "treatment=", tr_conf, " | tau_main=", tau_main, "\n",
           "flags file=", f_flags)
    }
    
    keep_idx <- match(keep_pid, pid_levels)
    stopifnot(all(!is.na(keep_idx)))
    
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    stopifnot(!is.null(post$alpha), !is.null(post$u), !is.null(post$b))
    alpha <- as.numeric(post$alpha)  # K
    u     <- post$u                  # K x N
    b     <- post$b                  # K x S
    
    K <- length(alpha)
    N <- length(pid_levels)
    S <- length(seq_levels)
    
    stopifnot(is.matrix(u), is.matrix(b))
    stopifnot(nrow(u) == K, ncol(u) == N)
    stopifnot(nrow(b) == K, ncol(b) == S)
    
    # ----------------------------
    # SEQUENCES confirmatory (subset-averaged)
    # mu_b_sub[k,s] = mean_{i in keep} inv_logit(alpha[k] + u[k,i] + b[k,s])
    # ----------------------------
    mu_b_sub <- matrix(NA_real_, nrow = K, ncol = S)
    
    for (k in seq_len(K)) {
      eta_i <- alpha[k] + u[k, keep_idx]      # |keep|
      eta_mat <- outer(eta_i, b[k, ], "+")    # |keep| x S
      mu_b_sub[k, ] <- colMeans(plogis(eta_mat))
    }
    
    n_trials_by_seq <- d_conf[, .(n_trials = .N), by = seq]
    setkey(n_trials_by_seq, seq)
    
    seq_tbl_c <- data.table(
      sequence    = seq_levels,
      mu_b_median = apply(mu_b_sub, 2, median),
      mu_b_mean   = apply(mu_b_sub, 2, mean),
      mu_b_q025   = apply(mu_b_sub, 2, quantile, probs = 0.025),
      mu_b_q975   = apply(mu_b_sub, 2, quantile, probs = 0.975)
    )
    
    seq_tbl_c[, n_trials := n_trials_by_seq[.(sequence), n_trials]]
    seq_tbl_c[, n_trials := as.integer(n_trials)]
    
    for (rho in rho_vec) {
      nm <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      thr <- 1 - rho
      seq_tbl_c[, (nm) := apply(mu_b_sub, 2, function(x) mean(x < thr))]
    }
    
    col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    seq_tbl_c[, underbet_label :=
                fifelse(get(col_main) >= 0.95, "strong",
                        fifelse(get(col_main) >= 0.80, "moderate",
                                fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
    
    setorder(seq_tbl_c, sequence)
    
    # ----------------------------
    # PARTICIPANTS confirmatory (keep participants only)
    # mu_i_sub[k,j] = mean_s inv_logit(alpha[k] + u[k,keep_j] + b[k,s])
    # ----------------------------
    mu_i_sub <- matrix(NA_real_, nrow = K, ncol = length(keep_idx))
    
    for (j in seq_along(keep_idx)) {
      i <- keep_idx[j]
      eta_mat <- b
      eta_mat <- sweep(eta_mat, 1, alpha + u[, i], "+")  # K x S
      mu_i_sub[, j] <- rowMeans(plogis(eta_mat))
    }
    
    n_trials_by_pid <- d_conf[, .(n_trials = .N), by = pid]
    setkey(n_trials_by_pid, pid)
    
    pid_tbl_c <- data.table(
      pid         = keep_pid,
      mu_i_median = apply(mu_i_sub, 2, median),
      mu_i_mean   = apply(mu_i_sub, 2, mean),
      mu_i_q025   = apply(mu_i_sub, 2, quantile, probs = 0.025),
      mu_i_q975   = apply(mu_i_sub, 2, quantile, probs = 0.975)
    )
    
    for (rho in rho_vec) {
      nm <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      thr <- 1 - rho
      pid_tbl_c[, (nm) := apply(mu_i_sub, 2, function(x) mean(x < thr))]
    }
    
    col_main_i <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    pid_tbl_c[, underbet_label :=
                fifelse(get(col_main_i) >= 0.95, "solid",
                        fifelse(get(col_main_i) >= 0.90, "likely",
                                fifelse(get(col_main_i) >= 0.75, "leaning", "neutral")))]
    
    pid_tbl_c[, n_trials := n_trials_by_pid[.(pid), n_trials]]
    pid_tbl_c[, n_trials := as.integer(n_trials)]
    setorder(pid_tbl_c, pid)
    
    fwrite(seq_tbl_c, f_seq_out)
    msg("Saved: ", f_seq_out)
    
    fwrite(pid_tbl_c, f_pid_out)
    msg("Saved: ", f_pid_out)
    
    outputs$confirmatory <- list(
      dataset = ds,
      confirmatory_treatment = tr_conf,
      tau_main = tau_main,
      rho_main = rho_main,
      subset_n = length(keep_pid),
      sequences_csv = f_seq_out,
      participants_csv = f_pid_out,
      sequences = seq_tbl_c,
      participants = pid_tbl_c
    )
  }
  
  invisible(outputs)
}

# Example:
# source(here::here("scripts","analysis","13_rq1_tables.R"))
# rq1_tables(cfg)