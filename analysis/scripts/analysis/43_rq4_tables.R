# ============================================================
# scripts/analysis/43_rq4_tables.R
# RQ4 tables (single pass; conditional on betting only)
# Naming:
#   rq4_<tr>_sequences.csv
#   rq4_<tr>_participants.csv
# ============================================================

library(data.table)
library(rstan)

rq4_tables <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  stopifnot(!is.null(design$seq), !is.null(design$seq$side_labels))
  
  # delta grid: main + sensitivities
  delta_vec <- NULL
  if (!is.null(design$rq4) && !is.null(design$rq4$delta)) delta_vec <- design$rq4$delta
  if (is.null(delta_vec)) delta_vec <- c(0.05, 0.03, 0.08)
  delta_vec <- as.numeric(delta_vec)
  stopifnot(length(delta_vec) >= 1L, all(is.finite(delta_vec)), all(delta_vec > 0))
  delta_main <- delta_vec[1]
  delta_nm_main <- gsub("\\.", "", sprintf("%.2f", delta_main))
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads))
  stopifnot(length(lab_tails) == 1L, nzchar(lab_tails))
  
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  master <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid","treat","seq","stake","side") %in% names(master)))
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[, side := as.character(side)]
  master[is.na(stake), stake := 0]
  
  # ---------- label helpers (prereg cutpoints + tie by max prob) ----------
  label_seq <- function(H, T) {
    # neutral cutoff is 0.50 for sequences
    if (max(H, T) < 0.50) return("neutral")
    if (H >= T) {
      if (H >= 0.95) return("strong_head")
      if (H >= 0.80) return("moderate_head")
      return("weak_head")   # here H >= 0.50
    } else {
      if (T >= 0.95) return("strong_tail")
      if (T >= 0.80) return("moderate_tail")
      return("weak_tail")   # here T >= 0.50
    }
  }
  
  label_pid <- function(H, T) {
    # neutral cutoff is 0.75 for participants
    if (max(H, T) < 0.75) return("neutral")
    if (H >= T) {
      if (H >= 0.95) return("solid_headish")
      if (H >= 0.90) return("likely_headish")
      return("leaning_headish")  # here H >= 0.75
    } else {
      if (T >= 0.95) return("solid_tailish")
      if (T >= 0.90) return("likely_tailish")
      return("leaning_tailish")  # here T >= 0.75
    }
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, "_full.rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, "_full.rds"))
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) next
    
    fit <- readRDS(f_fit)
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    stopifnot(length(pid_levels) >= 1L, length(seq_levels) >= 1L)
    
    post <- rstan::extract(fit)
    if (is.null(post$mu_h))  stop("RQ4 fit missing generated quantity mu_h: ", f_fit)
    if (is.null(post$hbar))  stop("RQ4 fit missing generated quantity hbar: ", f_fit)
    if (is.null(post$alpha)) stop("RQ4 fit missing parameter alpha: ", f_fit)
    if (is.null(post$u))     stop("RQ4 fit missing transformed parameter u: ", f_fit)
    if (is.null(post$beta))  stop("RQ4 fit missing transformed parameter beta: ", f_fit)
    
    mu_s_draws <- post$mu_h     # iters x S
    alpha      <- post$alpha    # iters
    u_mat      <- post$u        # iters x N
    b_mat      <- post$beta     # iters x S
    
    stopifnot(is.matrix(mu_s_draws), ncol(mu_s_draws) == length(seq_levels))
    stopifnot(is.matrix(u_mat), ncol(u_mat) == length(pid_levels))
    stopifnot(is.matrix(b_mat), ncol(b_mat) == length(seq_levels))
    stopifnot(length(alpha) == nrow(mu_s_draws),
              length(alpha) == nrow(u_mat),
              length(alpha) == nrow(b_mat))
    
    iters <- length(alpha)
    N <- length(pid_levels)
    S <- length(seq_levels)
    
    # baseline (prereg): use Stan-generated hbar = mean_i inv_logit(alpha + u_i)
    if (is.null(post$hbar)) stop("RQ4 fit missing generated quantity hbar: ", f_fit)
    hbar_draws <- as.numeric(post$hbar)   # iters
    stopifnot(length(hbar_draws) == length(alpha))
    
    # counts (betting trials only; align to fit levels)
    d_counts <- master[
      treat == tr &
        pid %in% pid_levels &
        seq %in% seq_levels &
        is.finite(stake) & stake > 0 &
        side %in% c(lab_heads, lab_tails)
    ]
    n_trials_by_seq <- d_counts[, .(n_trials = .N), by = seq]
    setkey(n_trials_by_seq, seq)
    
    # ------------------------------------------------------------
    # SEQUENCES
    # ------------------------------------------------------------
    seq_tbl <- data.table(
      sequence    = seq_levels,
      n_trials    = as.integer(n_trials_by_seq[.(seq_levels), n_trials]),
      mu_h_median = apply(mu_s_draws, 2, median),
      mu_h_mean   = apply(mu_s_draws, 2, mean),
      mu_h_q025   = apply(mu_s_draws, 2, quantile, probs = 0.025),
      mu_h_q975   = apply(mu_s_draws, 2, quantile, probs = 0.975),
      hbar_median = median(hbar_draws),
      hbar_mean   = mean(hbar_draws),
      hbar_q025   = as.numeric(quantile(hbar_draws, 0.025)),
      hbar_q975   = as.numeric(quantile(hbar_draws, 0.975))
    )
    seq_tbl[is.na(n_trials), n_trials := 0L]
    
    # store main H/T to use for labeling (delta_main)
    H_main <- numeric(S)
    T_main <- numeric(S)
    
    for (dlt in delta_vec) {
      nm <- gsub("\\.", "", sprintf("%.2f", dlt))
      nmH <- paste0("H_delta_", nm)
      nmT <- paste0("T_delta_", nm)
      nmA <- paste0("P_absdev05_delta_", nm)
      
      H <- apply(mu_s_draws, 2, function(x) mean(x > (hbar_draws + dlt)))
      T <- apply(mu_s_draws, 2, function(x) mean(x < (hbar_draws - dlt)))
      
      seq_tbl[, (nmH) := H]
      seq_tbl[, (nmT) := T]
      seq_tbl[, (nmA) := apply(mu_s_draws, 2, function(x) mean(abs(x - 0.5) > dlt))]
      
      if (isTRUE(all.equal(dlt, delta_main))) {
        H_main <- H
        T_main <- T
      }
    }
    
    seq_tbl[, direction_label := mapply(label_seq, H_main, T_main)]
    
    seq_tbl[, sequence := factor(sequence, levels = seq_levels)]
    setorder(seq_tbl, sequence)
    seq_tbl[, sequence := as.character(sequence)]
    
    f_seq_csv <- file.path(out_dir, paste0("rq4_", tr, "_sequences.csv"))
    if (!should_skip(
      paths = f_seq_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ4 sequences (", ds, "/", tr, ")")
    )) {
      fwrite(seq_tbl, f_seq_csv)
      msg("Saved: ", f_seq_csv)
    }
    
    # ------------------------------------------------------------
    # PARTICIPANTS
    # mu_i^h(k) = mean_s inv_logit(alpha_k + u_ki + beta_ks)
    # ------------------------------------------------------------
    mu_i_draws <- matrix(NA_real_, nrow = iters, ncol = N)
    for (k in seq_len(iters)) {
      eta <- matrix(alpha[k] + u_mat[k, ], nrow = N, ncol = S) +
        matrix(b_mat[k, ], nrow = N, ncol = S, byrow = TRUE)
      mu_i_draws[k, ] <- rowMeans(inv_logit(eta))
    }
    
    part_tbl <- data.table(
      pid         = pid_levels,
      mu_h_median = apply(mu_i_draws, 2, median),
      mu_h_mean   = apply(mu_i_draws, 2, mean),
      mu_h_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
      mu_h_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975),
      hbar_median = median(hbar_draws),
      hbar_mean   = mean(hbar_draws),
      hbar_q025   = as.numeric(quantile(hbar_draws, 0.025)),
      hbar_q975   = as.numeric(quantile(hbar_draws, 0.975))
    )
    
    # store main H/T to label
    H_i_main <- numeric(N)
    T_i_main <- numeric(N)
    
    for (dlt in delta_vec) {
      nm <- gsub("\\.", "", sprintf("%.2f", dlt))
      nmH <- paste0("H_delta_", nm)
      nmT <- paste0("T_delta_", nm)
      nmA <- paste0("P_absdev05_delta_", nm)
      
      H <- apply(mu_i_draws, 2, function(x) mean(x > (hbar_draws + dlt)))
      T <- apply(mu_i_draws, 2, function(x) mean(x < (hbar_draws - dlt)))
      
      part_tbl[, (nmH) := H]
      part_tbl[, (nmT) := T]
      part_tbl[, (nmA) := apply(mu_i_draws, 2, function(x) mean(abs(x - 0.5) > dlt))]
      
      if (isTRUE(all.equal(dlt, delta_main))) {
        H_i_main <- H
        T_i_main <- T
      }
    }
    
    part_tbl[, side_label := mapply(label_pid, H_i_main, T_i_main)]
    
    setorder(part_tbl, pid)
    
    f_pid_csv <- file.path(out_dir, paste0("rq4_", tr, "_participants.csv"))
    if (!should_skip(
      paths = f_pid_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ4 participants (", ds, "/", tr, ")")
    )) {
      fwrite(part_tbl, f_pid_csv)
      msg("Saved: ", f_pid_csv)
    }
    
    outputs[[tr]] <- list(
      sequences_csv = f_seq_csv,
      participants_csv = f_pid_csv,
      sequences = seq_tbl,
      participants = part_tbl
    )
  }
  
  invisible(outputs)
}