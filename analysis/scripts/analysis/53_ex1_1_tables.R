# ============================================================
# 53_ex1_1_tables.R
#
# EX1.I tables: sequences + participants (draw-level chi).
#
# For each posterior draw k and epsilon:
#   I_a(k,s) = 1{|theta_s^(k) - theta_a^(k)| < eps}
#   w_a(k,s) = (I_a(k,s) + eta) / (I_H + I_T + I_0 + 3*eta)
#   z_is^(k) = sign_is * (wH(k,s) - wT(k,s))
#   chi_s^(k) = mean_{trials on s} z_is^(k)
#   chi_i^(k) = mean_{trials of i} z_is^(k)
#
# Also reports collapsed scalar similarities d_a(s) and w_a(s).
#
# INPUTS
#   path_mod/rq4_fit_sequences_<tr>_full.rds
#   path_mod/rq4_seq_levels_<tr>_full.rds
#   path_mod/rq4_pid_levels_<tr>_full.rds
#   path_src/master_sequences.csv
#
# OUTPUTS (per treatment)
#   path_out/ex1_1_<tr>_sequences.csv
#   path_out/ex1_1_<tr>_participants.csv
#   path_mod/ex1_1_<tr>_sequences.rds
#   path_mod/ex1_1_<tr>_participants.rds
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex1_1_tables <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  lab_heads  <- as.character(design$seq$side_labels$heads)
  lab_tails  <- as.character(design$seq$side_labels$tails)
  pure_heads <- as.character(design$seq$anchor_labels$pure_heads)
  pure_tails <- as.character(design$seq$anchor_labels$pure_tails)
  
  eps_vec    <- as.numeric(design$ex1$eps)
  eps_main   <- eps_vec[1]
  eta        <- as.numeric(design$ex1$eta)
  delta_vec  <- as.numeric(design$ex1$delta)
  delta_main <- delta_vec[1]
  
  stopifnot(
    length(eps_vec) >= 1L, all(is.finite(eps_vec)), all(eps_vec > 0),
    length(eta) == 1L, is.finite(eta), eta > 0,
    length(delta_vec) >= 1L, all(is.finite(delta_vec)), all(delta_vec > 0)
  )
  
  # ---- Load trial data (once) ----
  f_master <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  dt <- fread(f_master, encoding = "UTF-8")
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  # Betting trials only; restrict to H/T side choices
  dt <- dt[is.finite(stake) & stake > 0 & side %in% c(lab_heads, lab_tails)]

  dt[, h       := as.integer(side == lab_heads)]
  dt[, sign_is := 2L * h - 1L]
  
  stopifnot(all(dt[side == lab_heads, sign_is] == 1L))
  stopifnot(all(dt[side == lab_tails, sign_is] == -1L))
  
  # ---- Label helpers (preregistered cutpoints) ----
  label_seq <- function(H, T) {
    if (max(H, T) < 0.50) return("neutral")
    if (H >= T) {
      if (H >= 0.95) return("strong_Hbar")
      if (H >= 0.80) return("moderate_Hbar")
      return("weak_Hbar")
    } else {
      if (T >= 0.95) return("strong_Tbar")
      if (T >= 0.80) return("moderate_Tbar")
      return("weak_Tbar")
    }
  }
  
  label_pid <- function(HH, G) {
    if (max(HH, G) < 0.75) return("neutral")
    if (HH >= G) {
      if (HH >= 0.95) return("solid_hot")
      if (HH >= 0.90) return("likely_hot")
      return("leaning_hot")
    } else {
      if (G >= 0.95) return("solid_gambler")
      if (G >= 0.90) return("likely_gambler")
      return("leaning_gambler")
    }
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ---- Load RQ4 posteriors ----
    f_fit <- file.path(path_mod, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_seq <- file.path(path_mod, paste0("rq4_seq_levels_",    tr, "_full.rds"))
    f_pid <- file.path(path_mod, paste0("rq4_pid_levels_",    tr, "_full.rds"))
    
    for (f in c(f_fit, f_seq, f_pid)) {
      if (!file.exists(f)) stop("EX1.1 tables: missing artifact: ", f)
    }
    
    fit        <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    
    post <- rstan::extract(fit)
    
    if (is.null(post$mu_h) || is.null(post$hbar)) {
      stop("EX1.1 tables: RQ4 fit missing 'mu_h' or 'hbar' for tr='", tr,
           "'.\nRegenerate RQ4 fits with the updated rq4_side.stan.")
    }
    
    thetaS <- post$mu_h           # K x S
    theta0 <- as.numeric(post$hbar)  # K
    K <- nrow(thetaS)
    S <- length(seq_levels)
    
    stopifnot(ncol(thetaS) == S, length(theta0) == K)
    
    idx_H <- match(pure_heads, seq_levels)
    idx_T <- match(pure_tails, seq_levels)
    if (is.na(idx_H)) stop("EX1.1 tables: anchor '", pure_heads, "' not in seq_levels for tr='", tr, "'.")
    if (is.na(idx_T)) stop("EX1.1 tables: anchor '", pure_tails, "' not in seq_levels for tr='", tr, "'.")
    
    thetaH <- thetaS[, idx_H]  # K
    thetaT <- thetaS[, idx_T]  # K
    
  
    
    # ---- Trial data for this treatment ----
    d <- dt[treat == tr]
    if (nrow(d) == 0L) {
      warning("EX1.1 tables: no betting trials for tr='", tr, "'. Skipping.")
      next
    }

    
    d[, sid := match(seq, seq_levels)]
    d <- d[!is.na(sid)]
    if (nrow(d) == 0L) next
    
    pid_levels <- sort(unique(d$pid))
    d[, pid_i := match(pid, pid_levels)]
    Np <- length(pid_levels)
    
    n_by_sid <- d[, .(n = .N), by = sid]
    n_sid <- integer(S)
    n_sid[n_by_sid$sid] <- n_by_sid$n
    
    n_by_pid <- d[, .(n = .N), by = pid_i]
    n_pid <- integer(Np)
    n_pid[n_by_pid$pid_i] <- n_by_pid$n
    
    sid_vec  <- as.integer(d$sid)
    pid_vec  <- as.integer(d$pid_i)
    sign_vec <- as.numeric(d$sign_is)
    
    # ---- Base tables ----
    seq_tbl  <- data.table(sequence = seq_levels, n_bets = as.integer(n_sid))
    part_tbl <- data.table(pid = pid_levels,      n_bets = as.integer(n_pid))
    
    chi_s_draws_main <- NULL
    chi_i_draws_main <- NULL
    
    # ---- Main loop: per epsilon ----
    for (eps in eps_vec) {
      
      suf_eps <- gsub("\\.", "", sprintf("%.2f", eps))
      
      # Draw-level indicator matrices (K x S)
      I_H <- abs(thetaS - matrix(thetaH, nrow = K, ncol = S)) < eps
      I_T <- abs(thetaS - matrix(thetaT, nrow = K, ncol = S)) < eps
      I_0 <- abs(thetaS - matrix(theta0,  nrow = K, ncol = S)) < eps
      
     
      
      # Collapsed scalar similarities (reported in seq_tbl)
      dH <- colMeans(I_H)
      dT <- colMeans(I_T)
      d0 <- colMeans(I_0)
      denom_sc <- dH + dT + d0 + 3 * eta
      
      seq_tbl[, (paste0("dH_eps_", suf_eps)) := dH]
      seq_tbl[, (paste0("dT_eps_", suf_eps)) := dT]
      seq_tbl[, (paste0("d0_eps_", suf_eps)) := d0]
      seq_tbl[, (paste0("wH_eps_", suf_eps)) := (dH + eta) / denom_sc]
      seq_tbl[, (paste0("wT_eps_", suf_eps)) := (dT + eta) / denom_sc]
      seq_tbl[, (paste0("w0_eps_", suf_eps)) := (d0 + eta) / denom_sc]
      
      # Draw-level weights (K x S) and weight difference
      denom  <- I_H + I_T + I_0 + 3 * eta
      wH     <- (I_H + eta) / denom
      wT     <- (I_T + eta) / denom
      diff_w <- wH - wT   # K x S
      
      # Draw loop: chi_s^(k) and chi_i^(k) via rowsum
      chi_s_draws <- matrix(NA_real_, nrow = K, ncol = S)
      chi_i_draws <- matrix(NA_real_, nrow = K, ncol = Np)
      
      for (k in seq_len(K)) {
        
        z_k <- sign_vec * diff_w[k, sid_vec]
        
        sum_s   <- rowsum(z_k, group = sid_vec, reorder = FALSE)
        sid_u   <- as.integer(rownames(sum_s))
        chi_s_draws[k, sid_u] <- as.numeric(sum_s[, 1]) / n_sid[sid_u]
        
        sum_i   <- rowsum(z_k, group = pid_vec, reorder = FALSE)
        pid_u   <- as.integer(rownames(sum_i))
        chi_i_draws[k, pid_u] <- as.numeric(sum_i[, 1]) / n_pid[pid_u]
        
        
      
        
      }
      
      # Sequence summaries
      seq_tbl[, (paste0("chi_median_eps_", suf_eps)) :=
                apply(chi_s_draws, 2, median,   na.rm = TRUE)]
      seq_tbl[, (paste0("chi_mean_eps_",   suf_eps)) :=
                apply(chi_s_draws, 2, mean,     na.rm = TRUE)]
      seq_tbl[, (paste0("chi_q025_eps_", suf_eps)) :=
                as.numeric(apply(chi_s_draws, 2, quantile, probs = 0.025, na.rm = TRUE))]
      seq_tbl[, (paste0("chi_q975_eps_", suf_eps)) :=
                as.numeric(apply(chi_s_draws, 2, quantile, probs = 0.975, na.rm = TRUE))]
      
      for (delta in delta_vec) {
        suf_d <- gsub("\\.", "", sprintf("%.2f", delta))
        seq_tbl[, (paste0("H_eps_",    suf_eps, "_delta_", suf_d)) :=
                  apply(chi_s_draws, 2, function(x) mean(x >  delta, na.rm = TRUE))]
        seq_tbl[, (paste0("T_eps_",    suf_eps, "_delta_", suf_d)) :=
                  apply(chi_s_draws, 2, function(x) mean(x < -delta, na.rm = TRUE))]
        seq_tbl[, (paste0("Pabs_eps_", suf_eps, "_delta_", suf_d)) :=
                  apply(chi_s_draws, 2, function(x) mean(abs(x) > delta, na.rm = TRUE))]
      }
      
      # Participant summaries at main eps only
      if (isTRUE(all.equal(eps, eps_main))) {
        
        part_tbl[, chi_median := apply(chi_i_draws, 2, median,   na.rm = TRUE)]
        part_tbl[, chi_mean   := apply(chi_i_draws, 2, mean,     na.rm = TRUE)]
        part_tbl[, chi_q025 := as.numeric(apply(chi_i_draws, 2, quantile, probs = 0.025, na.rm = TRUE))]
        part_tbl[, chi_q975 := as.numeric(apply(chi_i_draws, 2, quantile, probs = 0.975, na.rm = TRUE))]
        
        for (delta in delta_vec) {
          suf_d <- gsub("\\.", "", sprintf("%.2f", delta))
          part_tbl[, (paste0("HH_delta_",   suf_d)) :=
                     apply(chi_i_draws, 2, function(x) mean(x >  delta, na.rm = TRUE))]
          part_tbl[, (paste0("G_delta_",    suf_d)) :=
                     apply(chi_i_draws, 2, function(x) mean(x < -delta, na.rm = TRUE))]
          part_tbl[, (paste0("Pabs_delta_", suf_d)) :=
                     apply(chi_i_draws, 2, function(x) mean(abs(x) > delta, na.rm = TRUE))]
        }
        
        suf_dm <- gsub("\\.", "", sprintf("%.2f", delta_main))
        part_tbl[, class := mapply(label_pid,
                                   HH = .SD[[paste0("HH_delta_", suf_dm)]],
                                   G  = .SD[[paste0("G_delta_",  suf_dm)]])]
        
        chi_s_draws_main <- chi_s_draws
        chi_i_draws_main <- chi_i_draws
      }
    }
    
    # Sequence direction label (main eps + main delta)
    suf_em <- gsub("\\.", "", sprintf("%.2f", eps_main))
    suf_dm <- gsub("\\.", "", sprintf("%.2f", delta_main))
    seq_tbl[, direction_label := mapply(label_seq,
                                        H = .SD[[paste0("H_eps_", suf_em, "_delta_", suf_dm)]],
                                        T = .SD[[paste0("T_eps_", suf_em, "_delta_", suf_dm)]])]
    
    # ---- Write outputs ----
    f_seq_csv <- file.path(path_out, paste0("ex1_1_", tr, "_sequences.csv"))
    f_pid_csv <- file.path(path_out, paste0("ex1_1_", tr, "_participants.csv"))
    f_seq_rds <- file.path(path_mod, paste0("ex1_1_", tr, "_sequences.rds"))
    f_pid_rds <- file.path(path_mod, paste0("ex1_1_", tr, "_participants.rds"))
    
    if (!should_skip(c(f_seq_csv, f_pid_csv), cfg, "output",
                     paste0("EX1.1 tables (", tr, ")"))) {
      fwrite(seq_tbl,  f_seq_csv)
      fwrite(part_tbl, f_pid_csv)
      msg("Saved: ", f_seq_csv)
      msg("Saved: ", f_pid_csv)
    }
    
    if (!should_skip(c(f_seq_rds, f_pid_rds), cfg, "model",
                     paste0("EX1.1 RDS (", tr, ")"))) {
      saveRDS(
        list(dataset    = cfg$run$data_folder,
             treatment  = tr,
             seq_levels = seq_levels,
             eps_main   = eps_main,
             delta_main = delta_main,
             chi_draws  = chi_s_draws_main),
        f_seq_rds
      )
      saveRDS(
        list(dataset    = cfg$run$data_folder,
             treatment  = tr,
             pid_levels = pid_levels,
             eps_main   = eps_main,
             delta_main = delta_main,
             chi_draws  = chi_i_draws_main),
        f_pid_rds
      )
      msg("Saved: ", f_seq_rds)
      msg("Saved: ", f_pid_rds)
    }
    
    outputs[[tr]] <- list(sequences = seq_tbl, participants = part_tbl)
  }
  
  invisible(outputs)
}