# ============================================================
# scripts/analysis/53_ex1_1_tables.R
#
# EX1.I tables (single pass): sequences + participants
# (draw-level chi).
#
# Reads RQ4 posterior draws and constructs EX1.1 quantities:
#
# 1) From RQ4 (per treatment):
#    - theta_s^(k) = mu_h[k, s]  (population Heads-choice prob for sequence s)
#    - theta0^(k)  = hbar[k]     (population baseline from RQ4 generated quantities)
#    - anchors: pure_heads, pure_tails from cfg$design$seq$anchor_labels
#
# 2) Per epsilon (eps):
#    - draw-level similarity indicators:
#        I_H(k,s) = 1{|theta_s^(k) - theta_H^(k)| < eps}
#        I_T(k,s) = 1{|theta_s^(k) - theta_T^(k)| < eps}
#        I_0(k,s) = 1{|theta_s^(k) - theta0^(k) | < eps}
#    - draw-level weights (normalized per k,s):
#        w_a(k,s) = (I_a(k,s) + eta) / (I_H + I_T + I_0 + 3*eta)
#      and diff_w(k,s) = wH(k,s) - wT(k,s)
#
# 3) For each betting trial (i,s) with observed h_is in {0,1}:
#    - sign_is = 2*h_is - 1  (Heads=+1, Tails=-1)
#    - z_is^(k) = sign_is * diff_w(k, s)
#
# 4) Directional tendencies (posterior draws):
#    - chi_s^(k) = mean_{trials on s} z_is^(k)
#    - chi_i^(k) = mean_{trials of i} z_is^(k)
#
# 5) Report summaries + prereg probability masses and labels:
#    - sequences: H(.)=P(chi_s >  delta), T(.)=P(chi_s < -delta), Pabs(.)=P(|chi_s|>delta)
#      with delta grid from cfg$design$rq4$delta (defaults to 0.05,0.03,0.08)
#    - participants: HH(.)=P(chi_i >  delta), G(.)=P(chi_i < -delta), Pabs(.)=P(|chi_i|>delta)
#      and prereg cutpoint classes.
#
# Also reports (scalar) posterior similarity probabilities and weights:
#    d_a(s) = mean_k I_a(k,s), w_a(s) = normalize(d_a(s)+eta)
#
# Outputs per treatment:
#   data/clean/<ds>/output/ex1_1_<tr>_sequences.csv
#   data/clean/<ds>/output/ex1_1_<tr>_participants.csv
#   data/clean/<ds>/models/ex1_1_<tr>_sequences.rds
#   data/clean/<ds>/models/ex1_1_<tr>_participants.rds
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex1_1_tables <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  # ---- labels ----
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  
  pure_heads <- as.character(design$seq$anchor_labels$pure_heads)  # e.g., "HHHHHH"
  pure_tails <- as.character(design$seq$anchor_labels$pure_tails)  # e.g., "OOOOOO"
  
  # ---- eps grid (prereg) ----
  eps_vec  <- c(0.05, 0.03, 0.08)
  eps_main <- 0.05
  
  # ---- delta grid (prereg) ----
  delta_vec <- NULL
  if (!is.null(design$rq4) && !is.null(design$rq4$delta)) delta_vec <- design$rq4$delta
  if (is.null(delta_vec)) delta_vec <- c(0.05, 0.03, 0.08)
  delta_vec <- as.numeric(delta_vec)
  delta_main <- delta_vec[1]
  
  # ---- stabilizer (prereg) ----
  eta <- 1e-6
  
  # ---- paths ----
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  dt <- fread(infile, encoding = "UTF-8")
  
  # normalize schema (pipeline already ensures columns exist)
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  # betting trials only; keep only H/T
  dt <- dt[is.finite(stake) & stake > 0 & side %in% c(lab_heads, lab_tails)]
  dt[, h := as.integer(side == lab_heads)]    # 1=Heads, 0=Tails
  dt[, sign := 2L * h - 1L]                   # +1 Heads, -1 Tails
  
  # ---- label helpers (PREREG CUTPOINTS) ----
  label_seq <- function(H, T) {
    if (max(H, T) < 0.50) return("neutral")
    if (H >= T) {
      if (H >= 0.95) return("strong_Hbar")
      if (H >= 0.80) return("moderate_Hbar")
      return("weak_Hbar")   # H >= 0.50
    } else {
      if (T >= 0.95) return("strong_Tbar")
      if (T >= 0.80) return("moderate_Tbar")
      return("weak_Tbar")   # T >= 0.50
    }
  }
  
  label_pid <- function(HH, G) {
    if (max(HH, G) < 0.75) return("neutral")
    if (HH >= G) {
      if (HH >= 0.95) return("solid_hot")
      if (HH >= 0.90) return("likely_hot")
      return("leaning_hot")  # HH >= 0.75
    } else {
      if (G >= 0.95) return("solid_gambler")
      if (G >= 0.90) return("likely_gambler")
      return("leaning_gambler")  # G >= 0.75
    }
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ----------------------------
    # Load RQ4 fit + levels (per treatment)
    # ----------------------------
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, "_full.rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, "_full.rds"))
    
    fit <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    pid_levels_fit <- as.character(readRDS(f_pid))  # only for alignment diagnostics if needed
    
    post <- rstan::extract(fit)
    
    # Required from updated rq4_side.stan:
    # - mu_h: K x S
    # - hbar: K
    mu_h  <- post$mu_h
    hbar  <- as.numeric(post$hbar)
    
    K <- nrow(mu_h)
    S <- ncol(mu_h)
    
    # Anchors must exist in seq_levels
    idx_H <- match(pure_heads, seq_levels)
    idx_T <- match(pure_tails, seq_levels)
    
    thetaS <- mu_h                    # K x S
    thetaH <- thetaS[, idx_H]         # K
    thetaT <- thetaS[, idx_T]         # K
    theta0 <- hbar                    # K
    
    # ----------------------------
    # Observed betting trials for this treatment
    # ----------------------------
    d <- dt[treat == tr]
    if (nrow(d) == 0) next
    
    d[, sid := match(seq, seq_levels)]
    # drop sequences not in levels (should be none if pipeline consistent)
    d <- d[!is.na(sid)]
    if (nrow(d) == 0) next
    
    # participant index (within observed betting sample)
    pid_levels <- sort(unique(d$pid))
    d[, pid_i := match(pid, pid_levels)]
    Np <- length(pid_levels)
    
    # counts needed for means
    n_by_sid <- d[, .(n = .N), by = sid]
    n_sid <- integer(S); n_sid[] <- 0L
    n_sid[n_by_sid$sid] <- n_by_sid$n
    
    n_by_pid <- d[, .(n = .N), by = pid_i]
    n_pid <- integer(Np); n_pid[] <- 0L
    n_pid[n_by_pid$pid_i] <- n_by_pid$n
    
    # ----------------------------
    # Base output tables
    # ----------------------------
    seq_tbl  <- data.table(sequence = seq_levels, n_bets = as.integer(n_sid))
    part_tbl <- data.table(pid = pid_levels, n_bets = as.integer(n_pid))
    
    # placeholders for saved draws at main eps
    chi_s_draws_main <- NULL
    chi_i_draws_main <- NULL
    pid_levels_main  <- NULL
    
    # ----------------------------
    # For each eps: compute draw-level weights and chi draws
    # ----------------------------
    for (eps in eps_vec) {
      
      suf_eps <- gsub("\\.", "", sprintf("%.2f", eps))
      
      # draw-level similarity indicators (K x S)
      I_H <- abs(thetaS - matrix(thetaH, nrow = K, ncol = S)) < eps
      I_T <- abs(thetaS - matrix(thetaT, nrow = K, ncol = S)) < eps
      I_0 <- abs(thetaS - matrix(theta0, nrow = K, ncol = S)) < eps
      
      # scalar similarity probabilities (reported)
      dH <- colMeans(I_H)
      dT <- colMeans(I_T)
      d0 <- colMeans(I_0)
      
      denom_scalar <- dH + dT + d0 + 3 * eta
      wH_scalar <- (dH + eta) / denom_scalar
      wT_scalar <- (dT + eta) / denom_scalar
      w0_scalar <- (d0 + eta) / denom_scalar
      
      seq_tbl[, (paste0("dH_eps_", suf_eps)) := dH]
      seq_tbl[, (paste0("dT_eps_", suf_eps)) := dT]
      seq_tbl[, (paste0("d0_eps_", suf_eps)) := d0]
      seq_tbl[, (paste0("wH_eps_", suf_eps)) := wH_scalar]
      seq_tbl[, (paste0("wT_eps_", suf_eps)) := wT_scalar]
      seq_tbl[, (paste0("w0_eps_", suf_eps)) := w0_scalar]
      
      # draw-level weights (K x S)
      denom <- (I_H + I_T + I_0) + 3 * eta
      wH <- (I_H + eta) / denom
      wT <- (I_T + eta) / denom
      diff_w <- wH - wT   # K x S
      
      # ---- compute chi_s^(k) and chi_i^(k) efficiently via rowsum per draw ----
      chi_s_draws <- matrix(NA_real_, nrow = K, ncol = S)
      chi_i_draws <- matrix(NA_real_, nrow = K, ncol = Np)
      
      sid_vec  <- as.integer(d$sid)
      pid_vec  <- as.integer(d$pid_i)
      sign_vec <- as.numeric(d$sign)
      
      for (k in seq_len(K)) {
        
        z_k <- sign_vec * diff_w[k, sid_vec]
        
        # sequences
        sum_by_sid <- rowsum(z_k, group = sid_vec, reorder = FALSE)
        sid_uniq <- as.integer(rownames(sum_by_sid))
        chi_s_draws[k, sid_uniq] <- as.numeric(sum_by_sid[, 1]) / n_sid[sid_uniq]
        
        # participants
        sum_by_pid <- rowsum(z_k, group = pid_vec, reorder = FALSE)
        pid_uniq <- as.integer(rownames(sum_by_pid))
        chi_i_draws[k, pid_uniq] <- as.numeric(sum_by_pid[, 1]) / n_pid[pid_uniq]
      }
      
      # ---- sequence summaries for this eps ----
      seq_tbl[, (paste0("chi_median_eps_", suf_eps)) := apply(chi_s_draws, 2, median, na.rm = TRUE)]
      seq_tbl[, (paste0("chi_mean_eps_",   suf_eps)) := apply(chi_s_draws, 2, mean,   na.rm = TRUE)]
      seq_tbl[, (paste0("chi_q025_eps_",   suf_eps)) := apply(chi_s_draws, 2, quantile, probs = 0.025, na.rm = TRUE)]
      seq_tbl[, (paste0("chi_q975_eps_",   suf_eps)) := apply(chi_s_draws, 2, quantile, probs = 0.975, na.rm = TRUE)]
      
      for (delta in delta_vec) {
        suf_d <- gsub("\\.", "", sprintf("%.2f", delta))
        seq_tbl[, (paste0("H_eps_", suf_eps, "_delta_", suf_d)) :=
                  apply(chi_s_draws, 2, function(x) mean(x >  delta, na.rm = TRUE))]
        seq_tbl[, (paste0("T_eps_", suf_eps, "_delta_", suf_d)) :=
                  apply(chi_s_draws, 2, function(x) mean(x < -delta, na.rm = TRUE))]
        seq_tbl[, (paste0("Pabs_eps_", suf_eps, "_delta_", suf_d)) :=
                  apply(chi_s_draws, 2, function(x) mean(abs(x) > delta, na.rm = TRUE))]
      }
      
      # ---- participant summaries ONLY for MAIN eps (prereg tables focus on main eps) ----
      if (isTRUE(all.equal(eps, eps_main))) {
        part_tbl[, chi_median := apply(chi_i_draws, 2, median, na.rm = TRUE)]
        part_tbl[, chi_mean   := apply(chi_i_draws, 2, mean,   na.rm = TRUE)]
        part_tbl[, chi_q025   := apply(chi_i_draws, 2, quantile, probs = 0.025, na.rm = TRUE)]
        part_tbl[, chi_q975   := apply(chi_i_draws, 2, quantile, probs = 0.975, na.rm = TRUE)]
        
        for (delta in delta_vec) {
          suf_d <- gsub("\\.", "", sprintf("%.2f", delta))
          part_tbl[, (paste0("HH_delta_", suf_d)) :=
                     apply(chi_i_draws, 2, function(x) mean(x >  delta, na.rm = TRUE))]
          part_tbl[, (paste0("G_delta_",  suf_d)) :=
                     apply(chi_i_draws, 2, function(x) mean(x < -delta, na.rm = TRUE))]
          part_tbl[, (paste0("Pabs_delta_", suf_d)) :=
                     apply(chi_i_draws, 2, function(x) mean(abs(x) > delta, na.rm = TRUE))]
        }
        
        # prereg class uses MAIN delta (delta_main) and HH/G cutpoints
        suf_dm <- gsub("\\.", "", sprintf("%.2f", delta_main))
        HHm <- part_tbl[[paste0("HH_delta_", suf_dm)]]
        Gm  <- part_tbl[[paste0("G_delta_",  suf_dm)]]
        part_tbl[, class := mapply(label_pid, HH = HHm, G = Gm)]
        
        chi_s_draws_main <- chi_s_draws
        chi_i_draws_main <- chi_i_draws
        pid_levels_main  <- pid_levels
      }
    }
    
    # ----------------------------
    # Sequence direction_label (MAIN eps + MAIN delta)
    # ----------------------------
    suf_eps_main <- gsub("\\.", "", sprintf("%.2f", eps_main))
    suf_dm <- gsub("\\.", "", sprintf("%.2f", delta_main))
    
    H_main <- seq_tbl[[paste0("H_eps_", suf_eps_main, "_delta_", suf_dm)]]
    T_main <- seq_tbl[[paste0("T_eps_", suf_eps_main, "_delta_", suf_dm)]]
    seq_tbl[, direction_label := mapply(label_seq, H = H_main, T = T_main)]
    
    # ----------------------------
    # Write outputs
    # ----------------------------
    f_seq <- file.path(out_dir, paste0("ex1_1_", tr, "_sequences.csv"))
    f_pid <- file.path(out_dir, paste0("ex1_1_", tr, "_participants.csv"))
    
    fwrite(seq_tbl, f_seq)
    fwrite(part_tbl, f_pid)
    
    msg("Saved: ", f_seq)
    msg("Saved: ", f_pid)
    
    # --- save sequence chi draws for EX4 ---
    seq_rds <- list(
      dataset    = ds,
      treatment  = tr,
      seq_levels = seq_levels,
      eps_main   = eps_main,
      delta_main = delta_main,
      chi_draws  = chi_s_draws_main
    )
    
    f_seq_rds <- file.path(mod_dir, paste0("ex1_1_", tr, "_sequences.rds"))
    saveRDS(seq_rds, f_seq_rds)
    
    msg("Saved: ", f_seq_rds)
    
    # --- save participant chi draws for EX1.2 / EX4 ---
    part_rds <- list(
      dataset    = ds,
      treatment  = tr,
      pid_levels = pid_levels_main,
      eps_main   = eps_main,
      delta_main = delta_main,
      chi_draws  = chi_i_draws_main
    )
    
    f_pid_rds <- file.path(mod_dir, paste0("ex1_1_", tr, "_participants.rds"))
    saveRDS(part_rds, f_pid_rds)
    
    msg("Saved: ", f_pid_rds)
    
    outputs[[tr]] <- list(
      sequences    = seq_tbl,
      participants = part_tbl
    )
  }
  
  invisible(outputs)
}

# Example:
# ex1_1_tables(cfg)