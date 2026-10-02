# ============================================================
# 43_rq4_tables.R
#   RQ4 tables: sequence-level and participant-level posterior
#   summaries for side choice (Heads vs Tails conditional on betting)
#
# PURPOSE
#   Produces sequence-level and participant-level posterior summary
#   tables, a model summary table, and a sequences summary table.
#   Mirrors the structure of 13_rq1_tables.R.
#
# ESTIMANDS
#   Sequence level:
#     mu_h[s]  = E_{u~N(0,sigma_u)} logistic(alpha + u + beta_s)
#                population mean Heads probability for sequence s
#     hbar     = E_{u~N(0,sigma_u)} logistic(alpha + u)
#                population baseline (no sequence effect)
#     H_delta[s] = P(mu_h[s] > hbar + delta)  Head-biased
#     T_delta[s] = P(mu_h[s] < hbar - delta)  Tail-biased
#   Participant level:
#     mu_h_i[i] = logistic(alpha + u_i)
#                 participant baseline (= sequence-averaged tendency
#                 under sum-to-zero constraint on beta_s; from Stan GQ)
#     H_delta[i] = P(mu_h_i[i] > hbar + delta)
#     T_delta[i] = P(mu_h_i[i] < hbar - delta)
#
# DELTA GRID
#   cfg$design$rq4$delta = c(0.05, 0.03, 0.08)  first = main
#
# NO TAG LOOP
#   Single "full" tag per treatment (see 41_ preamble for rationale).
#
# CONFIG KEYS USED
#   cfg$run$treatment
#   cfg$run$overwrite_outputs
#   cfg$design$seq$side_labels
#   cfg$design$rq4$delta
#
# INPUTS
#   path_src/master_sequences.csv
#   path_mod/rq4_fit_sequences_<tr>_full.rds
#   path_mod/rq4_pid_levels_<tr>_full.rds
#   path_mod/rq4_seq_levels_<tr>_full.rds
#
# OUTPUTS
#   path_out/rq4_<tr>_sequences.csv
#   path_out/rq4_<tr>_sequences_summary.csv
#   path_out/rq4_<tr>_participants.csv
#   path_out/rq4_<tr>_model_summary.csv
#
# CALL ORDER
#   rq4_stan(cfg)        -- 41_
#   rq4_diagnostics(cfg) -- 42_
#   rq4_tables(cfg)      -- this script
#   rq4_figures(cfg)     -- 44_
# ============================================================

library(data.table)

rq4_tables <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  stopifnot(!is.null(design$seq$side_labels))
  
  # Delta grid: main (first) + sensitivity values
  delta_vec  <- as.numeric(design$rq4$delta %||% c(0.05, 0.03, 0.08))
  stopifnot(length(delta_vec) >= 1L, all(is.finite(delta_vec)), all(delta_vec > 0))
  delta_main <- delta_vec[1]
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(nzchar(lab_heads), nzchar(lab_tails))
  
  # ----------------------------
  # Paths
  # ----------------------------
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  
  master <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "stake", "side") %in% names(master)))
  
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[, side  := as.character(side)]
  master[is.na(stake), stake := 0]
  
  # ----------------------------
  # Direction labelling helpers
  # Mirrors RQ1 label_seq / label_pid but two-sided (Head or Tail bias)
  # ----------------------------
  
  # Sequence: H/T = P(mu_h[s] >/< hbar +/- delta_main)
  label_seq <- function(H, T) {
    if (max(H, T) < 0.50) return("neutral")
    if (H >= T) {
      if (H >= 0.95) "strong_head"
      else if (H >= 0.80) "moderate_head"
      else "weak_head"
    } else {
      if (T >= 0.95) "strong_tail"
      else if (T >= 0.80) "moderate_tail"
      else "weak_tail"
    }
  }
  
  # Participant: H/T = P(mu_h_i[i] >/< hbar +/- delta_main)
  label_pid <- function(H, T) {
    if (max(H, T) < 0.75) return("neutral")
    if (H >= T) {
      if (H >= 0.95) "solid_head"
      else if (H >= 0.90) "likely_head"
      else "leaning_head"
    } else {
      if (T >= 0.95) "solid_tail"
      else if (T >= 0.90) "likely_tail"
      else "leaning_tail"
    }
  }
  
  # Helper for model summary rows
  sum_draw <- function(x, nm) {
    data.table(
      parameter = nm,
      median    = median(x),
      mean      = mean(x),
      q025      = as.numeric(quantile(x, 0.025)),
      q975      = as.numeric(quantile(x, 0.975))
    )
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(path_mod, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_pid <- file.path(path_mod, paste0("rq4_pid_levels_",    tr, "_full.rds"))
    f_seq <- file.path(path_mod, paste0("rq4_seq_levels_",    tr, "_full.rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) {
      warning("RQ4 tables: missing artifacts for tr='", tr, "'. Skipping.")
      next
    }
    
    fit        <- readRDS(f_fit)
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    N          <- length(pid_levels)
    S          <- length(seq_levels)
    stopifnot(N >= 1L, S >= 2L)
    
    post <- rstan::extract(fit)
    
    # Validate required generated quantities
    for (gq in c("mu_h", "hbar", "mu_h_i", "alpha", "sigma_u", "sigma_s",
                 "gamma_drift")) {
      if (is.null(post[[gq]])) {
        stop("RQ4 tables: missing '", gq, "' in fit for tr='", tr,
             "'. Refit with updated rq4.stan.")
      }
    }
    
    mu_s_draws  <- post$mu_h      # iters x S  population mean per sequence
    hbar_draws  <- post$hbar      # iters      population baseline
    mu_i_draws  <- post$mu_h_i    # iters x N  participant baselines (from Stan GQ)
    iters       <- length(hbar_draws)
    
    stopifnot(is.matrix(mu_s_draws),  nrow(mu_s_draws) == iters, ncol(mu_s_draws) == S)
    stopifnot(is.matrix(mu_i_draws),  nrow(mu_i_draws) == iters, ncol(mu_i_draws) == N)
    
    # Betting trial counts (for n_trials column)
    d_counts <- master[
      treat == tr &
        pid  %in% pid_levels &
        seq  %in% seq_levels &
        is.finite(stake) & stake > 0 &
        side %in% c(lab_heads, lab_tails)
    ]
    n_by_seq <- d_counts[, .(n_trials = .N), by = seq]
    setkey(n_by_seq, seq)
    n_by_pid <- d_counts[, .(n_trials = .N), by = pid]
    setkey(n_by_pid, pid)
    
    # ========================================================
    # SEQUENCES
    # ========================================================
    seq_tbl <- data.table(
      sequence    = seq_levels,
      n_trials    = as.integer(n_by_seq[.(seq_levels), n_trials]),
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
    
    H_main_s <- numeric(S)
    T_main_s <- numeric(S)
    
    for (dlt in delta_vec) {
      nm  <- gsub("\\.", "", sprintf("%.2f", dlt))
      nmH <- paste0("H_delta_", nm)
      nmT <- paste0("T_delta_", nm)
      nmA <- paste0("P_absdev_delta_", nm)
      
      H <- apply(mu_s_draws, 2, function(x) mean(x > hbar_draws + dlt))
      T <- apply(mu_s_draws, 2, function(x) mean(x < hbar_draws - dlt))
      A <- apply(mu_s_draws, 2, function(x) mean(abs(x - 0.5) > dlt))
      
      seq_tbl[, (nmH) := H]
      seq_tbl[, (nmT) := T]
      seq_tbl[, (nmA) := A]
      
      if (isTRUE(all.equal(dlt, delta_main))) {
        H_main_s <- H
        T_main_s <- T
      }
    }
    
    seq_tbl[, direction_label := mapply(label_seq, H_main_s, T_main_s)]
    setorder(seq_tbl, sequence)
    
    f_seq_csv <- file.path(path_out, paste0("rq4_", tr, "_sequences.csv"))
    if (!should_skip(f_seq_csv, cfg, "output",
                     paste0("RQ4 sequences (", tr, ")"))) {
      fwrite(seq_tbl, f_seq_csv)
      msg("Saved: ", f_seq_csv)
    }
    
    # ---- Sequences summary ----
    f_seq_sum <- file.path(path_out, paste0("rq4_", tr, "_sequences_summary.csv"))
    if (!should_skip(f_seq_sum, cfg, "output",
                     paste0("RQ4 sequences summary (", tr, ")"))) {
      
      direction_levels <- c("strong_head", "moderate_head", "weak_head",
                            "neutral",
                            "weak_tail", "moderate_tail", "strong_tail")
      
      cnt <- seq_tbl[, .N, by = direction_label]
      cnt[, direction_label := factor(direction_label, levels = direction_levels)]
      setorder(cnt, direction_label)
      cnt[, pct := round(100 * N / nrow(seq_tbl), 1)]
      
      stats <- seq_tbl[, .(
        mu_mean   = round(mean(mu_h_mean),   3),
        mu_median = round(median(mu_h_mean), 3),
        mu_min    = round(min(mu_h_mean),    3),
        mu_max    = round(max(mu_h_mean),    3)
      ), by = direction_label]
      
      ids <- seq_tbl[, .(
        sequences = paste(sort(sequence), collapse = ", ")
      ), by = direction_label]
      
      seq_summary <- Reduce(function(a, b) merge(a, b, by = "direction_label"),
                            list(cnt, stats, ids))
      fwrite(seq_summary, f_seq_sum)
      msg("Saved: ", f_seq_sum)
    }
    
    # ========================================================
    # PARTICIPANTS
    # Uses mu_h_i[N] from Stan generated quantities directly.
    # mu_h_i[i] = logistic(alpha + u_i) = sequence-averaged
    # participant tendency under sum-to-zero constraint on beta_s.
    # ========================================================
    part_tbl <- data.table(
      pid         = pid_levels,
      n_trials    = as.integer(n_by_pid[.(pid_levels), n_trials]),
      mu_h_median = apply(mu_i_draws, 2, median),
      mu_h_mean   = apply(mu_i_draws, 2, mean),
      mu_h_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
      mu_h_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975),
      hbar_median = median(hbar_draws),
      hbar_mean   = mean(hbar_draws),
      hbar_q025   = as.numeric(quantile(hbar_draws, 0.025)),
      hbar_q975   = as.numeric(quantile(hbar_draws, 0.975))
    )
    part_tbl[is.na(n_trials), n_trials := 0L]
    
    H_main_i <- numeric(N)
    T_main_i <- numeric(N)
    
    for (dlt in delta_vec) {
      nm  <- gsub("\\.", "", sprintf("%.2f", dlt))
      nmH <- paste0("H_delta_", nm)
      nmT <- paste0("T_delta_", nm)
      nmA <- paste0("P_absdev_delta_", nm)
      
      H <- apply(mu_i_draws, 2, function(x) mean(x > hbar_draws + dlt))
      T <- apply(mu_i_draws, 2, function(x) mean(x < hbar_draws - dlt))
      A <- apply(mu_i_draws, 2, function(x) mean(abs(x - 0.5) > dlt))
      
      part_tbl[, (nmH) := H]
      part_tbl[, (nmT) := T]
      part_tbl[, (nmA) := A]
      
      if (isTRUE(all.equal(dlt, delta_main))) {
        H_main_i <- H
        T_main_i <- T
      }
    }
    
    part_tbl[, side_label := mapply(label_pid, H_main_i, T_main_i)]
    setorder(part_tbl, pid)
    
    f_pid_csv <- file.path(path_out, paste0("rq4_", tr, "_participants.csv"))
    if (!should_skip(f_pid_csv, cfg, "output",
                     paste0("RQ4 participants (", tr, ")"))) {
      fwrite(part_tbl, f_pid_csv)
      msg("Saved: ", f_pid_csv)
    }
    
    # ---- Participants summary ----
    f_pid_sum <- file.path(path_out, paste0("rq4_", tr, "_participants_summary.csv"))
    if (!should_skip(f_pid_sum, cfg, "output",
                     paste0("RQ4 participants summary (", tr, ")"))) {
      
      side_levels <- c("solid_head", "likely_head", "leaning_head",
                       "neutral",
                       "leaning_tail", "likely_tail", "solid_tail")
      
      cnt <- part_tbl[, .N, by = side_label]
      cnt[, side_label := factor(side_label, levels = side_levels)]
      setorder(cnt, side_label)
      cnt[, pct := round(100 * N / nrow(part_tbl), 1)]
      
      stats <- part_tbl[, .(
        mu_mean   = round(mean(mu_h_mean),   3),
        mu_median = round(median(mu_h_mean), 3),
        mu_min    = round(min(mu_h_mean),    3),
        mu_max    = round(max(mu_h_mean),    3)
      ), by = side_label]
      
      ids <- part_tbl[, .(
        pids = paste(sort(pid), collapse = ", ")
      ), by = side_label]
      
      pid_summary <- Reduce(function(a, b) merge(a, b, by = "side_label"),
                            list(cnt, stats, ids))
      
      fwrite(pid_summary, f_pid_sum)
      msg("Saved: ", f_pid_sum)
    }
    
    # ========================================================
    # MODEL SUMMARY
    # ========================================================
    f_mod_csv <- file.path(path_out, paste0("rq4_", tr, "_model_summary.csv"))
    if (!should_skip(f_mod_csv, cfg, "output",
                     paste0("RQ4 model summary (", tr, ")"))) {
      
      # drift_type in data tells us whether gamma_drift was active
      # Log it from the diagnostics CSV if available, else from fit data
      f_diag <- file.path(path_out, "rq4_diagnostics.csv")
      drift_type_str <- if (file.exists(f_diag)) {
        d_diag <- fread(f_diag)
        row    <- d_diag[treatment == tr]
        if (nrow(row) == 1L) as.character(row$drift) else "unknown"
      } else "unknown"
      
      mod_rows <- list(
        sum_draw(post$alpha,    "alpha (log-odds intercept)"),
        sum_draw(hbar_draws,    "hbar (population baseline Heads prob)"),
        sum_draw(post$sigma_u,  "sigma_u (between-participant SD)"),
        sum_draw(post$sigma_s,  "sigma_s (between-sequence SD)")
      )
      
      # Include drift parameter row only when drift was active
      if (!is.null(post$gamma_drift)) {
        mod_rows <- c(mod_rows,
                      list(sum_draw(post$gamma_drift,
                                    "gamma_drift (linear drift per block)")))
      }
      
      mod_tbl <- rbindlist(mod_rows)
      mod_tbl[, treatment    := tr]
      mod_tbl[, drift        := drift_type_str]
      mod_tbl[, n_trials     := nrow(d_counts)]
      mod_tbl[, n_participants := N]
      mod_tbl[, n_sequences  := S]
      
      fwrite(mod_tbl, f_mod_csv)
      msg("Saved: ", f_mod_csv)
    }
    
    outputs[[tr]] <- list(
      sequences_csv        = f_seq_csv,
      sequences_summary_csv = f_seq_sum,
      participants_csv     = f_pid_csv,
      model_summary_csv    = f_mod_csv,
      participants_summary_csv = f_pid_sum
    )
  }
  
  invisible(outputs)
}