# ============================================================
# scripts/analysis/23_rq2_tables.R
# RQ2 tables (exploratory + confirmatory) from *_full / *_conf fits
# Naming: rq2_<tr>_<type>_<unit>.csv
#   type = exploratory (tag full) | confirmatory (tag conf)
#   unit = sequences | participants
# ============================================================

library(data.table)
library(rstan)

rq2_tables <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  design <- cfg$design
  
  e <- as.numeric(design$seq$endowment)
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  
  rho_vec <- as.numeric(design$rq2$rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  rho_main <- rho_vec[1]
  rho_nm_main <- gsub("\\.", "", sprintf("%.2f", rho_main))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  tags <- c("full", "conf")
  outputs <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      # If you only produce *_conf for betting_normative==TRUE, skip early:
      if (tag == "conf" && !isTRUE(design$a_flags$betting_normative[[tr]])) next
      
      f_fit  <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, "_", tag, ".rds"))
      f_pid  <- file.path(mod_dir, paste0("rq2_pid_levels_",  tr, "_", tag, ".rds"))
      f_seq  <- file.path(mod_dir, paste0("rq2_seq_levels_",  tr, "_", tag, ".rds"))
      f_prep <- file.path(out_dir, paste0("rq2_prepared_", tr, "_", tag, ".csv"))
      
      if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq) || !file.exists(f_prep)) next
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      stopifnot(length(pid_levels) >= 1L, length(seq_levels) >= 1L)
      
      prep <- fread(f_prep, encoding = "UTF-8")
      
      # STRICT schema (prepared is already betting-only + min_bets-filtered + finite z)
      required_cols <- c("pid","seq","delta_bar","sd_star","n_bets")
      if (!all(required_cols %in% names(prep))) {
        stop("Prepared file missing required columns (", paste(setdiff(required_cols, names(prep)), collapse = ", "),
             "). Rerun rq2_stan(cfg) after deleting artifacts.")
      }
      
      prep[, pid := as.character(pid)]
      prep[, seq := as.character(seq)]
      
      # per-pid constants (should be constant within pid)
      pid_map <- prep[, .(
        delta_bar = delta_bar[1],
        sd_star   = sd_star[1],
        n_bets    = n_bets[1]
      ), by = pid]
      setkey(pid_map, pid)
      pid_map <- pid_map[.(pid_levels)]
      stopifnot(!anyNA(pid_map$pid))
      
      delta_bar_vec <- as.numeric(pid_map$delta_bar)
      sd_star_vec   <- as.numeric(pid_map$sd_star)
      n_bets_vec    <- as.integer(pid_map$n_bets)
      
      # n_trials by seq within this (tr,tag) dataset slice
      d_counts <- prep[pid %in% pid_levels]
      n_trials_by_seq <- d_counts[, .(n_trials = .N), by = seq]
      setkey(n_trials_by_seq, seq)
      
      # extract posterior
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      
      stopifnot(!is.null(post$alpha), !is.null(post$u), !is.null(post$b))
      alpha <- as.numeric(post$alpha)   # K
      u     <- post$u                   # K x N
      b     <- post$b                   # K x S (sequence effects)
      
      K <- length(alpha)
      N <- length(pid_levels)
      S <- length(seq_levels)
      
      stopifnot(is.matrix(u), is.matrix(b))
      stopifnot(nrow(u) == K, ncol(u) == N)
      stopifnot(nrow(b) == K, ncol(b) == S)
      
      # tag -> output type label
      type_nm <- if (tag == "full") "exploratory" else "confirmatory"
      
      # ============================================================
      # SEQUENCES: mu_a (mapped back to absolute scale, then /e)
      # mu_s[k,s] = mean_i [ delta_bar_i + (alpha[k] + u[k,i] + b[k,s]) * sd_star_i ] / e
      # ============================================================
      eta_base <- sweep(u, 1, alpha, "+")  # K x N
      mu_s_draws <- matrix(NA_real_, nrow = K, ncol = S)
      
      for (s in seq_len(S)) {
        eta_mat <- eta_base + b[, s]                 # K x N (recycled over cols)
        mapped  <- sweep(eta_mat, 2, sd_star_vec, "*")
        mapped  <- sweep(mapped,  2, delta_bar_vec, "+")
        mu_s_draws[, s] <- rowMeans(mapped) / e
      }
      
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
      
      # simple directional label using main rho
      colU_main <- paste0("U_a_rho_", rho_nm_main)
      colO_main <- paste0("O_a_rho_", rho_nm_main)
      seq_tbl[, calib_label :=
                fifelse(get(colU_main) >= 0.95, "under_strong",
                        fifelse(get(colU_main) >= 0.80, "under_moderate",
                                fifelse(get(colO_main) >= 0.95, "over_strong",
                                        fifelse(get(colO_main) >= 0.80, "over_moderate", "neutral"))))]
      
      setorder(seq_tbl, sequence)
      
      f_seq_csv <- file.path(out_dir, paste0("rq2_", tr, "_", type_nm, "_sequences.csv"))
      if (!should_skip(
        paths = f_seq_csv,
        cfg   = cfg,
        type  = "output",
        label = paste0("RQ2 sequences ", type_nm, " (", ds, "/", tr, ")")
      )) {
        fwrite(seq_tbl, f_seq_csv)
        msg("Saved: ", f_seq_csv)
      }
      
      # ============================================================
      # PARTICIPANTS: mu_a (participant mean deviation /e)
      # mu_i[k,i] = (delta_bar_i + (alpha[k] + u[k,i]) * sd_star_i) / e
      # (sequence effects average out due to sum-to-zero b)
      # ============================================================
      mu_i_draws <- (sweep(eta_base, 2, sd_star_vec, "*") + rep(delta_bar_vec, each = K)) / e
      stopifnot(is.matrix(mu_i_draws), nrow(mu_i_draws) == K, ncol(mu_i_draws) == N)
      
      part_tbl <- data.table(
        pid         = pid_levels,
        n_bets      = n_bets_vec,
        mu_a_median = apply(mu_i_draws, 2, median),
        mu_a_mean   = apply(mu_i_draws, 2, mean),
        mu_a_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
        mu_a_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975)
      )
      
      for (rho in rho_vec) {
        nmU <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        nmO <- paste0("O_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        part_tbl[, (nmU) := apply(mu_i_draws, 2, function(x) mean(x < -rho))]
        part_tbl[, (nmO) := apply(mu_i_draws, 2, function(x) mean(x >  rho))]
      }
      
      colU_i_main <- paste0("U_i_rho_", rho_nm_main)
      colO_i_main <- paste0("O_i_rho_", rho_nm_main)
      part_tbl[, calib_label :=
                 fifelse(get(colU_i_main) >= 0.95, "under_solid",
                         fifelse(get(colU_i_main) >= 0.90, "under_likely",
                                 fifelse(get(colO_i_main) >= 0.95, "over_solid",
                                         fifelse(get(colO_i_main) >= 0.90, "over_likely", "neutral"))))]
      
      setorder(part_tbl, pid)
      
      f_pid_csv <- file.path(out_dir, paste0("rq2_", tr, "_", type_nm, "_participants.csv"))
      if (!should_skip(
        paths = f_pid_csv,
        cfg   = cfg,
        type  = "output",
        label = paste0("RQ2 participants ", type_nm, " (", ds, "/", tr, ")")
      )) {
        fwrite(part_tbl, f_pid_csv)
        msg("Saved: ", f_pid_csv)
      }
      
      outputs[[paste(tr, tag, sep = "_")]] <- list(
        fit = f_fit, pid = f_pid, seq = f_seq, prep = f_prep,
        sequences_csv = f_seq_csv,
        participants_csv = f_pid_csv
      )
    }
  }
  
  invisible(outputs)
}
