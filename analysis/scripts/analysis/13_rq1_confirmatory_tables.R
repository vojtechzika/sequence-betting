# ============================================================
# scripts/analysis/13_rq1_confirmatory_tables.R
#
# RQ1 confirmatory tables (subset-averaged; confirmatory treatment only)
#
# What it does:
#   - Uses the FULL RQ1 Stan fit for the confirmatory treatment (defined as
#     the first element of cfg$run$treatment).
#   - Recomputes sequence-level mu_b[s] within each posterior draw, averaging
#     ONLY over the benchmark-consistent subset of participants (keep set).
#   - Produces:
#       (i) sequence table (subset-averaged mu_b by sequence)
#       (ii) participant table (subset participants only; mu_b_i by participant)
#
# Inputs:
#   models/rq1_fit_sequences_<tr>.rds
#   models/rq1_pid_levels_<tr>.rds
#   models/rq1_seq_levels_<tr>.rds
#   models/a_star_pid_flags_<tr>.rds   (from scripts/indices/05_pids_with_positive_a_star.R)
#   clean/<ds>/master_sequences.csv    (for n_trials)
#
# Config requirements:
#   cfg$run$treatment : character vector, first element is confirmatory treatment
#   cfg$design$rhos$rq1_rho : numeric vector, first is main rho
#   cfg$design$a_flags$tau : numeric vector, first is main tau (used here)
#
# Outputs:
#   output/rq1_sequences_confirmatory_<tr>.csv
#   output/rq1_participants_confirmatory_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

rq1_confirmatory_tables <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  ds <- as.character(cfg$run$dataset)
  
  # confirmatory treatment = first listed
  stopifnot(!is.null(cfg$run$treatment))
  tr_all <- as.character(cfg$run$treatment)
  stopifnot(length(tr_all) >= 1L, all(nzchar(tr_all)))
  tr <- tr_all[1]
  
  # rho thresholds
  stopifnot(!is.null(cfg$design$rhos), !is.null(cfg$design$rhos$rq1_rho))
  rho_vec <- as.numeric(cfg$design$rhos$rq1_rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  rho_main <- rho_vec[1]
  
  # tau thresholds (first is main)
  stopifnot(!is.null(cfg$design$a_flags), !is.null(cfg$design$a_flags$tau))
  tau_vec <- as.numeric(cfg$design$a_flags$tau)
  stopifnot(length(tau_vec) >= 1L, all(is.finite(tau_vec)), all(tau_vec > 0), all(tau_vec < 1))
  tau_main <- tau_vec[1]
  tau_nm <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- inputs: RQ1 artifacts ---
  f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
  f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, ".rds"))
  f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, ".rds"))
  stopifnot(file.exists(f_fit), file.exists(f_pid), file.exists(f_seq))
  
  # --- inputs: pid flags (source of truth for subset) ---
  f_flags <- file.path(mod_dir, paste0("a_star_pid_flags_", tr, ".rds"))
  stopifnot(file.exists(f_flags))
  
  # --- master for n_trials ---
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
  master[, y := as.integer(stake > 0)]
  
  d <- master[treat == tr]
  if (nrow(d) == 0) stop("RQ1 confirmatory: no rows for treat='", tr, "' in master_sequences.csv.")
  
  # --- outputs + overwrite control ---
  f_seq_out <- file.path(out_dir, paste0("rq1_sequences_confirmatory_", tr, ".csv"))
  f_pid_out <- file.path(out_dir, paste0("rq1_participants_confirmatory_", tr, ".csv"))
  
  if (should_skip(
    paths = c(f_seq_out, f_pid_out),
    cfg   = cfg,
    type  = "output",
    label = paste0("RQ1 confirmatory tables (", ds, "/", tr, ")")
  )) {
    return(invisible(NULL))
  }
  
  # --- levels ---
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
         "treatment=", tr, " | tau_main=", tau_main, "\n",
         "flags file=", f_flags)
  }
  
  keep_idx <- match(keep_pid, pid_levels)
  stopifnot(all(!is.na(keep_idx)))
  
  # --- extract posterior draws needed to recompute subset-based means ---
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
  
  # ============================================================
  # SEQUENCES (subset-averaged mu_b[s])
  # mu_b_sub[k,s] = mean_{i in subset} inv_logit(alpha[k] + u[k,i] + b[k,s])
  # ============================================================
  mu_b_sub <- matrix(NA_real_, nrow = K, ncol = S)
  
  for (k in seq_len(K)) {
    eta_i <- alpha[k] + u[k, keep_idx]          # length |subset|
    eta_mat <- outer(eta_i, b[k, ], "+")        # |subset| x S
    mu_b_sub[k, ] <- colMeans(plogis(eta_mat))  # length S
  }
  
  # n_trials per sequence (observed)
  n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
  setkey(n_trials_by_seq, seq)
  
  seq_tbl <- data.table(
    sequence    = seq_levels,
    mu_b_median = apply(mu_b_sub, 2, median),
    mu_b_mean   = apply(mu_b_sub, 2, mean),
    mu_b_q025   = apply(mu_b_sub, 2, quantile, probs = 0.025),
    mu_b_q975   = apply(mu_b_sub, 2, quantile, probs = 0.975)
  )
  
  seq_tbl[, n_trials := n_trials_by_seq[.(sequence), n_trials]]
  seq_tbl[, n_trials := as.integer(n_trials)]
  
  for (rho in rho_vec) {
    nm <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
    thr <- 1 - rho
    seq_tbl[, (nm) := apply(mu_b_sub, 2, function(x) mean(x < thr))]
  }
  
  col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
  seq_tbl[, underbet_label :=
            fifelse(get(col_main) >= 0.95, "strong",
                    fifelse(get(col_main) >= 0.80, "moderate",
                            fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
  
  setorder(seq_tbl, sequence)
  
  # ============================================================
  # PARTICIPANTS (subset participants only)
  # mu_i_sub[k,j] = mean_s inv_logit(alpha[k] + u[k,i_j] + b[k,s])
  # ============================================================
  mu_i_sub <- matrix(NA_real_, nrow = K, ncol = length(keep_idx))
  
  for (j in seq_along(keep_idx)) {
    i <- keep_idx[j]
    # K x S matrix
    eta_mat <- b
    eta_mat <- sweep(eta_mat, 1, alpha + u[, i], "+")
    mu_i_sub[, j] <- rowMeans(plogis(eta_mat))
  }
  
  n_trials_by_pid <- d[, .(n_trials = .N), by = pid]
  setkey(n_trials_by_pid, pid)
  
  pid_tbl <- data.table(
    pid         = keep_pid,
    mu_i_median = apply(mu_i_sub, 2, median),
    mu_i_mean   = apply(mu_i_sub, 2, mean),
    mu_i_q025   = apply(mu_i_sub, 2, quantile, probs = 0.025),
    mu_i_q975   = apply(mu_i_sub, 2, quantile, probs = 0.975)
  )
  
  for (rho in rho_vec) {
    nm <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
    thr <- 1 - rho
    pid_tbl[, (nm) := apply(mu_i_sub, 2, function(x) mean(x < thr))]
  }
  
  col_main_i <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
  pid_tbl[, underbet_label :=
            fifelse(get(col_main_i) >= 0.95, "solid",
                    fifelse(get(col_main_i) >= 0.90, "likely",
                            fifelse(get(col_main_i) >= 0.75, "leaning", "neutral")))]
  
  pid_tbl[, n_trials := n_trials_by_pid[.(pid), n_trials]]
  pid_tbl[, n_trials := as.integer(n_trials)]
  setorder(pid_tbl, pid)
  
  # --- write ---
  fwrite(seq_tbl, f_seq_out)
  msg("Saved: ", f_seq_out)
  
  fwrite(pid_tbl, f_pid_out)
  msg("Saved: ", f_pid_out)
  
  invisible(list(
    dataset = ds,
    confirmatory_treatment = tr,
    tau_main = tau_main,
    rho_main = rho_main,
    subset_n = length(keep_pid),
    sequences = seq_tbl,
    participants = pid_tbl
  ))
}

# Example:
# source(here::here("scripts","analysis","13_rq1_confirmatory_tables.R"))
# rq1_confirmatory_tables(cfg)