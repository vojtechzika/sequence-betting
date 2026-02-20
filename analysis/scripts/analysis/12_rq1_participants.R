# scripts/analysis/22_rq1_participants.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq1_participants <- function(dataset = "pilot",
                             treat_fn = "m25",
                             rho = 0.10,
                             seed = 12345) {
  
  # ----------------------------
  # Files
  # ----------------------------
  infile <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  # Fit + levels produced by 21_
  fit_candidates <- c(
    file.path(path_mod_ds(dataset), paste0("rq1_fit_sequences_", treat_fn, ".rds")),
    file.path(path_mod_ds(dataset), "rq1_sequences_fit.rds")
  )
  f_fit <- fit_candidates[file.exists(fit_candidates)][1]
  if (is.na(f_fit) || !nzchar(f_fit)) {
    stop("Could not find RQ1 fit file. Tried:\n", paste(fit_candidates, collapse = "\n"))
  }
  
  lvl_candidates_pid <- c(
    file.path(path_mod_ds(dataset), paste0("rq1_pid_levels_", dataset, "_", treat_fn, ".rds")),
    file.path(path_mod_ds(dataset), paste0("rq1_pid_levels_", treat_fn, ".rds")),
    file.path(path_mod_ds(dataset), "rq1_pid_levels.rds")
  )
  lvl_candidates_seq <- c(
    file.path(path_mod_ds(dataset), paste0("rq1_seq_levels_", dataset, "_", treat_fn, ".rds")),
    file.path(path_mod_ds(dataset), paste0("rq1_seq_levels_", treat_fn, ".rds")),
    file.path(path_mod_ds(dataset), "rq1_seq_levels.rds")
  )
  
  f_pid_levels <- lvl_candidates_pid[file.exists(lvl_candidates_pid)][1]
  f_seq_levels <- lvl_candidates_seq[file.exists(lvl_candidates_seq)][1]
  
  if (is.na(f_pid_levels) || !nzchar(f_pid_levels)) {
    stop("Could not find pid_levels file. Tried:\n", paste(lvl_candidates_pid, collapse = "\n"))
  }
  if (is.na(f_seq_levels) || !nzchar(f_seq_levels)) {
    stop("Could not find seq_levels file. Tried:\n", paste(lvl_candidates_seq, collapse = "\n"))
  }
  
  # ----------------------------
  # Load + schema checks
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "stake", "seq")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      "master_sequences.csv missing columns: ", paste(missing, collapse = ", "),
      "\nExpected at least: pid, treat, stake, seq"
    )
  }
  
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  
  # ----------------------------
  # Filter to FN (confirmatory)
  # ----------------------------
  d <- dt[treat == treat_fn]
  if (nrow(d) == 0) stop("No rows after filtering treat == '", treat_fn, "'")
  
  d[, y := as.integer(stake > 0)]
  
  # ----------------------------
  # Load levels (must match fit)
  # ----------------------------
  pid_levels <- readRDS(f_pid_levels)
  seq_levels <- readRDS(f_seq_levels)
  
  # Map ids using saved levels
  d[, pid_i := match(pid, pid_levels)]
  d[, sid_s := match(seq, seq_levels)]
  
  if (anyNA(d$pid_i)) {
    bad <- unique(d[is.na(pid_i), pid])
    stop("Found pids not present in saved pid_levels (21_). Example: ", paste(head(bad, 10), collapse = ", "))
  }
  if (anyNA(d$sid_s)) {
    bad <- unique(d[is.na(sid_s), seq])
    stop("Found seq values not present in saved seq_levels (21_). Example: ", paste(head(bad, 10), collapse = ", "))
  }
  
  # ----------------------------
  # Load fit + extract draws
  # ----------------------------
  fit <- readRDS(f_fit)
  post <- rstan::extract(fit)
  
  alpha_draws <- post$alpha  # iters
  u_draws     <- post$u      # iters x N
  b_draws     <- post$b      # iters x S
  
  iters <- length(alpha_draws)
  N <- length(pid_levels)
  S <- length(seq_levels)
  
  stopifnot(is.matrix(u_draws), is.matrix(b_draws))
  stopifnot(nrow(u_draws) == iters, nrow(b_draws) == iters)
  stopifnot(ncol(u_draws) == N, ncol(b_draws) == S)
  
  # ----------------------------
  # Participant-level mu_i^b draws: mean over sequences
  # mu_i^b(m) = mean_s inv_logit(alpha(m) + u_i(m) + b_s(m))
  # ----------------------------
  mu_i_draws <- matrix(NA_real_, nrow = iters, ncol = N)
  
  for (i in 1:N) {
    eta <- alpha_draws + u_draws[, i]            # iters
    eta_mat <- sweep(b_draws, 1, eta, FUN = "+") # iters x S
    mu_i_draws[, i] <- rowMeans(plogis(eta_mat)) # iters
  }
  
  # ----------------------------
  # Summaries per participant
  # ----------------------------
  thr <- 1 - rho
  
  out <- data.table(
    pid         = pid_levels,
    mu_i_median = apply(mu_i_draws, 2, median),
    mu_i_mean   = apply(mu_i_draws, 2, mean),
    mu_i_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
    mu_i_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975),
    U_i_rho     = apply(mu_i_draws, 2, function(x) mean(x < thr))
  )
  
  out[, underbet_label :=
        fifelse(U_i_rho >= 0.95, "solid",
                fifelse(U_i_rho >= 0.90, "likely",
                        fifelse(U_i_rho >= 0.75, "leaning", "neutral")))]
  
  # Add n_trials per participant (from data)
  n_trials_by_pid <- d[, .(n_trials = .N), by = pid]
  setkey(n_trials_by_pid, pid)
  out[, n_trials := n_trials_by_pid[.(pid), n_trials]]
  out[, n_trials := as.integer(n_trials)]
  
  setorder(out, pid)
  
  # ----------------------------
  # Save
  # ----------------------------
  out_csv <- file.path(
    path_out_ds(dataset),
    paste0("rq1_participants_", dataset, "_", treat_fn, ".csv")
  )
  fwrite(out, out_csv)
  msg("Saved RQ1 participant table:", out_csv)
  
  invisible(list(
    participants = out,
    mu_i_draws = mu_i_draws,
    fit_file = f_fit,
    pid_levels_file = f_pid_levels,
    seq_levels_file = f_seq_levels
  ))
}

# Example:
# rq1_participants("pilot", treat_fn = "m25", rho = 0.10)
# rq1_participants("main",  treat_fn = "m25", rho = 0.10)