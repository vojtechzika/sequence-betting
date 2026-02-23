# scripts/analysis/12_rq1_participants.R

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq1_participants <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  run    <- cfg$run
  design <- cfg$design
  model  <- cfg$model
  
  ds   <- as.character(run$dataset)
  seed <- as.integer(run$seed)
  
  # ----------------------------
  # Settings from cfg
  # ----------------------------
  stopifnot(!is.null(design$rhos), !is.null(design$rhos$rq1_rho))
  rho_vec <- as.numeric(design$rhos$rq1_rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  
  rho_main <- rho_vec[1]
  thr <- 1 - rho_main
  
  # treatments to run (already resolved against merged.csv)
  treat_list <- as.character(cfg$plan$by)
  stopifnot(length(treat_list) >= 1L, all(nzchar(treat_list)))
  
  # ----------------------------
  # Input data
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  required <- c("pid", "treat", "stake", "seq")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      "master_sequences.csv missing columns: ", paste(missing, collapse = ", "),
      "\nExpected at least: ", paste(required, collapse = ", ")
    )
  }
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, y     := as.integer(!is.na(stake) & stake > 0)]
  
  # ----------------------------
  # Helper: compute participant table from a saved fit
  # ----------------------------
  run_one <- function(tr) {
    
    # Fit + levels produced by 11_rq1_sequences.R (treatment-specific)
    f_fit        <- file.path(path_mod_ds(ds), paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(path_mod_ds(ds), paste0("rq1_pid_levels_",  tr, ".rds"))
    f_seq_levels <- file.path(path_mod_ds(ds), paste0("rq1_seq_levels_",  tr, ".rds"))
    
    if (!file.exists(f_fit)) {
      stop("RQ1 fit file not found: ", f_fit, "\nRun rq1_sequences(cfg) first.")
    }
    if (!file.exists(f_pid_levels)) {
      stop("pid_levels file not found: ", f_pid_levels, "\nRun rq1_sequences(cfg) first.")
    }
    if (!file.exists(f_seq_levels)) {
      stop("seq_levels file not found: ", f_seq_levels, "\nRun rq1_sequences(cfg) first.")
    }
    
    d <- dt[treat == tr]
    if (nrow(d) == 0) stop("No rows for treat='", tr, "' in master_sequences.csv (dataset='", ds, "').")
    
    pid_levels <- readRDS(f_pid_levels)
    seq_levels <- readRDS(f_seq_levels)
    
    # map using saved levels (must match fit)
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    if (anyNA(d$pid_i)) {
      bad <- unique(d[is.na(pid_i), pid])
      stop("Found pids not present in saved pid_levels for treat='", tr, "'. Example: ", paste(head(bad, 10), collapse = ", "))
    }
    if (anyNA(d$sid_s)) {
      bad <- unique(d[is.na(sid_s), seq])
      stop("Found seq values not present in saved seq_levels for treat='", tr, "'. Example: ", paste(head(bad, 10), collapse = ", "))
    }
    
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    # draws
    alpha_draws <- post$alpha  # iters
    u_draws     <- post$u      # iters x N
    b_draws     <- post$b      # iters x S
    
    iters <- length(alpha_draws)
    N <- length(pid_levels)
    S <- length(seq_levels)
    
    stopifnot(is.matrix(u_draws), is.matrix(b_draws))
    stopifnot(nrow(u_draws) == iters, nrow(b_draws) == iters)
    stopifnot(ncol(u_draws) == N, ncol(b_draws) == S)
    
    # participant-level mu_i^b draws: mean over sequences
    mu_i_draws <- matrix(NA_real_, nrow = iters, ncol = N)
    for (i in 1:N) {
      eta <- alpha_draws + u_draws[, i]             # iters
      eta_mat <- sweep(b_draws, 1, eta, FUN = "+")  # iters x S
      mu_i_draws[, i] <- rowMeans(plogis(eta_mat))
    }
    
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
    
    # n_trials per participant (from filtered data)
    n_trials_by_pid <- d[, .(n_trials = .N), by = pid]
    setkey(n_trials_by_pid, pid)
    out[, n_trials := n_trials_by_pid[.(pid), n_trials]]
    out[, n_trials := as.integer(n_trials)]
    
    setorder(out, pid)
    
    out_csv <- file.path(path_out_ds(ds), paste0("rq1_participants_", tr, ".csv"))
    fwrite(out, out_csv)
    
    msg(
      "Saved: ", out_csv, "\n",
      "  Inputs:\n",
      "    ", f_fit, "\n",
      "    ", f_pid_levels, "\n",
      "    ", f_seq_levels
    )
    
    invisible(list(
      participants = out,
      mu_i_draws = mu_i_draws,
      fit_file = f_fit,
      pid_levels_file = f_pid_levels,
      seq_levels_file = f_seq_levels,
      out_csv = out_csv
    ))
  }
  
  # ----------------------------
  # Run for each selected treatment
  # ----------------------------
  outputs <- list()
  out_files <- character(0)
  
  for (tr in treat_list) {
    res <- run_one(tr)
    outputs[[tr]] <- res
    out_files <- c(out_files, res$out_csv)
  }
  
  msg("\nRQ1 participant outputs produced:\n  ", paste(out_files, collapse = "\n  "))
  
  invisible(outputs)
}

# Example:
# rq1_participants(cfg)