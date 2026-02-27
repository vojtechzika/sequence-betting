# ============================================================
# scripts/analysis/12_rq1_fit_check.R
#
# RQ1: Bernoulli adequacy check (PPC-style; sequence-wise dispersion)
#
# Goal (per prereg):
#   Detect systematic lack of fit of Bernoulli likelihood for b_is,
#   specifically over/under-dispersion in sequence-wise betting frequencies.
#
# Method:
#   - For each treatment, read:
#       models/rq1_fit_sequences_<tr>.rds
#       models/rq1_pid_levels_<tr>.rds
#       models/rq1_seq_levels_<tr>.rds
#     plus master_sequences.csv.
#   - Compute observed per-sequence betting rates p_obs[s].
#   - Posterior predictive:
#       y_rep[t] ~ Bernoulli_logit(alpha + u[pid[t]] + b[sid[t]])
#     Aggregate to p_rep[k,s] and compute dispersion statistic:
#       D = var_s(p_s)
#   - Posterior predictive p-value:
#       p = P(D_rep >= D_obs)
#   - Print:
#       "RQ1 Bernoulli PPC: adequate" if p in [0.05, 0.95]
#       else "inadequate (flag beta-binomial robustness)"
#
# Output:
#   data/clean/<ds>/output/rq1_fit_check_<tr>.csv
#     with D_obs, ppc_p, interval, K_used, T, S, N
# ============================================================

library(data.table)
library(rstan)

rq1_fit_check <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$plan))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # --- run controls ---
  # Number of posterior draws used for PPC (subsample for speed)
  K_ppc <- cfg$run$rq1_ppc_k
  if (is.null(K_ppc)) K_ppc <- 500L
  
  K_ppc <- as.integer(K_ppc)
  stopifnot(
    length(K_ppc) == 1L,
    is.finite(K_ppc),
    K_ppc >= 50L
  )
  
  # p-value interval for "adequate"
  p_lo <- cfg$run$rq1_ppc_p_lo
  p_hi <- cfg$run$rq1_ppc_p_hi
  if (is.null(p_lo)) p_lo <- 0.05
  if (is.null(p_hi)) p_hi <- 0.95
  p_lo <- as.numeric(p_lo); p_hi <- as.numeric(p_hi)
  stopifnot(is.finite(p_lo), is.finite(p_hi), p_lo > 0, p_hi < 1, p_lo < p_hi)
  
  # RNG seed for PPC subsampling / simulation
  seed <- cfg$run$seed
  if (is.null(seed)) seed <- 12345L
  seed <- as.integer(seed)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- load master once ---
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  dt <- fread(f_master, encoding = "UTF-8")
  
  required <- c("pid", "treat", "stake", "seq")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) stop("master_sequences.csv missing columns: ", paste(missing, collapse = ", "))
  
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[is.na(stake), stake := 0]
  dt[, y := as.integer(stake > 0)]  # RQ1 definition
  
  # helper: safe mean by group (ensures integer counts)
  seq_rate <- function(d) {
    # returns data.table(sequence, n_trials, p_obs)
    d[, .(n_trials = .N, p_obs = mean(y)), by = seq]
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_", tr, ".rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) {
      warning("RQ1 fit check: missing Stan artifacts for tr='", tr, "'. Skipping.")
      next
    }
    
    d <- dt[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ1 fit check: no rows in master for tr='", tr, "'. Skipping.")
      next
    }
    
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    
    # (Optional but helpful) enforce full design: 64 sequences expected
    # If your experiment always includes 64 sequences per treatment, keep this strict.
    # Otherwise, comment it out.
    stopifnot(length(seq_levels) == 64L)
    
    # align indices EXACTLY as used in Stan fitting
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    if (anyNA(d$pid_i)) stop("RQ1 fit check: pid mismatch vs pid_levels in tr='", tr, "'.")
    if (anyNA(d$sid_s)) stop("RQ1 fit check: seq mismatch vs seq_levels in tr='", tr, "'.")
    
    N <- length(pid_levels)
    S <- length(seq_levels)
    Tn <- nrow(d)
    
    # observed per-sequence betting rates (aligned to seq_levels)
    obs_tbl <- seq_rate(d)
    setkey(obs_tbl, seq)
    p_obs <- obs_tbl[.(seq_levels), p_obs]
    n_trials <- obs_tbl[.(seq_levels), n_trials]
    p_obs[is.na(p_obs)] <- 0  # should not happen if all sequences exist
    n_trials[is.na(n_trials)] <- 0L
    
    # dispersion statistic: variance across sequences of observed rates
    D_obs <- stats::var(p_obs)
    
    # extract posterior draws
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    if (is.null(post$alpha) || is.null(post$u) || is.null(post$b)) {
      stop("RQ1 fit check: fit missing alpha/u/b in tr='", tr, "'.")
    }
    
    alpha <- as.numeric(post$alpha) # K_all
    u <- post$u                     # K_all x N
    b <- post$b                     # K_all x S
    K_all <- length(alpha)
    
    stopifnot(is.matrix(u), is.matrix(b))
    stopifnot(nrow(u) == K_all, ncol(u) == N)
    stopifnot(nrow(b) == K_all, ncol(b) == S)
    
    # choose K draws for PPC
    if (is.na(K_ppc)) {
      K_use <- K_all
      k_idx <- seq_len(K_all)
    } else {
      K_use <- min(K_ppc, K_all)
      set.seed(seed + 1001L)
      k_idx <- sort(sample.int(K_all, K_use, replace = FALSE))
    }
    
    # precompute trial-level indices
    pid_i <- d$pid_i
    sid_s <- d$sid_s
    
    # Posterior predictive dispersion draws
    # We simulate y_rep per draw, then compute per-seq mean and variance across seq.
    # Implementation note: we compute per-seq sums via accumarray-like split.
    D_rep <- numeric(K_use)
    
    # For speed: split trial indices by sequence once
    idx_by_s <- split(seq_len(Tn), sid_s)
    
    # RNG for simulation
    set.seed(seed + 2000L)
    
    for (jj in seq_len(K_use)) {
      k <- k_idx[jj]
      
      # linear predictor per trial
      eta_t <- alpha[k] + u[k, pid_i] + b[k, sid_s]
      p_t <- plogis(eta_t)
      
      # simulate y_rep
      y_rep <- rbinom(Tn, size = 1L, prob = p_t)
      
      # per-sequence means
      p_rep_s <- numeric(S)
      for (s in seq_len(S)) {
        idx <- idx_by_s[[as.character(s)]]
        # idx should exist for every s if all sequences are present
        if (is.null(idx) || length(idx) == 0L) {
          p_rep_s[s] <- 0
        } else {
          p_rep_s[s] <- mean(y_rep[idx])
        }
      }
      
      D_rep[jj] <- stats::var(p_rep_s)
    }
    
    # PPC p-value: P(D_rep >= D_obs)
    ppc_p <- mean(D_rep >= D_obs)
    
    adequate <- (ppc_p >= p_lo && ppc_p <= p_hi)
    
    if (adequate) {
      msg("RQ1 Bernoulli PPC (tr=", tr, "): adequate | p=", sprintf("%.3f", ppc_p),
          " | D_obs=", signif(D_obs, 4))
    } else {
      warning("RQ1 Bernoulli PPC (tr=", tr, "): INADEQUATE (flag beta-binomial robustness) | p=",
              sprintf("%.3f", ppc_p), " | D_obs=", signif(D_obs, 4))
    }
    
    tbl <- data.table(
      dataset = ds,
      treatment = tr,
      N = N,
      S = S,
      T = Tn,
      K_all = K_all,
      K_used = K_use,
      D_obs = D_obs,
      D_rep_median = stats::median(D_rep),
      D_rep_q025 = stats::quantile(D_rep, 0.025),
      D_rep_q975 = stats::quantile(D_rep, 0.975),
      ppc_p = ppc_p,
      adequate = adequate,
      p_interval_lo = p_lo,
      p_interval_hi = p_hi
    )
    
    f_out <- file.path(out_dir, paste0("rq1_fit_check_", tr, ".csv"))
    
    if (should_skip(
      paths = f_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ1 fit check (", ds, "/", tr, ")")
    )) {
      next
    }
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}

# Example:
# source(here::here("scripts","analysis","12_rq1_fit_check.R"))
# rq1_fit_check(cfg)