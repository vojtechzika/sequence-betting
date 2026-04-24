# ============================================================
# 12_rq1_diagnostics.R
#
# PURPOSE
#   Posterior predictive check for RQ1 Bernoulli likelihood adequacy.
#   Detects over/under-dispersion in sequence-wise betting frequencies.
#   Flags Beta-Binomial robustness if PPC is inadequate.
#
# METHOD
#   For each treatment x tag:
#   - Compute observed per-sequence betting rates p_obs[s]
#   - Draw K_ppc posterior predictive replications
#   - Compute dispersion statistic D = var_s(p_s) for observed and replicated
#   - PPC p-value: P(D_rep >= D_obs)
#   - Adequate if p in [p_lo, p_hi]
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/rq1_fit_sequences_<tr>_<tag>.rds
#   path_mod/rq1_pid_levels_<tr>_<tag>.rds
#   path_mod/rq1_seq_levels_<tr>_<tag>.rds
#   path_mod/drift_decisions.rds
#
# OUTPUT
#   path_out/rq1_diagnostics.csv
#
# NOTES
#   - Adequate flag controls whether rq1_stan(robustness=TRUE) fits BB model
#   - PPC interval [p_lo, p_hi] from cfg$model$ppc$rq1_interval
# ============================================================

rq1_diagnostics <- function(cfg) {
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  stopifnot(!is.null(cfg$model), !is.null(cfg$model$ppc))
  
  K_ppc <- as.integer(cfg$model$ppc$rq1_k)
  stopifnot(length(K_ppc) == 1L, is.finite(K_ppc), K_ppc >= 50L)
  
  p_lo <- as.numeric(cfg$model$ppc$rq1_interval[1])
  p_hi <- as.numeric(cfg$model$ppc$rq1_interval[2])
  stopifnot(is.finite(p_lo), is.finite(p_hi), p_lo > 0, p_hi < 1, p_lo < p_hi)
  
  f_out   <- file.path(path_out, "rq1_diagnostics.csv")
  f_drift <- file.path(path_mod, "drift_decisions.rds")
  
  if (should_skip(f_out, cfg, "output", "RQ1 diagnostics")) return(invisible(NULL))
  
  stopifnot(file.exists(f_drift))
  drift_decisions <- readRDS(f_drift)
  
  f_master <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  dt <- fread(f_master, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "seq", "block") %in% names(dt)))
  
  dt[, pid     := as.character(pid)]
  dt[, treat   := as.character(treat)]
  dt[, seq     := as.character(seq)]
  dt[, stake   := as.numeric(stake)]
  dt[is.na(stake), stake := 0]
  dt[, y       := as.integer(stake > 0)]
  dt[, block   := as.integer(block)]
  dt[, block_c := as.numeric(block) - 2.5]
  
  tags     <- c("full", "confirmatory")
  all_rows <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      f_fit <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_", tag, ".rds"))
      f_pid <- file.path(path_mod, paste0("rq1_pid_levels_",    tr, "_", tag, ".rds"))
      f_seq <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_", tag, ".rds"))
      
      if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) {
        warning("RQ1 diagnostics: missing Stan artifacts for tr='", tr,
                "', tag='", tag, "'. Skipping.")
        next
      }
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      stopifnot(length(seq_levels) == 64L)
      
      d <- dt[treat == tr & pid %in% pid_levels]
      if (nrow(d) == 0) next
      
      d[, pid_i := match(pid, pid_levels)]
      d[, sid_s := match(seq, seq_levels)]
      stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
      
      N  <- length(pid_levels)
      S  <- length(seq_levels)
      Tn <- nrow(d)
      
      # ---- Drift config ----
      key           <- paste(tr, tag, "bet", sep = "_")
      drift_res     <- drift_decisions[[key]]
      include_drift <- !is.null(drift_res) && isTRUE(drift_res$drift)
      drift_type    <- if (!is.null(drift_res)) drift_res$drift_type else "none"
      
      # ---- Observed dispersion ----
      obs_tbl <- d[, .(p_obs = mean(y)), by = seq]
      setkey(obs_tbl, seq)
      p_obs <- obs_tbl[.(seq_levels), p_obs]
      p_obs[is.na(p_obs)] <- 0
      D_obs <- stats::var(p_obs)
      
      # ---- Load posterior ----
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      stopifnot(!is.null(post$alpha), !is.null(post$u), !is.null(post$b),
                !is.null(post$gamma_drift))
      
      alpha       <- as.numeric(post$alpha)
      u           <- post$u
      b           <- post$b
      gamma_drift <- as.numeric(post$gamma_drift)
      K_all       <- length(alpha)
      
      stopifnot(is.matrix(u), nrow(u) == K_all, ncol(u) == N)
      stopifnot(is.matrix(b), nrow(b) == K_all, ncol(b) == S)
      
      K_use <- min(K_ppc, K_all)
      set.seed(seed + 1001L)
      k_idx <- sort(sample.int(K_all, K_use, replace = FALSE))
      
      pid_i    <- d$pid_i
      sid_s    <- d$sid_s
      block_c  <- d$block_c
      idx_by_s <- split(seq_len(Tn), sid_s)
      
      # ---- PPC replications ----
      set.seed(seed + 2000L)
      D_rep <- numeric(K_use)
      
      for (jj in seq_len(K_use)) {
        k   <- k_idx[jj]
        eta <- alpha[k] + u[k, pid_i] + b[k, sid_s]
        if (include_drift) eta <- eta + gamma_drift[k] * block_c
        y_rep <- rbinom(Tn, size = 1L, prob = plogis(eta))
        
        p_rep_s <- numeric(S)
        for (s in seq_len(S)) {
          idx <- idx_by_s[[as.character(s)]]
          p_rep_s[s] <- if (length(idx) == 0L) 0 else mean(y_rep[idx])
        }
        D_rep[jj] <- stats::var(p_rep_s)
      }
      
      ppc_p    <- mean(D_rep >= D_obs)
      adequate <- (ppc_p >= p_lo && ppc_p <= p_hi)
      
      if (adequate) {
        msg("RQ1 PPC (", tr, "/", tag, "): adequate | p=",
            sprintf("%.3f", ppc_p), " | D_obs=", signif(D_obs, 4),
            " | drift=", drift_type)
      } else {
        warning("RQ1 PPC (", tr, "/", tag,
                "): INADEQUATE -- Beta-Binomial robustness will be fitted | p=",
                sprintf("%.3f", ppc_p), " | D_obs=", signif(D_obs, 4),
                " | drift=", drift_type)
      }
      
      all_rows[[paste(tr, tag, sep = "_")]] <- data.table(
        treatment    = tr,
        tag          = tag,
        drift        = drift_type,
        N            = N,
        S            = S,
        T            = Tn,
        K_all        = K_all,
        K_used       = K_use,
        D_obs        = D_obs,
        D_rep_median = stats::median(D_rep),
        D_rep_q025   = stats::quantile(D_rep, 0.025),
        D_rep_q975   = stats::quantile(D_rep, 0.975),
        ppc_p        = ppc_p,
        adequate     = adequate,
        p_lo         = p_lo,
        p_hi         = p_hi
      )
    }
  }
  
  if (length(all_rows) > 0) {
    fwrite(rbindlist(all_rows), f_out)
    msg("Saved: ", f_out)
  }
  
  invisible(all_rows)
}