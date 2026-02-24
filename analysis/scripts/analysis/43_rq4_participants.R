# ============================================================
# scripts/analysis/43_rq4_participants.R
# Participant table from RQ4 fit (sequence-averaged Heads tendency)
#
# Prereg-compliant:
# - mu_hi^(k) = mean_s inv_logit(alpha^(k) + u_i^(k) + beta_s^(k))
# - hbar^(k)  = inv_logit(alpha^(k))
# - H_i(delta) = P(mu_hi > hbar + delta)
# - T_i(delta) = P(mu_hi < hbar - delta)
# - P(|mu_hi - 0.5| > delta)
# - Descriptive labels: solid/likely/leaning Headish/Tailish, else neutral
#
# Writes: data/clean/<ds>/output/rq4_participants_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

rq4_participants <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # delta settings (main + sensitivity)
  # ----------------------------
  # Prefer prereg-like vector if you store it, else fall back to (0.05, 0.03, 0.08).
  delta_vec <- NULL
  if (!is.null(cfg$design$rhos)) {
    if (!is.null(cfg$design$rhos$rq4_delta)) delta_vec <- cfg$design$rhos$rq4_delta
    if (!is.null(cfg$design$rhos$rq4_rho))   delta_vec <- cfg$design$rhos$rq4_rho  # tolerate alt name
  }
  if (is.null(delta_vec)) delta_vec <- c(0.05, 0.03, 0.08)
  
  delta_vec <- as.numeric(delta_vec)
  stopifnot(length(delta_vec) >= 1L, all(is.finite(delta_vec)), all(delta_vec > 0))
  delta_main <- delta_vec[1]
  
  # ----------------------------
  # paths
  # ----------------------------
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # label helper per prereg thresholds
  classify_one <- function(H, T) {
    # Headish
    if (H >= 0.95) return("solid_headish")
    if (H >= 0.90) return("likely_headish")
    if (H >= 0.75) return("leaning_headish")
    # Tailish
    if (T >= 0.95) return("solid_tailish")
    if (T >= 0.90) return("likely_tailish")
    if (T >= 0.75) return("leaning_tailish")
    # Neutral
    "neutral"
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, ".rds"))
    if (!file.exists(f_fit)) stop("RQ4 fit not found: ", f_fit)
    if (!file.exists(f_pid)) stop("RQ4 pid_levels not found: ", f_pid)
    
    fit <- readRDS(f_fit)
    pid_levels <- readRDS(f_pid)
    pid_levels <- as.character(pid_levels)
    
    post <- rstan::extract(fit)
    
    # Required for prereg participant summaries
    if (is.null(post$alpha)) stop("RQ4 fit missing parameter alpha.")
    if (is.null(post$u))     stop("RQ4 fit missing transformed parameter u (participant effects).")
    if (is.null(post$beta))  stop("RQ4 fit missing transformed parameter beta (sequence effects).")
    
    alpha <- post$alpha                 # iters
    u_mat <- post$u                     # iters x N
    b_mat <- post$beta                  # iters x S
    stopifnot(ncol(u_mat) == length(pid_levels))
    stopifnot(length(alpha) == nrow(u_mat))
    stopifnot(nrow(b_mat) == length(alpha))
    
    iters <- length(alpha)
    N <- length(pid_levels)
    S <- ncol(b_mat)
    
    # ----------------------------------------
    # mu_hi draws: mu_hi^(k) = mean_s inv_logit(alpha^(k) + u_i^(k) + beta_s^(k))
    # and baseline hbar^(k) = inv_logit(alpha^(k))
    # ----------------------------------------
    mu_hi <- matrix(NA_real_, nrow = iters, ncol = N)
    hbar  <- inv_logit(alpha)
    
    # Loop over draws (iters) is clean and fast enough (pilot-scale).
    for (k in seq_len(iters)) {
      # vector length S: alpha_k + beta_k,s
      ab <- alpha[k] + b_mat[k, ]
      # for each i: mean_s logistic(ab_s + u_ki)
      for (i in seq_len(N)) {
        mu_hi[k, i] <- mean(inv_logit(ab + u_mat[k, i]))
      }
    }
    
    # ----------------------------------------
    # Posterior summaries for mu_hi
    # ----------------------------------------
    tbl <- data.table(
      pid            = pid_levels,
      mu_h_median    = apply(mu_hi, 2, median),
      mu_h_mean      = apply(mu_hi, 2, mean),
      mu_h_q025      = apply(mu_hi, 2, quantile, probs = 0.025),
      mu_h_q975      = apply(mu_hi, 2, quantile, probs = 0.975)
    )
    
    # Also include baseline hbar summaries (same for all participants)
    tbl[, hbar_median := median(hbar)]
    tbl[, hbar_mean   := mean(hbar)]
    tbl[, hbar_q025   := as.numeric(quantile(hbar, probs = 0.025))]
    tbl[, hbar_q975   := as.numeric(quantile(hbar, probs = 0.975))]
    
    # ----------------------------------------
    # Prereg probabilities + labels (main delta)
    # H_i(delta) = P(mu_hi > hbar + delta)
    # T_i(delta) = P(mu_hi < hbar - delta)
    # P(|mu_hi - 0.5| > delta)
    # ----------------------------------------
    H_main <- numeric(N)
    T_main <- numeric(N)
    A_main <- numeric(N)
    
    for (i in seq_len(N)) {
      H_main[i] <- mean(mu_hi[, i] > (hbar + delta_main))
      T_main[i] <- mean(mu_hi[, i] < (hbar - delta_main))
      A_main[i] <- mean(abs(mu_hi[, i] - 0.5) > delta_main)
    }
    
    tbl[, `:=`(
      delta_main       = delta_main,
      H_main           = H_main,
      T_main           = T_main,
      P_absdev05_main  = A_main,
      side_label       = vapply(seq_len(N), function(i) classify_one(H_main[i], T_main[i]), character(1))
    )]
    
    # ----------------------------------------
    # Sensitivity deltas (if provided)
    # Adds: H_delta_*, T_delta_*, P_absdev05_delta_*
    # ----------------------------------------
    if (length(delta_vec) > 1L) {
      for (dlt in delta_vec) {
        nm <- gsub("\\.", "", sprintf("%.2f", dlt))
        
        H <- numeric(N)
        T <- numeric(N)
        A <- numeric(N)
        
        for (i in seq_len(N)) {
          H[i] <- mean(mu_hi[, i] > (hbar + dlt))
          T[i] <- mean(mu_hi[, i] < (hbar - dlt))
          A[i] <- mean(abs(mu_hi[, i] - 0.5) > dlt)
        }
        
        tbl[, (paste0("H_delta_", nm)) := H]
        tbl[, (paste0("T_delta_", nm)) := T]
        tbl[, (paste0("P_absdev05_delta_", nm)) := A]
      }
    }
    
    setorder(tbl, pid)
    
    f_out <- file.path(out_dir, paste0("rq4_participants_", tr, ".csv"))
    fwrite(tbl, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}