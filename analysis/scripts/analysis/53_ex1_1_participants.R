# ============================================================
# scripts/analysis/53_ex1_1_participants.R
#
# EX1.I: Participant-level directional tendencies (GHI)
#
# chi_i = E_s[z_is]
# HH(.) = P(chi_i >  delta)
# G(.)  = P(chi_i < -delta)
#
# Classification:
#   solid / likely / leaning hot-handish
#   solid / likely / leaning gamblerish
#   neutral otherwise
#
# Uses anchor-weight construction from prereg:
#   d_a(s) = P(|theta_s - theta_a| < eps)
#   w_a(s) = normalize(d_a(s) + eta)
# and z_is = sign(h_is) * (wH(s) - wT(s))  (w0 cancels)
#
# Outputs per treatment:
#  - data/clean/<ds>/output/ex1_1_participants_<tr>.csv
#  - data/clean/<ds>/models/ex1_1_participants_<tr>.rds  (for EX1.2)
# ============================================================

source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

ex1_1_participants <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Labels from design
  # ----------------------------
  stopifnot(!is.null(design$seq), !is.null(design$seq$side_labels), !is.null(design$seq$anchor_labels))
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, length(lab_tails) == 1L, nzchar(lab_heads), nzchar(lab_tails), lab_heads != lab_tails)
  
  pure_heads <- as.character(design$seq$anchor_labels$pure_heads)
  pure_tails <- as.character(design$seq$anchor_labels$pure_tails)
  stopifnot(length(pure_heads) == 1L, length(pure_tails) == 1L, nzchar(pure_heads), nzchar(pure_tails), pure_heads != pure_tails)
  
  # ----------------------------
  # Main thresholds from design (with safe defaults)
  # ----------------------------
  delta <- 0.05
  if (!is.null(design$rhos) && !is.null(design$rhos$ex1_delta_main)) {
    delta <- as.numeric(design$rhos$ex1_delta_main)
  }
  stopifnot(length(delta) == 1L, is.finite(delta), delta > 0)
  
  eps <- 0.05
  if (!is.null(design$rhos) && !is.null(design$rhos$ex1_eps_main)) {
    eps <- as.numeric(design$rhos$ex1_eps_main)
  }
  stopifnot(length(eps) == 1L, is.finite(eps), eps > 0)
  
  eta <- 1e-6
  if (!is.null(design$rhos) && !is.null(design$rhos$ex1_eta)) {
    eta <- as.numeric(design$rhos$ex1_eta)
  }
  stopifnot(length(eta) == 1L, is.finite(eta), eta > 0)
  
  base_mc <- 2000L
  if (!is.null(cfg$run$ex1_base_mc)) base_mc <- as.integer(cfg$run$ex1_base_mc)
  stopifnot(length(base_mc) == 1L, base_mc >= 200L)
  
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  req <- c("pid", "treat", "seq", "stake", "side")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("master_sequences.csv missing: ", paste(miss, collapse = ", "))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    if (!file.exists(f_fit)) stop("EX1.1 participants: missing RQ4 fit: ", f_fit)
    if (!file.exists(f_seq)) stop("EX1.1 participants: missing RQ4 seq levels: ", f_seq)
    
    fit <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    
    post <- rstan::extract(fit)
    if (is.null(post$mu_h))    stop("EX1.1 participants: RQ4 fit missing mu_h.")
    if (is.null(post$alpha))   stop("EX1.1 participants: RQ4 fit missing alpha.")
    if (is.null(post$sigma_u)) stop("EX1.1 participants: RQ4 fit missing sigma_u.")
    
    mu_h    <- post$mu_h        # K x S
    alpha   <- as.numeric(post$alpha)
    sigma_u <- as.numeric(post$sigma_u)
    
    K <- nrow(mu_h)
    S <- ncol(mu_h)
    stopifnot(length(alpha) == K, length(sigma_u) == K)
    
    idx_H <- match(pure_heads, seq_levels)
    idx_T <- match(pure_tails, seq_levels)
    if (is.na(idx_H) || is.na(idx_T)) {
      stop(
        "EX1.1 participants: anchor labels not found in rq4_seq_levels for tr='", tr, "'.\n",
        "pure_heads=", pure_heads, " pure_tails=", pure_tails
      )
    }
    
    theta_H <- mu_h[, idx_H]
    theta_T <- mu_h[, idx_T]
    
    # ------------------------------------------------------------
    # Baseline theta0^(k) = E_{u~N(0,sigma_u^(k))}[ logistic(alpha^(k)+u) ]
    # (Monte Carlo; purely numerical approximation parameter = base_mc)
    # Use common random numbers across k for stability.
    # ------------------------------------------------------------
    set.seed(seed + 881L)
    z_std <- rnorm(base_mc)  # common across k
    
    # theta0[k] = mean( logistic(alpha[k] + sigma_u[k] * z_std) )
    theta0 <- vapply(seq_len(K), function(k) {
      mean(inv_logit(alpha[k] + sigma_u[k] * z_std))
    }, numeric(1))
    
    stopifnot(all(is.finite(theta0)), all(theta0 > 0), all(theta0 < 1))
    
    # ------------------------------------------------------------
    # Betting trials only (align with RQ4 / EX1.1 definition)
    # ------------------------------------------------------------
    d <- dt[treat == tr & is.finite(stake) & stake > 0]
    d <- d[side %in% c(lab_heads, lab_tails)]
    
    if (nrow(d) == 0) {
      warning("EX1.1 participants: no betting trials after filtering for tr='", tr, "'. Skipping.")
      next
    }
    
    d[, h := as.integer(side == lab_heads)]
    d[, sign := 2L * h - 1L]  # Heads=+1, Tails=-1
    
    # Map sequences into [1..S]
    d[, sid := match(seq, seq_levels)]
    if (anyNA(d$sid)) stop("EX1.1 participants: some seq not found in rq4_seq_levels for tr='", tr, "'.")
    
    # participant index
    pid_levels <- sort(unique(d$pid))
    d[, pid_i := match(pid, pid_levels)]
    Np <- length(pid_levels)
    
    # ------------------------------------------------------------
    # Similarities and weights (draw-by-draw, sequence-by-sequence)
    # d_a(k,s) = 1{|theta_s^(k) - theta_a^(k)| < eps}
    # w_a(k,s) = normalize(d_a(k,s) + eta)
    # ------------------------------------------------------------
    # These are K x S logical matrices
    dH <- abs(mu_h - matrix(theta_H, nrow = K, ncol = S)) < eps
    dT <- abs(mu_h - matrix(theta_T, nrow = K, ncol = S)) < eps
    d0 <- abs(mu_h - matrix(theta0,  nrow = K, ncol = S)) < eps
    
    denom <- (dH + dT + d0) + 3 * eta
    wH <- (dH + eta) / denom
    wT <- (dT + eta) / denom
    diff_w <- wH - wT  # K x S
    
    # ------------------------------------------------------------
    # chi_draws[k,i] = mean_s z_is^(k) over betting trials
    # and z_is^(k) = sign(h_is) * diff_w[k, sid]
    # ------------------------------------------------------------
    chi_draws <- matrix(NA_real_, nrow = K, ncol = Np)
    
    for (j in seq_len(Np)) {
      idx <- which(d$pid_i == j)
      if (length(idx) == 0) next
      
      sid_j  <- d$sid[idx]
      sign_j <- d$sign[idx]  # length = n_trials_i
      
      # vectorized over trials for each k
      # chi[k] = mean( sign_j * diff_w[k, sid_j] )
      chi_draws[, j] <- vapply(seq_len(K), function(k) {
        mean(sign_j * diff_w[k, sid_j])
      }, numeric(1))
    }
    
    if (!all(is.finite(chi_draws))) {
      stop("EX1.1 participants: non-finite chi_draws produced for tr='", tr, "'.")
    }
    
    # ------------------------------------------------------------
    # Summaries + classification (prereg cutpoints)
    # ------------------------------------------------------------
    tbl <- data.table(
      pid        = pid_levels,
      chi_median = apply(chi_draws, 2, median),
      chi_mean   = apply(chi_draws, 2, mean),
      chi_q025   = apply(chi_draws, 2, quantile, probs = 0.025),
      chi_q975   = apply(chi_draws, 2, quantile, probs = 0.975)
    )
    
    tbl[, HH   := apply(chi_draws, 2, function(x) mean(x >  delta))]
    tbl[, G    := apply(chi_draws, 2, function(x) mean(x < -delta))]
    tbl[, Pabs := apply(chi_draws, 2, function(x) mean(abs(x) > delta))]
    
    tbl[, class :=
          fifelse(HH >= 0.95, "solid_hot",
                  fifelse(HH >= 0.90, "likely_hot",
                          fifelse(HH >= 0.75, "leaning_hot",
                                  fifelse(G  >= 0.95, "solid_gambler",
                                          fifelse(G  >= 0.90, "likely_gambler",
                                                  fifelse(G  >= 0.75, "leaning_gambler",
                                                          "neutral"))))))]
    
    setorder(tbl, pid)
    
    # ------------------------------------------------------------
    # Save CSV + RDS (RDS is required for EX1.2)
    # ------------------------------------------------------------
    f_csv <- file.path(out_dir, paste0("ex1_1_participants_", tr, ".csv"))
    fwrite(tbl, f_csv)
    msg("Saved: ", f_csv)
    
    f_rds <- file.path(mod_dir, paste0("ex1_1_participants_", tr, ".rds"))
    saveRDS(list(
      dataset = ds,
      treatment = tr,
      pid_levels = pid_levels,
      K = K,
      delta = delta,
      eps = eps,
      eta = eta,
      base_mc = base_mc,
      chi_draws = chi_draws,
      table = tbl
    ), f_rds)
    msg("Saved: ", f_rds)
    
    outputs[[tr]] <- list(csv = f_csv, rds = f_rds, table = tbl)
  }
  
  invisible(outputs)
}

# Example:
# ex1_1_participants(cfg)