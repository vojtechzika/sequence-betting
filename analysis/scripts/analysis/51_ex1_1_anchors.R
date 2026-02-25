# ============================================================
# scripts/analysis/51_ex1_1_anchors.R
#
# EX1.I.1 Anchors + distributional similarity (NO sequence summaries)
#
# Implements prereg Section "Distributional Similarity and the GHI" up to:
#  - theta_s^(k) from RQ4 posterior draws (mu_h)
#  - anchors: HHHHHH, TTTTTT
#  - baseline theta0^(k) = E_{u~N(0,sigma_u^(k))}[ logistic(alpha^(k) + u) ]
#  - similarity: d_a(s) = P(|theta_s - theta_a| < eps | data)
#  - weights: w_a(s) = (d_a(s)+eta) / (d_H + d_T + d_0 + 3eta)
#
# Outputs per treatment:
#  - data/clean/<ds>/output/ex1_1_anchors_<tr>.csv
#  - data/clean/<ds>/models/ex1_1_anchors_<tr>.rds   (for downstream scripts)
#
# Notes:
# - Does NOT compute trial-level z_is, chi_s, chi_i, or any labels.
# ============================================================

library(data.table)
library(rstan)

ex1_1_anchors <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Prereg constants
  # ----------------------------
  eta <- 1e-6
  
  eps_main <- 0.05
  eps_vec  <- c(0.05, 0.03, 0.08)
  
  # MC size for logistic-normal expectation in theta0^(k)
  # (not specified in prereg; purely numerical approximation parameter)
  L_mc <- cfg$run$ex1_L_mc
  if (is.null(L_mc)) L_mc <- 4000L
  L_mc <- as.integer(L_mc)
  stopifnot(length(L_mc) == 1L, L_mc >= 500L)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # theta0^(k) = E_u logistic(alpha^(k) + u), u~N(0, sigma_u^(k))
  # computed by MC with COMMON random numbers across k for stability
  theta0_mc <- function(alpha, sigma_u, z_std) {
    # alpha: scalar, sigma_u: scalar, z_std: length L_mc (standard normal)
    mean(inv_logit(alpha + sigma_u * z_std))
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    
    if (!file.exists(f_fit)) stop("EX1.1: Missing RQ4 fit: ", f_fit)
    if (!file.exists(f_seq)) stop("EX1.1: Missing RQ4 seq levels: ", f_seq)
    
    fit <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    stopifnot(length(seq_levels) >= 1L)
    
    post <- rstan::extract(fit)
    
    # We take theta_s^(k) from generated quantity mu_h in RQ4 fit:
    # theta_s^(k) = mu_s^{h(k)} = E_i[pi_is^{h(k)}]
    if (is.null(post$mu_h)) stop("EX1.1: RQ4 fit missing generated quantity mu_h.")
    mu_h <- post$mu_h  # K x S
    
    if (is.null(post$alpha))   stop("EX1.1: RQ4 fit missing alpha.")
    if (is.null(post$sigma_u)) stop("EX1.1: RQ4 fit missing sigma_u.")
    
    alpha   <- as.numeric(post$alpha)     # K
    sigma_u <- as.numeric(post$sigma_u)   # K
    
    K <- length(alpha)
    S <- length(seq_levels)
    stopifnot(length(sigma_u) == K)
    stopifnot(nrow(mu_h) == K, ncol(mu_h) == S)
    
    thetaS <- mu_h
    
    # ------------------------------------------------------------
    # Anchors from design config (NO hardcoding)
    # ------------------------------------------------------------
    stopifnot(!is.null(cfg$design$seq$anchor_labels))
    stopifnot(is.list(cfg$design$seq$anchor_labels))
    
    lab_H <- as.character(cfg$design$seq$anchor_labels$pure_heads)
    lab_T <- as.character(cfg$design$seq$anchor_labels$pure_tails)
    
    stopifnot(length(lab_H) == 1L, nzchar(lab_H))
    stopifnot(length(lab_T) == 1L, nzchar(lab_T))
    stopifnot(lab_H != lab_T)
    
    # Validate anchors exist in seq_levels
    sH <- match(lab_H, seq_levels)
    sT <- match(lab_T, seq_levels)
    
    if (is.na(sH) || is.na(sT)) {
      stop(
        "EX1.1: Anchor sequences defined in design$seq$anchor_labels ",
        "not found in rq4_seq_levels for treatment='", tr, "'.\n",
        "Expected anchors: ", lab_H, " and ", lab_T
      )
    }
    
    thetaH <- thetaS[, sH]
    thetaT <- thetaS[, sT]
    
    # Baseline theta0^(k) as per prereg:
    # theta0^(k) = E_{u~N(0,sigma_u^(k))}[ logistic(alpha^(k) + u) ]
    set.seed(seed)
    z_std <- rnorm(L_mc)
    
    theta0 <- vapply(seq_len(K), function(k) {
      theta0_mc(alpha = alpha[k], sigma_u = sigma_u[k], z_std = z_std)
    }, numeric(1))
    
    stopifnot(all(is.finite(theta0)), all(theta0 > 0), all(theta0 < 1))
    
    # ----------------------------
    # Similarities d_a(s) and weights w_a(s)
    #   d_a(s) = (1/K) sum_k 1{|theta_s^(k) - theta_a^(k)| < eps}
    #   w_a(s) = (d_a(s) + eta) / (d_H + d_T + d_0 + 3eta)
    # ----------------------------
    tbl <- data.table(sequence = seq_levels)
    
    for (eps in eps_vec) {
      
      stopifnot(is.finite(eps), eps > 0)
      
      # Indicators K x S
      I_H <- abs(thetaS - matrix(thetaH, nrow = K, ncol = S)) < eps
      I_T <- abs(thetaS - matrix(thetaT, nrow = K, ncol = S)) < eps
      I_0 <- abs(thetaS - matrix(theta0, nrow = K, ncol = S)) < eps
      
      dH <- colMeans(I_H)
      dT <- colMeans(I_T)
      d0 <- colMeans(I_0)
      
      denom <- dH + dT + d0 + 3 * eta
      
      wH <- (dH + eta) / denom
      wT <- (dT + eta) / denom
      w0 <- (d0 + eta) / denom
      
      # sanity: weights sum to 1 (up to numerical tolerance)
      if (max(abs((wH + wT + w0) - 1)) > 1e-10) {
        stop("EX1.1: weight normalization failed numerically (eps=", eps, ", tr=", tr, ").")
      }
      
      suf <- gsub("\\.", "", sprintf("%.2f", eps))
      
      tbl[, (paste0("dH_eps_", suf)) := dH]
      tbl[, (paste0("dT_eps_", suf)) := dT]
      tbl[, (paste0("d0_eps_", suf)) := d0]
      
      tbl[, (paste0("wH_eps_", suf)) := wH]
      tbl[, (paste0("wT_eps_", suf)) := wT]
      tbl[, (paste0("w0_eps_", suf)) := w0]
    }
    
    # Keep a clear pointer to which eps is "main"
    suf_main <- gsub("\\.", "", sprintf("%.2f", eps_main))
    tbl[, eps_main := eps_main]
    tbl[, wH_main := get(paste0("wH_eps_", suf_main))]
    tbl[, wT_main := get(paste0("wT_eps_", suf_main))]
    tbl[, w0_main := get(paste0("w0_eps_", suf_main))]
    
    setorder(tbl, sequence)
    
    # Save CSV for inspection
    f_csv <- file.path(out_dir, paste0("ex1_1_anchors_", tr, ".csv"))
    fwrite(tbl, f_csv)
    msg("Saved: ", f_csv)
    
    # Save RDS for downstream scripts (stores draws + anchors + baseline)
    out_rds <- list(
      dataset = ds,
      treatment = tr,
      seq_levels = seq_levels,
      K = K,
      eps_vec = eps_vec,
      eps_main = eps_main,
      eta = eta,
      L_mc = L_mc,
      theta_draws = list(
        thetaS = thetaS,    # K x S  (mu_h draws)
        thetaH = thetaH,    # K
        thetaT = thetaT,    # K
        theta0 = theta0     # K
      ),
      table = tbl
    )
    
    f_rds <- file.path(mod_dir, paste0("ex1_1_anchors_", tr, ".rds"))
    saveRDS(out_rds, f_rds)
    msg("Saved: ", f_rds)
    
    outputs[[tr]] <- out_rds
  }
  
  invisible(outputs)
}

# Example:
# ex1_1_anchors(cfg)