# ============================================================
# scripts/analysis/51_ex1_1_anchors.R
#
# EX1.I.1 Anchors + distributional similarity (NO sequence summaries)
#
# Implements prereg "Distributional Similarity and the GHI" up to:
#  - theta_s^(k) from RQ4 posterior draws (mu_h)
#  - anchors: pure_heads, pure_tails from cfg$design$seq$anchor_labels
#  - baseline theta0^(k) from RQ4 posterior draws (hbar)
#  - similarity: d_a(s) = P(|theta_s - theta_a| < eps | data)
#  - weights: w_a(s) = (d_a(s)+eta) / (d_H + d_T + d_0 + 3eta)
#
# Outputs per treatment:
#  - data/clean/<ds>/output/ex1_1_anchors_<tr>.csv
#  - data/clean/<ds>/models/ex1_1_anchors_<tr>.rds
#
# Notes:
# - Does NOT compute trial-level z_is, chi_s, chi_i, or any labels.
# - Uses ONLY Stan-generated hbar for theta0.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex1_1_anchors <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)  # kept for reproducibility even if not used here
  
  tr_vec <- cfg$run$treatment
  
  # ----------------------------
  # Prereg constants
  # ----------------------------
  eta <- 1e-6
  
  eps_main <- 0.05
  eps_vec  <- c(0.05, 0.03, 0.08)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    
    fit <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    
    post <- rstan::extract(fit)
    
    # theta_s^(k) from generated quantity mu_h (population-based)
    thetaS <- post$mu_h  # K x S
    
    # baseline theta0^(k) MUST come from RQ4 generated quantity hbar
    if (is.null(post$hbar)) {
      stop(
        "EX1.1: RQ4 fit missing generated quantity hbar.\n",
        "This script uses ONLY Stan-generated hbar for theta0; regenerate RQ4 fits with the updated rq4_side.stan."
      )
    }
    theta0 <- as.numeric(post$hbar)  # K
    
    K <- length(theta0)
    S <- length(seq_levels)
    
    # ------------------------------------------------------------
    # Anchors from design config (NO hardcoding)
    # ------------------------------------------------------------
    lab_H <- as.character(cfg$design$seq$anchor_labels$pure_heads)
    lab_T <- as.character(cfg$design$seq$anchor_labels$pure_tails)
    
    sH <- match(lab_H, seq_levels)
    sT <- match(lab_T, seq_levels)
    
    thetaH <- thetaS[, sH]
    thetaT <- thetaS[, sT]
    
    # ----------------------------
    # Similarities d_a(s) and weights w_a(s)
    # ----------------------------
    tbl <- data.table(sequence = seq_levels)
    
    for (eps in eps_vec) {
      
      I_H <- abs(sweep(thetaS, 1, thetaH, "-")) < eps
      I_T <- abs(sweep(thetaS, 1, thetaT, "-")) < eps
      I_0 <- abs(sweep(thetaS, 1, theta0, "-")) < eps
      
      dH <- colMeans(I_H)
      dT <- colMeans(I_T)
      d0 <- colMeans(I_0)
      
      denom <- dH + dT + d0 + 3 * eta
      
      wH <- (dH + eta) / denom
      wT <- (dT + eta) / denom
      w0 <- (d0 + eta) / denom
      
      suf <- gsub("\\.", "", sprintf("%.2f", eps))
      
      tbl[, (paste0("dH_eps_", suf)) := dH]
      tbl[, (paste0("dT_eps_", suf)) := dT]
      tbl[, (paste0("d0_eps_", suf)) := d0]
      
      tbl[, (paste0("wH_eps_", suf)) := wH]
      tbl[, (paste0("wT_eps_", suf)) := wT]
      tbl[, (paste0("w0_eps_", suf)) := w0]
    }
    
    # main eps pointers
    suf_main <- gsub("\\.", "", sprintf("%.2f", eps_main))
    tbl[, eps_main := eps_main]
    tbl[, wH_main := get(paste0("wH_eps_", suf_main))]
    tbl[, wT_main := get(paste0("wT_eps_", suf_main))]
    tbl[, w0_main := get(paste0("w0_eps_", suf_main))]
    
    setorder(tbl, sequence)
    
    # Save CSV
    f_csv <- file.path(out_dir, paste0("ex1_1_anchors_", tr, ".csv"))
    fwrite(tbl, f_csv)
    msg("Saved: ", f_csv)
    
    # Save RDS for downstream
    out_rds <- list(
      dataset = ds,
      treatment = tr,
      seq_levels = seq_levels,
      K = K,
      eps_vec = eps_vec,
      eps_main = eps_main,
      eta = eta,
      theta_draws = list(
        thetaS = thetaS,   # K x S (mu_h draws)
        thetaH = thetaH,   # K
        thetaT = thetaT,   # K
        theta0 = theta0    # K (hbar draws)
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