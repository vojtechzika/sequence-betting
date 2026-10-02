# ============================================================
# 51_ex1_1_anchors.R
#
# EX1.I  Collapsed anchor similarity table
#
# Computes the scalar posterior similarity probabilities d_a(s)
# and normalized weights w_a(s) for display and the EX4 crosswalk.
# Draw-level propagation (chi_s, chi_i) is handled in 53_.
#
# OUTPUTS (per treatment)
#   path_out/ex1_1_anchors_<tr>.csv
#   path_mod/ex1_1_anchors_<tr>.rds
# ============================================================

ex1_1_anchors <- function(cfg) {
  
  eps_vec  <- as.numeric(cfg$design$ex1$eps)
  eps_main <- eps_vec[1]
  eta      <- as.numeric(cfg$design$ex1$eta)
  
  stopifnot(
    length(eps_vec) >= 1L, all(is.finite(eps_vec)), all(eps_vec > 0),
    length(eta) == 1L, is.finite(eta), eta > 0
  )
  
  lab_H <- as.character(cfg$design$seq$anchor_labels$pure_heads)
  lab_T <- as.character(cfg$design$seq$anchor_labels$pure_tails)
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(path_mod, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_seq <- file.path(path_mod, paste0("rq4_seq_levels_",    tr, "_full.rds"))
    
    fit        <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    
    post <- rstan::extract(fit)
    
    if (is.null(post$mu_h) || is.null(post$hbar)) {
      stop(
        "EX1.1: RQ4 fit missing 'mu_h' or 'hbar' for tr='", tr,
        "'.\nRegenerate RQ4 fits with the updated rq4_side.stan."
      )
    }
    
    thetaS <- post$mu_h              # K x S
    theta0 <- as.numeric(post$hbar)  # K
    
    sH <- match(lab_H, seq_levels)
    sT <- match(lab_T, seq_levels)
    
    if (is.na(sH)) stop("EX1.1: anchor '", lab_H, "' not in seq_levels for tr='", tr, "'.")
    if (is.na(sT)) stop("EX1.1: anchor '", lab_T, "' not in seq_levels for tr='", tr, "'.")
    
    thetaH <- thetaS[, sH]  # K
    thetaT <- thetaS[, sT]  # K
    
    # ---- Collapsed similarity probabilities and weights ----
    tbl <- data.table(sequence = seq_levels)
    
    for (eps in eps_vec) {
      
      I_H <- abs(sweep(thetaS, 1, thetaH, "-")) < eps
      I_T <- abs(sweep(thetaS, 1, thetaT, "-")) < eps
      I_0 <- abs(sweep(thetaS, 1, theta0,  "-")) < eps
      
      dH <- colMeans(I_H)
      dT <- colMeans(I_T)
      d0 <- colMeans(I_0)
      
      denom <- dH + dT + d0 + 3 * eta
      
      suf <- gsub("\\.", "", sprintf("%.2f", eps))
      tbl[, (paste0("dH_eps_", suf)) := dH]
      tbl[, (paste0("dT_eps_", suf)) := dT]
      tbl[, (paste0("d0_eps_", suf)) := d0]
      tbl[, (paste0("wH_eps_", suf)) := (dH + eta) / denom]
      tbl[, (paste0("wT_eps_", suf)) := (dT + eta) / denom]
      tbl[, (paste0("w0_eps_", suf)) := (d0 + eta) / denom]
    }
    
    suf_main <- gsub("\\.", "", sprintf("%.2f", eps_main))
    tbl[, eps_main := eps_main]
    tbl[, wH_main := get(paste0("wH_eps_", suf_main))]
    tbl[, wT_main := get(paste0("wT_eps_", suf_main))]
    tbl[, w0_main := get(paste0("w0_eps_", suf_main))]
    setorder(tbl, sequence)
    
    # ---- Outputs ----
    f_csv <- file.path(path_out, paste0("ex1_1_anchors_", tr, ".csv"))
    f_rds <- file.path(path_mod, paste0("ex1_1_anchors_", tr, ".rds"))
    
    if (!should_skip(f_csv, cfg, "output", paste0("EX1.1 anchors (", tr, ")"))) {
      fwrite(tbl, f_csv)
      msg("Saved: ", f_csv)
    }
    
    out_rds <- list(
      treatment  = tr,
      seq_levels = seq_levels,
      eps_vec    = eps_vec,
      eps_main   = eps_main,
      eta        = eta,
      theta_draws = list(
        thetaS = thetaS,
        thetaH = thetaH,
        thetaT = thetaT,
        theta0 = theta0
      ),
      table = tbl
    )
    
    if (!should_skip(f_rds, cfg, "model", paste0("EX1.1 anchors RDS (", tr, ")"))) {
      saveRDS(out_rds, f_rds)
      msg("Saved: ", f_rds)
    }
    
    outputs[[tr]] <- out_rds
  }
  
  invisible(outputs)
}