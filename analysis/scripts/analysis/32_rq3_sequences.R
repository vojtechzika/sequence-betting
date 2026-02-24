# scripts/analysis/32_rq3_sequences.R
# Sequence table from rq3 fit.

library(data.table)
library(rstan)

rq3_sequences <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  stopifnot(!is.null(cfg$design$rhos), !is.null(cfg$design$rhos$rq3_rho))
  rho_vec <- as.numeric(cfg$design$rhos$rq3_rho)
  stopifnot(length(rho_vec) >= 1L,
            all(is.finite(rho_vec)),
            all(rho_vec > 0),
            all(rho_vec < 1))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq3_seq_levels_", tr, ".rds"))
    
    if (!file.exists(f_fit)) stop("RQ3 fit not found: ", f_fit)
    if (!file.exists(f_seq)) stop("RQ3 seq_levels not found: ", f_seq)
    
    fit <- readRDS(f_fit)
    seq_levels <- as.character(readRDS(f_seq))
    
    post <- rstan::extract(fit)
    if (is.null(post$mu_c)) stop("RQ3 fit does not contain generated quantity mu_c.")
    
    mu_draws <- post$mu_c  # iters x S
    
    stopifnot(is.matrix(mu_draws))
    stopifnot(ncol(mu_draws) == length(seq_levels))
    
    tbl <- data.table(
      sequence     = seq_levels,
      mu_c_median  = apply(mu_draws, 2, median),
      mu_c_mean    = apply(mu_draws, 2, mean),
      mu_c_q025    = apply(mu_draws, 2, quantile, probs = 0.025),
      mu_c_q975    = apply(mu_draws, 2, quantile, probs = 0.975)
    )
    
    # L_s(rho) = P(mu_s^c > rho)
    for (rho in rho_vec) {
      nm <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      tbl[, (nm) := apply(mu_draws, 2, function(x) mean(x > rho))]
    }
    
    # Main rho (first element) for categorical label
    rho_main <- rho_vec[1]
    colL <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    
    tbl[, loss_label :=
          fifelse(get(colL) >= 0.95, "strong",
                  fifelse(get(colL) >= 0.80, "moderate",
                          fifelse(get(colL) >= 0.50, "weak", "neutral")))]
    
    # Deterministic ordering aligned with seq_levels
    tbl[, sequence := factor(sequence, levels = seq_levels)]
    setorder(tbl, sequence)
    tbl[, sequence := as.character(sequence)]
    
    f_out <- file.path(out_dir, paste0("rq3_sequences_", tr, ".csv"))
    fwrite(tbl, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}