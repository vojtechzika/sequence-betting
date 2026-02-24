# ============================================================
# scripts/analysis/42_rq4_sequences.R
# Sequence table from RQ4 fit (side choices conditional on betting)
# Writes: data/clean/<ds>/output/rq4_sequences_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

rq4_sequences <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  
  # delta grid: main + sensitivities
  delta_vec <- c(0.05, 0.03, 0.08)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    if (!file.exists(f_fit)) stop("RQ4 fit not found: ", f_fit)
    if (!file.exists(f_seq)) stop("RQ4 seq_levels not found: ", f_seq)
    
    fit <- readRDS(f_fit)
    seq_levels <- readRDS(f_seq)
    
    post <- rstan::extract(fit)
    if (is.null(post$mu_h)) stop("RQ4 fit missing generated quantity mu_h.")
    if (is.null(post$hbar)) stop("RQ4 fit missing generated quantity hbar.")
    
    mu_draws <- post$mu_h  # iters x S
    hbar_draws <- post$hbar  # iters
    
    tbl <- data.table(
      sequence    = as.character(seq_levels),
      mu_h_median = apply(mu_draws, 2, median),
      mu_h_mean   = apply(mu_draws, 2, mean),
      mu_h_q025   = apply(mu_draws, 2, quantile, probs = 0.025),
      mu_h_q975   = apply(mu_draws, 2, quantile, probs = 0.975)
    )
    
    # baseline summaries
    tbl[, hbar_median := median(hbar_draws)]
    tbl[, hbar_mean   := mean(hbar_draws)]
    tbl[, hbar_q025   := quantile(hbar_draws, 0.025)]
    tbl[, hbar_q975   := quantile(hbar_draws, 0.975)]
    
    # H_s(delta) and T_s(delta)
    for (delta in delta_vec) {
      nmH <- paste0("H_delta_", gsub("\\.", "", sprintf("%.2f", delta)))
      nmT <- paste0("T_delta_", gsub("\\.", "", sprintf("%.2f", delta)))
      
      # per-sequence posterior prob against per-iter baseline
      tbl[, (nmH) := apply(mu_draws, 2, function(x) mean(x > (hbar_draws + delta)))]
      tbl[, (nmT) := apply(mu_draws, 2, function(x) mean(x < (hbar_draws - delta)))]
    }
    
    # 7-point directional label (main delta = 0.05) using 0.95/0.80 cutpoints
    dmain <- 0.05
    colH <- paste0("H_delta_", gsub("\\.", "", sprintf("%.2f", dmain)))
    colT <- paste0("T_delta_", gsub("\\.", "", sprintf("%.2f", dmain)))
    
    tbl[, direction_label :=
          fifelse(get(colT) >= 0.95, "strong_tail",
                  fifelse(get(colT) >= 0.80, "moderate_tail",
                          fifelse(get(colT) >= 0.50, "weak_tail",
                                  fifelse(get(colH) >= 0.50, "weak_head",
                                          fifelse(get(colH) >= 0.80, "moderate_head",
                                                  fifelse(get(colH) >= 0.95, "strong_head", "neutral"))))))]
    
    setorder(tbl, sequence)
    
    f_out <- file.path(out_dir, paste0("rq4_sequences_", tr, ".csv"))
    fwrite(tbl, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}