# ============================================================
# scripts/analysis/12_rq1_sequences.R
#   - post-processes RQ1 sequence summaries per treatment
#   - reads ONLY: fit + saved levels
#   - writes ONLY: rq1_sequences_<ds>_<tr>.csv
#   - safeguard: if output CSV exists, skip that treatment
# ============================================================

library(data.table)
library(rstan)

rq1_sequences <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  run    <- cfg$run
  design <- cfg$design
  
  ds <- as.character(run$dataset)
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Rho thresholds (from design)
  # ----------------------------
  stopifnot(!is.null(design$rhos), !is.null(design$rhos$rq1_rho))
  rho_vec <- as.numeric(design$rhos$rq1_rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  
  # ----------------------------
  # Input data (for n_trials)
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
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, y := as.integer(!is.na(stake) & stake > 0)]
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit        <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq1_seq_levels_", tr, ".rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid_levels) || !file.exists(f_seq_levels)) {
      stop(
        "Missing RQ1 Stan artifacts for dataset='", ds, "', treatment='", tr, "'.\n",
        "Expected:\n",
        "  ", f_fit, "\n",
        "  ", f_pid_levels, "\n",
        "  ", f_seq_levels, "\n",
        "Run rq1_stan(cfg) first."
      )
    }
    
    f_out_csv <- file.path(out_dir, paste0("rq1_sequences_", tr, ".csv"))
    if (file.exists(f_out_csv)) {
      warning(
        "RQ1 sequence table already exists for dataset='", ds, "', treatment='", tr, "'.\n",
        "Post-processing was NOT executed.\n",
        "To rerun, delete:\n  ", f_out_csv
      )
      next
    }
    
    # subset for n_trials
    d <- dt[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ1 sequences: No rows for treat='", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    
    pid_levels <- readRDS(f_pid_levels)
    seq_levels <- readRDS(f_seq_levels)
    
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    mu_draws <- post$mu_b  # iterations x S
    
    # n_trials per sequence (from observed trials)
    n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
    setkey(n_trials_by_seq, seq)
    
    seq_tbl <- data.table(
      sequence    = seq_levels,
      mu_b_median = apply(mu_draws, 2, median),
      mu_b_mean   = apply(mu_draws, 2, mean),
      mu_b_q025   = apply(mu_draws, 2, quantile, probs = 0.025),
      mu_b_q975   = apply(mu_draws, 2, quantile, probs = 0.975)
    )
    
    seq_tbl[, n_trials := n_trials_by_seq[.(sequence), n_trials]]
    seq_tbl[, n_trials := as.integer(n_trials)]
    
    for (rho in rho_vec) {
      nm <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      seq_tbl[, (nm) := apply(mu_draws, 2, function(x) mean(x < (1 - rho)))]
    }
    
    rho_main <- rho_vec[1]
    col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    seq_tbl[, underbet_label :=
              fifelse(get(col_main) >= 0.95, "strong",
                      fifelse(get(col_main) >= 0.80, "moderate",
                              fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
    
    setorder(seq_tbl, sequence)
    
    fwrite(seq_tbl, f_out_csv)
    msg("Saved:", f_out_csv)
    
    outputs[[tr]] <- list(
      sequences = seq_tbl,
      out_csv = f_out_csv,
      fit_file = f_fit,
      pid_levels_file = f_pid_levels,
      seq_levels_file = f_seq_levels
    )
  }
  
  invisible(outputs)
}