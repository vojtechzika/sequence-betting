# ============================================================
# scripts/analysis/14_rq1_exploratory_tables.R
#   RQ1: Post-process FULL-sample outputs per treatment
#   - Sequence-level summaries: mu_b (from Stan generated quantities)
#   - Participant-level summaries: mu_b_i (from Stan generated quantities)
#
# Reads ONLY:
#   models/rq1_fit_sequences_<tr>.rds
#   models/rq1_pid_levels_<tr>.rds
#   models/rq1_seq_levels_<tr>.rds
#   clean/<ds>/master_sequences.csv (for n_trials)
#
# Writes ONLY:
#   output/rq1_sequences_exploratory_<tr>.csv
#   output/rq1_participants_exploratory_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

rq1_exploratory_tables <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$plan))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Rho thresholds (from design)
  # ----------------------------
  stopifnot(!is.null(cfg$design$rhos), !is.null(cfg$design$rhos$rq1_rho))
  rho_vec <- as.numeric(cfg$design$rhos$rq1_rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  rho_main <- rho_vec[1]
  
  # ----------------------------
  # Load master (for n_trials)
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  required <- c("pid", "treat", "stake", "seq")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) stop("master_sequences.csv missing columns: ", paste(missing, collapse = ", "))
  
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[is.na(stake), stake := 0]
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, ".rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) {
      stop(
        "Missing RQ1 Stan artifacts for dataset='", ds, "', treatment='", tr, "'.\n",
        "Expected:\n  ", f_fit, "\n  ", f_pid, "\n  ", f_seq, "\n",
        "Run rq1_stan(cfg) first."
      )
    }
    
    d <- dt[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ1 post: No rows for treat='", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    
    # levels
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    stopifnot(length(seq_levels) == 64L)
    
    # fit + posterior
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    # STRICT: no fallback
    if (is.null(post$mu_b))   stop("RQ1 post: fit missing 'mu_b' (rerun RQ1 Stan).")
    if (is.null(post$mu_b_i)) stop("RQ1 post: fit missing 'mu_b_i' (rerun RQ1 Stan).")
    
    mu_b_draws  <- post$mu_b    # iters x S
    mu_bi_draws <- post$mu_b_i  # iters x N
    
    stopifnot(is.matrix(mu_b_draws),  ncol(mu_b_draws)  == length(seq_levels))
    stopifnot(is.matrix(mu_bi_draws), ncol(mu_bi_draws) == length(pid_levels))
    
    # ----------------------------
    # SEQUENCES table
    # ----------------------------
    f_seq_csv <- file.path(out_dir, paste0("rq1_sequences_exploratory_", tr, ".csv"))
    if (should_skip(
      paths = f_seq_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ1 sequences exploratory (", ds, "/", tr, ")")
    )) {
      seq_tbl <- NULL
    } else {
      
      n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
      setkey(n_trials_by_seq, seq)
      
      seq_tbl <- data.table(
        sequence    = seq_levels,
        mu_b_median = apply(mu_b_draws, 2, median),
        mu_b_mean   = apply(mu_b_draws, 2, mean),
        mu_b_q025   = apply(mu_b_draws, 2, quantile, probs = 0.025),
        mu_b_q975   = apply(mu_b_draws, 2, quantile, probs = 0.975)
      )
      
      seq_tbl[, n_trials := n_trials_by_seq[.(sequence), n_trials]]
      seq_tbl[, n_trials := as.integer(n_trials)]
      
      for (rho in rho_vec) {
        nm <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        thr <- 1 - rho
        seq_tbl[, (nm) := apply(mu_b_draws, 2, function(x) mean(x < thr))]
      }
      
      col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
      seq_tbl[, underbet_label :=
                fifelse(get(col_main) >= 0.95, "strong",
                        fifelse(get(col_main) >= 0.80, "moderate",
                                fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
      
      setorder(seq_tbl, sequence)
      fwrite(seq_tbl, f_seq_csv)
      msg("Saved: ", f_seq_csv)
    }
    
    # ----------------------------
    # PARTICIPANTS table
    # ----------------------------
    f_pid_csv <- file.path(out_dir, paste0("rq1_participants_exploratory_", tr, ".csv"))
    if (should_skip(
      paths = f_pid_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ1 participants exploratory (", ds, "/", tr, ")")
    )) {
      pid_tbl <- NULL
    } else {
      
      n_trials_by_pid <- d[, .(n_trials = .N), by = pid]
      setkey(n_trials_by_pid, pid)
      
      pid_tbl <- data.table(
        pid         = pid_levels,
        mu_i_median = apply(mu_bi_draws, 2, median),
        mu_i_mean   = apply(mu_bi_draws, 2, mean),
        mu_i_q025   = apply(mu_bi_draws, 2, quantile, probs = 0.025),
        mu_i_q975   = apply(mu_bi_draws, 2, quantile, probs = 0.975)
      )
      
      for (rho in rho_vec) {
        nm <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        thr <- 1 - rho
        pid_tbl[, (nm) := apply(mu_bi_draws, 2, function(x) mean(x < thr))]
      }
      
      col_main <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
      pid_tbl[, underbet_label :=
                fifelse(get(col_main) >= 0.95, "solid",
                        fifelse(get(col_main) >= 0.90, "likely",
                                fifelse(get(col_main) >= 0.75, "leaning", "neutral")))]
      
      pid_tbl[, n_trials := n_trials_by_pid[.(pid), n_trials]]
      pid_tbl[, n_trials := as.integer(n_trials)]
      
      setorder(pid_tbl, pid)
      fwrite(pid_tbl, f_pid_csv)
      msg("Saved: ", f_pid_csv)
    }
    
    outputs[[tr]] <- list(
      sequences_csv = f_seq_csv,
      participants_csv = f_pid_csv,
      sequences = seq_tbl,
      participants = pid_tbl
    )
  }
  
  invisible(outputs)
}

# Example:
# source(here::here("scripts","analysis","14_rq1_exploratory_tables.R"))
# rq1_exploratory_tables(cfg)