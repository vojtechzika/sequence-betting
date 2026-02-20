# scripts/analysis/21_rq1_sequences.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq1_sequences <- function(dataset = "pilot",
                          treat_fn = "m25",
                          rho_vec = c(0.10, 0.08, 0.12),
                          seed = 12345) {
  
  # ----------------------------
  # Files
  # ----------------------------
  infile <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- here::here("stan", "rq1_bets.stan")
  stopifnot(file.exists(stan_file))
  
  # Normalize line endings (prevents hash-mismatch churn)
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # ----------------------------
  # Load + schema checks
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "stake", "seq")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      "master_sequences.csv missing columns: ", paste(missing, collapse = ", "),
      "\nExpected at least: pid, treat, stake, seq"
    )
  }
  
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  
  # ----------------------------
  # Filter to FN (confirmatory)
  # ----------------------------
  d <- dt[treat == treat_fn]
  if (nrow(d) == 0) stop("No rows after filtering treat == '", treat_fn, "'")
  
  # bet indicator
  d[, y := as.integer(stake > 0)]
  
  # ----------------------------
  # Map ids
  # ----------------------------
  pid_levels <- sort(unique(d$pid))
  seq_levels <- sort(unique(d$seq))  # alphabetical
  
  if (length(seq_levels) != 64) {
    warning(
      "RQ1: expected 64 sequences but found ", length(seq_levels),
      ". Proceeding anyway."
    )
  }
  
  d[, pid_i := match(pid, pid_levels)]
  d[, sid_s := match(seq, seq_levels)]
  
  data_list <- list(
    N   = length(pid_levels),
    S   = length(seq_levels),
    T   = nrow(d),
    pid = d$pid_i,
    sid = d$sid_s,
    y   = d$y
  )
  
  # ----------------------------
  # Sampling settings
  # ----------------------------
  if (dataset == "pilot") {
    iter_val <- 4000; warmup_val <- 2000; chains_val <- 4
    adapt_delta_val <- 0.99; treedepth_val <- 15
  } else if (dataset == "main") {
    iter_val <- 4000; warmup_val <- 2000; chains_val <- 4
    adapt_delta_val <- 0.99; treedepth_val <- 15
  } else {
    stop("dataset must be 'pilot' or 'main'")
  }
  
  msg(
    "RQ1 (", treat_fn, ") Stan settings:",
    " iter=", iter_val,
    " warmup=", warmup_val,
    " chains=", chains_val,
    " adapt_delta=", adapt_delta_val,
    " treedepth=", treedepth_val
  )
  
  sm <- rstan::stan_model(stan_file)
  
  fit <- rstan::sampling(
    sm,
    data = data_list,
    iter = iter_val,
    warmup = warmup_val,
    chains = chains_val,
    seed = seed,
    control = list(
      adapt_delta = adapt_delta_val,
      max_treedepth = treedepth_val
    )
  )
  
  # ----------------------------
  # Save fit + levels (REQUIRED CHANGE)
  # ----------------------------
  f_fit <- file.path(path_mod_ds(dataset), paste0("rq1_fit_sequences_", treat_fn, ".rds"))
  saveRDS(fit, f_fit)
  msg("Saved RQ1 fit:", f_fit)
  
  saveRDS(
    pid_levels,
    file.path(path_mod_ds(dataset), paste0("rq1_pid_levels_", dataset, "_", treat_fn, ".rds"))
  )
  saveRDS(
    seq_levels,
    file.path(path_mod_ds(dataset), paste0("rq1_seq_levels_", dataset, "_", treat_fn, ".rds"))
  )
  
  # ----------------------------
  # Extract sequence-level mu_b draws
  # ----------------------------
  post <- rstan::extract(fit)
  mu_draws <- post$mu_b  # iterations x S
  
  # Trials per sequence
  n_trials_by_seq <- d[, .(n_trials = .N), by = seq]
  setkey(n_trials_by_seq, seq)
  
  # Summary table
  seq_tbl <- data.table(
    sequence    = seq_levels,
    mu_b_median = apply(mu_draws, 2, median),
    mu_b_mean   = apply(mu_draws, 2, mean),
    mu_b_q025   = apply(mu_draws, 2, quantile, probs = 0.025),
    mu_b_q975   = apply(mu_draws, 2, quantile, probs = 0.975)
  )
  
  seq_tbl[, n_trials := n_trials_by_seq[.(sequence), n_trials]]
  seq_tbl[, n_trials := as.integer(n_trials)]
  
  # U_b(s) for each rho
  for (rho in rho_vec) {
    nm <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
    seq_tbl[, (nm) := apply(mu_draws, 2, function(x) mean(x < (1 - rho)))]
  }
  
  # Labeling for rho = 0.10
  rho_main <- 0.10
  col_main <- paste0("U_b_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
  seq_tbl[, underbet_label :=
            fifelse(get(col_main) >= 0.95, "strong",
                    fifelse(get(col_main) >= 0.80, "moderate",
                            fifelse(get(col_main) >= 0.50, "weak", "neutral")))]
  
  setorder(seq_tbl, sequence)
  
  # Save table (treatment-specific filename is safer)
  out_csv <- file.path(
    path_out_ds(dataset),
    paste0("rq1_sequences_", dataset, "_", treat_fn, ".csv")
  )
  fwrite(seq_tbl, out_csv)
  msg("Saved RQ1 sequence table:", out_csv)
  
  invisible(list(
    fit = fit,
    sequences = seq_tbl,
    pid_levels = pid_levels,
    seq_levels = seq_levels
  ))
}

# Example:
# rq1_sequences("pilot", treat_fn = "m25")
# rq1_sequences("main",  treat_fn = "m25")