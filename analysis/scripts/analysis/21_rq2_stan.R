# scripts/analysis/21_rq2_stan.R
# ------------------------------------------------------------
# RQ2 (stakes): Preprocess + Stan in ONE step
#
# What this does (per dataset x treatment):
# 1) Load master_sequences.csv
# 2) Keep betting trials only: stake > 0
# 3) Merge cached a*_i draws (indices) and take median a*_i (pipeline test)
# 4) Compute Δa_is = a_is - a*_i
# 5) Within-participant normalization on betting trials:
#      z_is = (Δa_is - mean_i(Δa)) / s*_i,  s*_i = max(sd_i(Δa), sd_floor)
# 6) Exclude participants with < rq2_min_bets betting trials
# 7) Fit Stan model stan/rq2_stakes.stan
# 8) Save:
#    - fit:        data/clean/<ds>/models/rq2_fit_sequences_<tr>.rds
#    - pid levels: data/clean/<ds>/models/rq2_pid_levels_<tr>.rds
#    - seq levels: data/clean/<ds>/models/rq2_seq_levels_<tr>.rds
#    - prepared:   data/clean/<ds>/output/rq2_prepared_<tr>.csv
#
# NOTE:
# - No dataset name in filenames.
# - This script does NOT create sequence/participant summary tables (22_/23_).
# ------------------------------------------------------------

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq2_stan <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  run    <- cfg$run
  design <- cfg$design
  model  <- cfg$model
  
  ds   <- as.character(run$dataset)
  seed <- as.integer(run$seed)
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Design params
  # ----------------------------
  stopifnot(!is.null(design$seq), !is.null(design$seq$endowment))
  e <- as.numeric(design$seq$endowment)
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  
  stopifnot(!is.null(design$exclusion), !is.null(design$exclusion$rq2_min_bets))
  min_bets <- as.integer(design$exclusion$rq2_min_bets)
  stopifnot(length(min_bets) == 1L, min_bets >= 1L)
  
  stopifnot(!is.null(design$rhos), !is.null(design$rhos$rq2_sd_floor))
  sd_floor <- as.numeric(design$rhos$rq2_sd_floor)
  stopifnot(length(sd_floor) == 1L, is.finite(sd_floor), sd_floor > 0)
  
  # ----------------------------
  # Stan sampling settings
  # ----------------------------
  stopifnot(!is.null(model$stan), !is.null(model$stan$rq2), !is.null(model$stan$rq2[[ds]]))
  st <- model$stan$rq2[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  msg("RQ2 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed)
  
  # ----------------------------
  # Files / dirs
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- here::here("stan", "rq2_stakes.stan")
  stopifnot(file.exists(stan_file))
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Pre-check artifacts BEFORE compiling Stan
  # (fit + levels + prepared csv)
  # ----------------------------
  to_run <- character(0)
  
  for (tr in tr_vec) {
    
    f_fit        <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq2_seq_levels_", tr, ".rds"))
    f_prep_csv   <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    
    if (file.exists(f_fit) || file.exists(f_pid_levels) || file.exists(f_seq_levels) || file.exists(f_prep_csv)) {
      warning(
        "RQ2 results already exist for dataset='", ds, "', treatment='", tr, "'.\n",
        "Stan was NOT executed.\n",
        "To rerun, delete (as applicable):\n",
        "  ", f_fit, "\n",
        "  ", f_pid_levels, "\n",
        "  ", f_seq_levels, "\n",
        "  ", f_prep_csv
      )
    } else {
      to_run <- c(to_run, tr)
    }
  }
  
  if (length(to_run) == 0L) {
    msg("RQ2: No treatments require estimation. Nothing to run.")
    return(invisible(NULL))
  }
  
  # ----------------------------
  # Load input once
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "seq", "stake")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      "master_sequences.csv missing columns: ", paste(missing, collapse = ", "),
      "\nExpected at least: ", paste(required, collapse = ", ")
    )
  }
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  
  # intensive margin only
  dt <- dt[!is.na(stake) & stake > 0]
  if (nrow(dt) == 0) stop("RQ2: No betting trials (stake > 0) in master_sequences.csv for dataset='", ds, "'.")
  
  # ----------------------------
  # Compile Stan once (only if needed)
  # ----------------------------
  sm <- rstan::stan_model(stan_file)
  
  # ----------------------------
  # Run per treatment
  # ----------------------------
  outputs <- list()
  
  for (tr in to_run) {
    
    # output paths (NO dataset in filename)
    f_fit        <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq2_seq_levels_", tr, ".rds"))
    f_prep_csv   <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    
    # ----------------------------
    # Load cached a* draws (indices)
    # ----------------------------
    f_astar <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    if (!file.exists(f_astar)) {
      stop(
        "a* cache not found for dataset='", ds, "', treatment='", tr, "':\n  ", f_astar, "\n",
        "Run cache_a_star_from_r_draws() in indices first."
      )
    }
    
    ast <- readRDS(f_astar)
    stopifnot(is.list(ast), all(c("pid", "a_star_draws") %in% names(ast)))
    stopifnot(is.matrix(ast$a_star_draws))
    
    pid_astar <- as.character(ast$pid)
    a_star_med <- apply(ast$a_star_draws, 2, median)
    a_map <- data.table(pid = pid_astar, a_star = as.numeric(a_star_med))
    
    # ----------------------------
    # Filter to treatment + merge a*
    # ----------------------------
    d <- dt[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ2: No betting trials after filtering treat == '", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    
    d <- merge(d, a_map, by = "pid", all.x = TRUE)
    if (anyNA(d$a_star)) {
      bad <- unique(d[is.na(a_star), pid])
      stop("Missing a* for some pids in treatment='", tr, "'. Example: ", paste(head(bad, 10), collapse = ", "))
    }
    
    # Δa = a - a*
    d[, delta_a := stake - a_star]
    
    # participant stats on betting trials
    pid_stats <- d[, .(
      n_bets    = .N,
      delta_bar = mean(delta_a, na.rm = TRUE),
      delta_sd  = sd(delta_a,   na.rm = TRUE)
    ), by = pid]
    
    pid_stats[, sd_star := pmax(delta_sd, sd_floor)]
    pid_stats[, keep := (n_bets >= min_bets)]
    
    keep_pid <- pid_stats[keep == TRUE, pid]
    d <- d[pid %in% keep_pid]
    
    if (nrow(d) == 0) {
      warning("RQ2: All participants excluded by rq2_min_bets=", min_bets, " for treatment='", tr, "'. Skipping.")
      next
    }
    
    # attach constants + z
    d <- merge(d, pid_stats[, .(pid, n_bets, delta_bar, sd_star)], by = "pid", all.x = TRUE)
    stopifnot(!anyNA(d$delta_bar), !anyNA(d$sd_star))
    stopifnot(all(is.finite(d$delta_bar)), all(is.finite(d$sd_star)), all(d$sd_star > 0))
    
    d[, z := (delta_a - delta_bar) / sd_star]
    d <- d[is.finite(z)]
    if (nrow(d) == 0) stop("RQ2: All z became non-finite after standardization for treatment='", tr, "'.")
    
    # ----------------------------
    # Save prepared CSV (for 22_/23_)
    # ----------------------------
    fwrite(
      d[, .(pid, treat, seq, stake, a_star, delta_a, delta_bar, sd_star, z, n_bets)],
      f_prep_csv
    )
    msg("Saved: ", f_prep_csv)
    
    # ----------------------------
    # Build Stan data (must match rq2_stakes.stan)
    # ----------------------------
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
    
    # participant-level constants aligned to pid_levels
    pid_const <- unique(d[, .(pid, delta_bar, sd_star)])
    setkey(pid_const, pid)
    pid_const <- pid_const[.(pid_levels)]
    stopifnot(!anyNA(pid_const$pid))
    
    data_list <- list(
      N         = length(pid_levels),
      S         = length(seq_levels),
      T         = nrow(d),
      pid       = as.integer(d$pid_i),
      sid       = as.integer(d$sid_s),
      z         = as.vector(d$z),
      delta_bar = as.vector(pid_const$delta_bar),
      s_star    = as.vector(pid_const$sd_star),
      e         = as.numeric(e)
    )
    
    msg("Fitting RQ2 Stan for treatment=", tr,
        " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T)
    
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
    
    # Do NOT save empty fits
    if (length(rstan::get_sampler_params(fit, inc_warmup = FALSE)) == 0) {
      stop("RQ2: Stan produced no samples for treatment='", tr, "'. Fit will NOT be saved.")
    }
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels, f_pid_levels)
    saveRDS(seq_levels, f_seq_levels)
    
    msg("Saved RQ2 fit: ", f_fit)
    msg("Saved RQ2 pid levels: ", f_pid_levels)
    msg("Saved RQ2 seq levels: ", f_seq_levels)
    
    outputs[[tr]] <- list(
      fit_file = f_fit,
      pid_levels_file = f_pid_levels,
      seq_levels_file = f_seq_levels,
      prepared_csv = f_prep_csv,
      N = data_list$N,
      S = data_list$S,
      T = data_list$T
    )
  }
  
  invisible(outputs)
}

# Example:
# rq2_stan(cfg)