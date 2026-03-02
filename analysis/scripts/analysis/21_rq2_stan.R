# ------------------------------------------------------------
# RQ2 (stakes): Preprocess + Stan in ONE step
#
# What this does (per dataset x treatment):
# 1) Load master_sequences.csv
# 2) Keep betting trials only: stake > 0
# 3) Merge cached a*_i draws and use posterior mean a*_i
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
# - Full sample is used for fitting.
# - Confirmatory subset is handled in post-processing (23_).
# ------------------------------------------------------------

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq2_stan <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  # ----------------------------
  # Design parameters
  # ----------------------------
  e <- as.numeric(design$seq$endowment)
  min_bets <- as.integer(design$rq2$min_bets)
  sd_floor <- as.numeric(design$rq2$sd_floor)
  
  stopifnot(length(e) == 1L, e > 0)
  stopifnot(length(min_bets) == 1L, min_bets >= 1L)
  stopifnot(length(sd_floor) == 1L, sd_floor > 0)
  
  # ----------------------------
  # Stan sampling settings
  # ----------------------------
  st <- model$stan$rq2[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  msg("RQ2 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed)
  
  # ----------------------------
  # Paths
  # ----------------------------
  infile   <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stan_file <- here::here("stan", "rq2_stakes.stan")
  
  stopifnot(file.exists(infile))
  stopifnot(file.exists(stan_file))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Load input
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "seq", "stake")
  stopifnot(all(required %in% names(dt)))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  
  # Intensive margin only
  dt <- dt[!is.na(stake) & stake > 0]
  if (nrow(dt) == 0) stop("RQ2: No betting trials (stake > 0).")
  
  # Compile Stan once
  sm <- rstan::stan_model(stan_file)
  
  outputs <- list()
  
  # ----------------------------
  # Loop over treatments
  # ----------------------------
  for (tr in tr_vec) {
    
    f_fit        <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq2_seq_levels_", tr, ".rds"))
    f_prep_csv   <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    
    if (should_skip(
      paths = c(f_fit, f_pid_levels, f_seq_levels),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ2 Stan (", ds, "/", tr, ")")
    )) next
    
    if (should_skip(
      paths = c(f_prep_csv),
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ2 Stan (", ds, "/", tr, ")")
    )) next
    
    # ----------------------------
    # Load cached a*_i draws
    # ----------------------------
    f_astar <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_astar))
    
    ast <- readRDS(f_astar)
    stopifnot(is.matrix(ast$a_star_draws))
    
    a_star_mean <- colMeans(ast$a_star_draws)
    a_map <- data.table(
      pid = as.character(ast$pid),
      a_star = as.numeric(a_star_mean)
    )
    
    # ----------------------------
    # Filter + merge a*
    # ----------------------------
    d <- dt[treat == tr]
    if (nrow(d) == 0) next
    
    d <- merge(d, a_map, by = "pid", all.x = TRUE)
    stopifnot(!anyNA(d$a_star))
    
    d[, delta_a := stake - a_star]
    
    # Participant stats
    # Participant stats
    pid_stats <- d[, .(
      n_bets    = .N,
      delta_bar = mean(delta_a),
      delta_sd  = sd(delta_a)
    ), by = pid]
    
    pid_stats[, sd_star := pmax(delta_sd, sd_floor)]
    pid_stats[, keep := (n_bets >= min_bets)]
    
    d <- d[pid %in% pid_stats[keep == TRUE, pid]]
    if (nrow(d) == 0) next
    
    # attach participant constants INCLUDING n_bets
    d <- merge(
      d,
      pid_stats[, .(pid, n_bets, delta_bar, sd_star)],
      by = "pid",
      all.x = TRUE
    )
    
    d[, z := (delta_a - delta_bar) / sd_star]
    d <- d[is.finite(z)]
    
    # WRITE prepared (clean schema, required by 23_)
    fwrite(
      d[, .(pid, treat, seq, stake, a_star, delta_a,
            delta_bar, sd_star, z, n_bets)],
      f_prep_csv
    )
    msg("Saved: ", f_prep_csv)
    
    # ----------------------------
    # Build Stan data
    # ----------------------------
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    pid_const <- unique(d[, .(pid, delta_bar, sd_star)])
    setkey(pid_const, pid)
    pid_const <- pid_const[.(pid_levels)]
    
    data_list <- list(
      N         = length(pid_levels),
      S         = length(seq_levels),
      T         = nrow(d),
      pid       = as.integer(d$pid_i),
      sid       = as.integer(d$sid_s),
      z         = as.vector(d$z),
      delta_bar = as.vector(pid_const$delta_bar),
      s_star    = as.vector(pid_const$sd_star),
      e         = as.integer(e)
    )
    
    msg("Fitting RQ2 Stan for treatment=", tr,
        " | N=", data_list$N,
        " | S=", data_list$S,
        " | T=", data_list$T)
    
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
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels, f_pid_levels)
    saveRDS(seq_levels, f_seq_levels)
    
    msg("Saved RQ2 fit: ", f_fit)
    
    outputs[[tr]] <- list(
      fit_file = f_fit,
      prepared_csv = f_prep_csv,
      N = data_list$N,
      S = data_list$S,
      T = data_list$T
    )
  }
  
  invisible(outputs)
}