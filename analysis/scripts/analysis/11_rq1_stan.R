rq1_stan <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  design <- cfg$design
  model  <- cfg$model
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stan_file <- here::here("stan", "rq1_bets.stan")
  stopifnot(file.exists(infile), file.exists(stan_file))
  
  # Stan settings
  st <- model$stan$rq1[[ds]]
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load master
  dt <- data.table::fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid","treat","stake","seq") %in% names(dt)))
  
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, y := as.integer(!is.na(stake) & as.numeric(stake) > 0)]
  
  # Compile once
  sm <- rstan::stan_model(stan_file)
  
  # Helper: fit one dataset slice and write artifacts
  fit_one <- function(d, tr, tag) {
    
    f_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, "_", tag, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, "_", tag, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, "_", tag, ".rds"))
    
    if (should_skip(
      paths = c(f_fit, f_pid, f_seq),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ1 Stan (", ds, "/", tr, "/", tag, ")")
    )) return(invisible(NULL))
    
    if (nrow(d) == 0) return(invisible(NULL))
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    stopifnot(length(seq_levels) == 64L)
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    data_list <- list(
      N   = length(pid_levels),
      S   = length(seq_levels),
      T   = nrow(d),
      pid = as.integer(d$pid_i),
      sid = as.integer(d$sid_s),
      y   = as.integer(d$y)
    )
    
    msg("RQ1 Stan: fitting tr=", tr, " tag=", tag,
        " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T)
    
    fit <- rstan::sampling(
      sm,
      data = data_list,
      iter = iter_val,
      warmup = warmup_val,
      chains = chains_val,
      seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels, f_pid)
    saveRDS(seq_levels, f_seq)
    
    msg("Saved: ", f_fit)
    invisible(list(fit = f_fit, pid = f_pid, seq = f_seq))
  }
  
  # ----------------------------
  # PASS 1: full sample (per treatment)
  # ----------------------------
  for (tr in tr_vec) {
    d_full <- dt[treat == tr]
    if (nrow(d_full) == 0) next
    fit_one(d_full, tr, tag = "full")
  }
  
  # ----------------------------
  # PASS 2: confirmatory subset (only where betting_normative == TRUE)
  # If confirmatory pid set == full pid set, copy full artifacts -> conf.
  # ----------------------------
  tau_main <- as.numeric(design$a_flags$tau[1])   # first = main
  tau_nm <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  for (tr in tr_vec) {
    
    # run CONF only where benchmark is "betting is optimal"
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    
    f_rds_flags <- file.path(mod_dir, paste0("a_star_pid_flags_", tr, ".rds"))
    stopifnot(file.exists(f_rds_flags))
    
    flags <- readRDS(f_rds_flags)
    stopifnot(is.list(flags), !is.null(flags$pid_sets))
    
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    
    keep_pid <- as.character(flags$pid_sets[[nm_keep]])
    
    # full pid set in this treatment (based on master)
    full_pid <- sort(unique(dt[treat == tr, pid]))
    keep_pid <- sort(unique(keep_pid))
    
    # file targets
    f_fit_full <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, "_full.rds"))
    f_pid_full <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, "_full.rds"))
    f_seq_full <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, "_full.rds"))
    
    f_fit_conf <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, "_conf.rds"))
    f_pid_conf <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, "_conf.rds"))
    f_seq_conf <- file.path(mod_dir, paste0("rq1_seq_levels_",  tr, "_conf.rds"))
    
    # if identical, copy full -> conf (unless conf artifacts should be skipped)
    if (identical(full_pid, keep_pid)) {
      
      if (should_skip(
        paths = c(f_fit_conf, f_pid_conf, f_seq_conf),
        cfg   = cfg,
        type  = "model",
        label = paste0("RQ1 Stan (", ds, "/", tr, "/conf copy)")
      )) next
      
      stopifnot(file.exists(f_fit_full), file.exists(f_pid_full), file.exists(f_seq_full))
      
      ok1 <- file.copy(f_fit_full, f_fit_conf, overwrite = TRUE)
      ok2 <- file.copy(f_pid_full, f_pid_conf, overwrite = TRUE)
      ok3 <- file.copy(f_seq_full, f_seq_conf, overwrite = TRUE)
      stopifnot(ok1, ok2, ok3)
      
      msg("RQ1 Stan: conf subset == full for tr=", tr, " -> copied full artifacts to conf.")
      next
    }
    
    # otherwise: fit on subset
    d_conf <- dt[treat == tr & pid %in% keep_pid]
    if (nrow(d_conf) == 0) next
    fit_one(d_conf, tr, tag = "conf")
  }
  
  invisible(TRUE)
}