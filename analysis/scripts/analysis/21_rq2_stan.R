library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq2_stan <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  design <- cfg$design
  model  <- cfg$model
  
  # ----------------------------
  # Design parameters
  # ----------------------------
  e        <- as.numeric(design$seq$endowment)
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
  infile    <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stan_file <- here::here("stan", "rq2_stakes.stan")
  stopifnot(file.exists(infile), file.exists(stan_file))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Load master once
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid","treat","seq","stake") %in% names(dt)))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  
  # Intensive margin only
  dt <- dt[!is.na(stake) & stake > 0]
  if (nrow(dt) == 0) stop("RQ2: No betting trials (stake > 0).")
  
  # Compile once
  sm <- rstan::stan_model(stan_file)
  
  # ----------------------------
  # Helper: prepare + fit + write artifacts
  # ----------------------------
  fit_one <- function(d0, tr, tag) {
    
    f_fit <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, "_", tag, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq2_pid_levels_",  tr, "_", tag, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq2_seq_levels_",  tr, "_", tag, ".rds"))
    f_prep <- file.path(out_dir, paste0("rq2_prepared_", tr, "_", tag, ".csv"))
    
    if (should_skip(
      paths = c(f_fit, f_pid, f_seq),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ2 Stan (", ds, "/", tr, "/", tag, ")")
    )) return(invisible(NULL))
    
    if (should_skip(
      paths = f_prep,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ2 prepared (", ds, "/", tr, "/", tag, ")")
    )) return(invisible(NULL))
    
    if (nrow(d0) == 0) return(invisible(NULL))
    
    # ---- a* posterior mean per pid (treatment-specific) ----
    f_astar <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_astar))
    
    ast <- readRDS(f_astar)
    stopifnot(is.list(ast), all(c("pid","a_star_draws") %in% names(ast)))
    stopifnot(is.matrix(ast$a_star_draws))
    
    a_star_mean <- colMeans(ast$a_star_draws)
    a_map <- data.table(
      pid    = as.character(ast$pid),
      a_star = as.numeric(a_star_mean)
    )
    
    # ---- merge a* and compute delta ----
    d <- merge(d0, a_map, by = "pid", all.x = TRUE)
    stopifnot(!anyNA(d$a_star))
    
    d[, delta_a := stake - a_star]
    
    # ---- within-pid stats and min bets filter ----
    pid_stats <- d[, .(
      n_bets    = .N,
      delta_bar = mean(delta_a),
      delta_sd  = sd(delta_a)
    ), by = pid]
    
    pid_stats[, sd_star := pmax(delta_sd, sd_floor)]
    pid_stats[, keep := (n_bets >= min_bets)]
    
    # ---- save excluded participants (< min_bets) ----
    excluded_pid <- pid_stats[keep == FALSE, .(pid, n_bets)]
    f_excl <- file.path(out_dir, paste0("rq2_excluded_minbets_", tr, "_", tag, ".csv"))
    
    if (!should_skip(
      paths = f_excl,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ2 excluded (<min_bets) (", ds, "/", tr, "/", tag, ")")
    )) {
      fwrite(excluded_pid, f_excl)
      msg("Saved: ", f_excl)
    }
    
    # ---- apply filter ----
    d <- d[pid %in% pid_stats[keep == TRUE, pid]]
    if (nrow(d) == 0) return(invisible(NULL))
    
    d <- merge(
      d,
      pid_stats[, .(pid, n_bets, delta_bar, sd_star)],
      by = "pid",
      all.x = TRUE
    )
    
    d[, z := (delta_a - delta_bar) / sd_star]
    d <- d[is.finite(z)]
    if (nrow(d) == 0) return(invisible(NULL))
    
    # ---- write prepared (required by tables later) ----
    fwrite(
      d[, .(pid, treat, seq, stake, a_star, delta_a, delta_bar, sd_star, z, n_bets)],
      f_prep
    )
    msg("Saved: ", f_prep)
    
    # ---- build Stan data ----
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
    
    msg("RQ2 Stan: fitting tr=", tr, " tag=", tag,
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
    
    invisible(list(fit = f_fit, pid = f_pid, seq = f_seq, prep = f_prep))
  }
  
  # ----------------------------
  # PASS 1: full sample fits
  # ----------------------------
  for (tr in tr_vec) {
    d_full <- dt[treat == tr]
    if (nrow(d_full) == 0) next
    fit_one(d_full, tr, tag = "full")
  }
  
  # ----------------------------
  # PASS 2: confirmatory subset fits
  #   only where design$a_flags$betting_normative[[tr]] == TRUE
  #   copy full -> conf if subset == full (after RQ2 preprocessing)
  # ----------------------------
  tau_main <- as.numeric(design$a_flags$tau[1])
  stopifnot(is.finite(tau_main), tau_main > 0, tau_main < 1)
  tau_nm <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  for (tr in tr_vec) {
    
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    
    f_flags <- file.path(mod_dir, paste0("a_star_pid_flags_", tr, ".rds"))
    stopifnot(file.exists(f_flags))
    
    flags <- readRDS(f_flags)
    stopifnot(is.list(flags), !is.null(flags$pid_sets))
    
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    
    keep_pid_raw <- as.character(flags$pid_sets[[nm_keep]])
    
    # --- full artifacts (must exist, since PASS 1 ran) ---
    f_fit_full  <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, "_full.rds"))
    f_pid_full  <- file.path(mod_dir, paste0("rq2_pid_levels_",  tr, "_full.rds"))
    f_seq_full  <- file.path(mod_dir, paste0("rq2_seq_levels_",  tr, "_full.rds"))
    f_prep_full <- file.path(out_dir, paste0("rq2_prepared_", tr, "_full.csv"))
    
    stopifnot(file.exists(f_fit_full), file.exists(f_pid_full), file.exists(f_seq_full), file.exists(f_prep_full))
    
    # IMPORTANT: compare to the *preprocessed* full pid set (min_bets + finite z already applied)
    full_pid_prepped <- sort(as.character(readRDS(f_pid_full)))
    
    # confirmatory pid set after preprocessing is just intersection
    keep_pid_prepped <- sort(intersect(full_pid_prepped, keep_pid_raw))
    
    # conf file targets
    f_fit_conf  <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, "_conf.rds"))
    f_pid_conf  <- file.path(mod_dir, paste0("rq2_pid_levels_",  tr, "_conf.rds"))
    f_seq_conf  <- file.path(mod_dir, paste0("rq2_seq_levels_",  tr, "_conf.rds"))
    f_prep_conf <- file.path(out_dir, paste0("rq2_prepared_", tr, "_conf.csv"))
    
    # if identical, copy full -> conf (unless conf artifacts should be skipped)
    if (identical(full_pid_prepped, keep_pid_prepped)) {
      
      if (should_skip(
        paths = c(f_fit_conf, f_pid_conf, f_seq_conf),
        cfg   = cfg,
        type  = "model",
        label = paste0("RQ2 Stan (", ds, "/", tr, "/conf copy)")
      )) next
      
      if (should_skip(
        paths = f_prep_conf,
        cfg   = cfg,
        type  = "output",
        label = paste0("RQ2 prepared (", ds, "/", tr, "/conf copy)")
      )) next
      
      ok1 <- file.copy(f_fit_full,  f_fit_conf,  overwrite = TRUE)
      ok2 <- file.copy(f_pid_full,  f_pid_conf,  overwrite = TRUE)
      ok3 <- file.copy(f_seq_full,  f_seq_conf,  overwrite = TRUE)
      ok4 <- file.copy(f_prep_full, f_prep_conf, overwrite = TRUE)
      stopifnot(ok1, ok2, ok3, ok4)
      
      msg("RQ2 Stan: conf subset == full for tr=", tr, " -> copied full artifacts to conf.")
      next
    }
    
    # otherwise: fit on subset (raw keep set; preprocessing inside fit_one will drop min_bets etc.)
    d_conf <- dt[treat == tr & pid %in% keep_pid_raw]
    if (nrow(d_conf) == 0) next
    
    fit_one(d_conf, tr, tag = "conf")
  }
  
  invisible(TRUE)
}