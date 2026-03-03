library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq3_stan <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Run-level model switch
  # ----------------------------
  rq3_model <- as.character(model$rq3$model_type)
  if (is.na(rq3_model) || !nzchar(rq3_model)) rq3_model <- "hurdle_gamma"
  stopifnot(rq3_model %in% c("hurdle_gamma", "gaussian"))
  
  # Draw subsample size for E_draws approximation
  rq3_krep <- as.integer(model$ppc$rq3_k[[ds]])
  stopifnot(length(rq3_krep) == 1L, rq3_krep >= 10L)
  
  # ----------------------------
  # Design constants
  # ----------------------------
  e     <- as.numeric(design$seq$endowment)
  xmin  <- as.numeric(design$seq$xmin)
  p_win <- as.numeric(design$seq$coin_prob)
  m_map <- design$seq$treatments
  
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  stopifnot(length(xmin) == 1L, is.finite(xmin), xmin > 0)
  stopifnot(length(p_win) == 1L, is.finite(p_win), p_win >= 0, p_win <= 1)
  stopifnot(is.list(m_map), length(m_map) >= 1L)
  
  # ----------------------------
  # Stan settings
  # ----------------------------
  st <- model$stan$rq3[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  msg("RQ3 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed,
      " | rq3_model=", rq3_model,
      " | rq3_krep=", rq3_krep)
  
  # ----------------------------
  # Paths
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- switch(
    rq3_model,
    hurdle_gamma = here::here("stan", "rq3_hurdle_gamma.stan"),
    gaussian     = here::here("stan", "rq3_gaussian.stan")
  )
  stopifnot(file.exists(stan_file))
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Load master once
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  req <- c("pid", "treat", "seq", "stake")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[is.na(stake), stake := 0]
  
  # ----------------------------
  # Load shared r draws (pooled)
  # ----------------------------
  f_r_all <- file.path(mod_dir, "mpl_r_draws.rds")
  stopifnot(file.exists(f_r_all))
  
  r_all <- readRDS(f_r_all)
  stopifnot(is.list(r_all), all(c("pid", "r_draws") %in% names(r_all)))
  stopifnot(is.matrix(r_all$r_draws))
  
  pid_r_levels <- as.character(r_all$pid)
  r_draws_all  <- r_all$r_draws
  K_all <- nrow(r_draws_all)
  stopifnot(K_all >= 10L)
  
  # subsample draws (fixed per run)
  Krep <- min(rq3_krep, K_all)
  set.seed(seed)
  k_idx <- sort(sample.int(K_all, Krep, replace = FALSE))
  msg("RQ3 propagation: using Krep=", Krep, " draws out of K_all=", K_all)
  
  # ----------------------------
  # Compile Stan once
  # ----------------------------
  sm <- rstan::stan_model(stan_file)
  
  # ----------------------------
  # Helper: compute y + fit + write artifacts for (tr,tag)
  # ----------------------------
  fit_one <- function(d0, tr, tag) {
    
    f_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, "_", tag, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq3_pid_levels_",  tr, "_", tag, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq3_seq_levels_",  tr, "_", tag, ".rds"))
    
    if (should_skip(
      paths = c(f_fit, f_pid, f_seq),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ3 Stan (", ds, "/", tr, "/", tag, ")")
    )) return(invisible(NULL))
    
    if (nrow(d0) == 0) return(invisible(NULL))
    
    stopifnot(tr %in% names(m_map))
    m <- as.numeric(m_map[[tr]])
    stopifnot(is.finite(m), m > 1)
    
    # ---- load treatment-specific a* draws (same pattern as RQ2) ----
    f_a <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_a))
    
    a_obj <- readRDS(f_a)
    stopifnot(is.list(a_obj), all(c("pid", "a_star_draws") %in% names(a_obj)))
    stopifnot(is.matrix(a_obj$a_star_draws))
    stopifnot(nrow(a_obj$a_star_draws) == K_all)
    
    pid_a_levels <- as.character(a_obj$pid)
    a_star_all   <- a_obj$a_star_draws
    
    # levels for fit
    pid_levels_fit <- sort(unique(as.character(d0$pid)))
    seq_levels_fit <- sort(unique(as.character(d0$seq)))
    
    d <- copy(d0)
    d[, pid := as.character(pid)]
    d[, seq := as.character(seq)]
    d[, pid_i := match(pid, pid_levels_fit)]
    d[, sid_s := match(seq, seq_levels_fit)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
    
    # map pid -> draw columns
    idx_r <- match(pid_levels_fit, pid_r_levels)
    idx_a <- match(pid_levels_fit, pid_a_levels)
    if (anyNA(idx_r)) stop("RQ3: pid missing from r cache (mpl_r_draws.rds) in tr='", tr, "' tag='", tag, "'.")
    if (anyNA(idx_a)) stop("RQ3: pid missing from a* cache (a_star_draws_", tr, ".rds) in tag='", tag, "'.")
    
    r_draws <- r_draws_all[k_idx, idx_r, drop = FALSE]  # Krep x N_fit
    a_draws <- a_star_all[k_idx, idx_a, drop = FALSE]   # Krep x N_fit
    
    # compute y_t = E_k[max(CE(a*)-CE(a),0)] / e
    Tn <- nrow(d)
    y  <- numeric(Tn)
    
    ii    <- d$pid_i
    a_obs <- as.numeric(d$stake)
    
    for (t in seq_len(Tn)) {
      i <- ii[t]
      
      r_vec <- r_draws[, i]
      a_opt <- a_draws[, i]
      a_act <- rep(a_obs[t], Krep)
      
      ce_opt <- ce_stake_vec_r(a = a_opt, r = r_vec, m = m, e = e, p_win = p_win, xmin = xmin)
      ce_act <- ce_stake_vec_r(a = a_act, r = r_vec, m = m, e = e, p_win = p_win, xmin = xmin)
      
      y[t] <- mean(pmax(ce_opt - ce_act, 0)) / e
    }
    
    if (!all(is.finite(y))) stop("RQ3: non-finite y encountered for tr='", tr, "' tag='", tag, "'.")
    
    # Stan data
    if (rq3_model == "hurdle_gamma") {
      data_list <- list(
        N       = length(pid_levels_fit),
        S       = length(seq_levels_fit),
        T       = Tn,
        pid     = as.integer(d$pid_i),
        sid     = as.integer(d$sid_s),
        y       = as.vector(pmax(y, 0)),
        is_zero = as.integer(y <= 0)
      )
    } else {
      data_list <- list(
        N   = length(pid_levels_fit),
        S   = length(seq_levels_fit),
        T   = Tn,
        pid = as.integer(d$pid_i),
        sid = as.integer(d$sid_s),
        y   = as.vector(pmax(y, 0))
      )
    }
    
    msg("RQ3 Stan: fitting tr=", tr, " tag=", tag,
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
    saveRDS(pid_levels_fit, f_pid)
    saveRDS(seq_levels_fit, f_seq)
    
    msg("Saved: ", f_fit)
    msg("Saved: ", f_pid)
    msg("Saved: ", f_seq)
    
    invisible(list(fit = f_fit, pid = f_pid, seq = f_seq))
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
  # PASS 2: confirmatory subset fits (tau-based a_star_pid_flags)
  #   only where design$a_flags$betting_normative[[tr]] == TRUE
  #   copy full -> conf if subset == full (after FULL pid_levels)
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
    
    # --- full artifacts must exist (PASS 1 ran) ---
    f_fit_full <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, "_full.rds"))
    f_pid_full <- file.path(mod_dir, paste0("rq3_pid_levels_",  tr, "_full.rds"))
    f_seq_full <- file.path(mod_dir, paste0("rq3_seq_levels_",  tr, "_full.rds"))
    stopifnot(file.exists(f_fit_full), file.exists(f_pid_full), file.exists(f_seq_full))
    
    # IMPORTANT: compare to the *actual* full pid set used in FULL fit
    full_pid_prepped <- sort(as.character(readRDS(f_pid_full)))
    keep_pid_prepped <- sort(intersect(full_pid_prepped, keep_pid_raw))
    
    # conf file targets
    f_fit_conf <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, "_conf.rds"))
    f_pid_conf <- file.path(mod_dir, paste0("rq3_pid_levels_",  tr, "_conf.rds"))
    f_seq_conf <- file.path(mod_dir, paste0("rq3_seq_levels_",  tr, "_conf.rds"))
    
    if (identical(full_pid_prepped, keep_pid_prepped)) {
      
      if (should_skip(
        paths = c(f_fit_conf, f_pid_conf, f_seq_conf),
        cfg   = cfg,
        type  = "model",
        label = paste0("RQ3 Stan (", ds, "/", tr, "/conf copy)")
      )) next
      
      ok1 <- file.copy(f_fit_full, f_fit_conf, overwrite = TRUE)
      ok2 <- file.copy(f_pid_full, f_pid_conf, overwrite = TRUE)
      ok3 <- file.copy(f_seq_full, f_seq_conf, overwrite = TRUE)
      stopifnot(ok1, ok2, ok3)
      
      msg("RQ3 Stan: conf subset == full for tr=", tr, " -> copied full artifacts to conf.")
      next
    }
    
    d_conf <- dt[treat == tr & pid %in% keep_pid_raw]
    if (nrow(d_conf) == 0) next
    
    fit_one(d_conf, tr, tag = "conf")
  }
  
  invisible(TRUE)
}