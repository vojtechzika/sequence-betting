# ============================================================
# 21_rq2_stan.R
#
# PURPOSE
#   Fits hierarchical stake deviation models for RQ2 (intensive
#   margin of betting). Always fits both primary and alternative
#   models. Sensitivity analyses for confirmatory subset only.
#   Model selection handled downstream in rq2_tables().
#
# CONFIRMATORY SUBSET FOR RQ2:
#   Participants for whom betting is EU-optimal (P(a* > 0) >= tau).
#   Controlled by cfg$design$a_flags$betting_normative and
#   cfg$design$a_flags$tau.
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/a_star_draws_<tr>.rds       -- optimal stakes
#   path_mod/a_star_pid_flags_<tr>.rds   -- confirmatory subset
#   path_mod/drift_decisions.rds         -- drift adjustment
#   stan/rq2_primary.stan                -- Gaussian
#   stan/rq2_alternative.stan            -- one-inflated Beta-Binomial
#
# OUTPUT
#   path_mod/rq2_fit_sequences_<tr>_<tag>.rds           -- primary
#   path_mod/rq2_fit_sequences_<tr>_<tag>_alt.rds       -- alternative
#   path_mod/rq2_fit_sequences_<tr>_<tag>_floor*.rds    -- sensitivity
#   path_mod/rq2_fit_sequences_<tr>_<tag>_mad.rds       -- sensitivity
#   path_mod/rq2_prepared_<tr>_<tag>[_alt|_floor*|_mad].rds
#   path_out/rq2_excluded_minbets_<tr>_<tag>[...].csv
#
# TAGS
#   full          -- all participants in treatment
#   confirmatory  -- normative betters only
#
# SENSITIVITY (confirmatory only):
#   _floor<n> -- alternative sd_floor
#   _mad      -- MAD dispersion measure
#
# CALL ORDER IN PIPELINE:
#   rq2_stan(cfg)         -- primary + alternative + sensitivity
#   rq2_diagnostics(cfg)  -- determines selected model
#   rq2_tables(cfg)       -- produces tables from selected model
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq2_stan <- function(cfg) {
  
  # ============================================================
  # MODEL SPECIFICATION
  #   primary     : Gaussian on standardized stake deviations
  #   alternative : one-inflated Beta-Binomial on raw stakes
  #   sensitivity : alternative sd floors and MAD (confirmatory only)
  # ============================================================
  models <- list(
    list(label     = "primary",
         suffix    = "",
         stan_file = here::here("stan", "rq2_primary.stan"),
         desc      = "Gaussian"),
    list(label     = "alternative",
         suffix    = "_alt",
         stan_file = here::here("stan", "rq2_alternative.stan"),
         desc      = "Beta-Binomial")
  )
  # ============================================================
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  model  <- cfg$model
  
  e                  <- as.numeric(design$seq$endowment)
  min_bets           <- as.integer(design$rq2$min_bets)
  sd_floor           <- as.numeric(design$rq2$sd_floor)
  sd_floor_sens_low  <- as.numeric(design$rq2$sd_floor_sens_low)
  sd_floor_sens_high <- as.numeric(design$rq2$sd_floor_sens_high)
  
  stopifnot(length(e) == 1L, e > 0)
  stopifnot(length(min_bets) == 1L, min_bets >= 1L)
  stopifnot(length(sd_floor) == 1L, sd_floor > 0)
  
  for (m in models) stopifnot(file.exists(m$stan_file))
  
  infile  <- file.path(path_src, "master_sequences.csv")
  f_drift <- file.path(path_mod, "drift_decisions.rds")
  stopifnot(file.exists(infile), file.exists(f_drift))
  
  drift_decisions <- readRDS(f_drift)
  
  st              <- model$stan$rq2
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "stake", "block") %in% names(dt)))
  
  dt[, pid     := as.character(pid)]
  dt[, treat   := as.character(treat)]
  dt[, seq     := as.character(seq)]
  dt[, stake   := as.numeric(stake)]
  dt[, block_c := as.numeric(block) - 2.5]
  
  dt <- dt[!is.na(stake) & stake > 0]
  if (nrow(dt) == 0) stop("RQ2: No betting trials (stake > 0).")
  
  # ---- Get drift config ----
  get_drift_cfg <- function(tr, tag) {
    key <- paste(tr, tag, "stake", sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) return(list(drift = FALSE, type = "none"))
    list(drift = isTRUE(res$drift), type = res$drift_type)
  }
  
  # ---- Prepare z scores ----
  prepare_z <- function(d0, tr, floor_val, dispersion = "sd") {
    
    f_astar <- file.path(path_mod, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_astar))
    
    ast <- readRDS(f_astar)
    stopifnot(is.list(ast), all(c("pid", "a_star_draws") %in% names(ast)))
    
    a_map <- data.table(
      pid    = as.character(ast$pid),
      a_star = as.numeric(colMeans(ast$a_star_draws))
    )
    
    d <- merge(d0, a_map, by = "pid", all.x = TRUE)
    
    n_missing <- sum(is.na(d$a_star))
    if (n_missing > 0) {
      msg("RQ2 prepare_z: dropping ", n_missing,
          " trials with no a* match (likely inconsistent HL participants)")
      d <- d[!is.na(a_star)]
    }
    if (nrow(d) == 0) return(NULL)
    
    d[, delta_a := stake - a_star]
    
    pid_stats <- d[, {
      dbar <- mean(delta_a)
      disp <- if (dispersion == "mad") {
        median(abs(delta_a - median(delta_a)))
      } else {
        sd(delta_a)
      }
      .(n_bets = .N, delta_bar = dbar, delta_sd = disp)
    }, by = pid]
    
    pid_stats[, sd_star := pmax(delta_sd, floor_val)]
    pid_stats[, keep    := (n_bets >= min_bets)]
    
    d <- d[pid %in% pid_stats[keep == TRUE, pid]]
    if (nrow(d) == 0) return(NULL)
    
    d <- merge(d, pid_stats[, .(pid, n_bets, delta_bar, sd_star)],
               by = "pid", all.x = TRUE)
    
    d[, z := (delta_a - delta_bar) / sd_star]
    d <- d[is.finite(z)]
    
    stopifnot(all(c("pid", "treat", "seq", "stake", "a_star",
                    "delta_a", "delta_bar", "sd_star", "z", "n_bets") %in% names(d)))
    d
  }
  
  # ---- Fit one subset for one model ----
  fit_one <- function(d0, tr, tag, m, floor_val = sd_floor, dispersion = "sd") {
    
    suffix <- m$suffix
    
    f_fit  <- file.path(path_mod, paste0("rq2_fit_sequences_", tr, "_", tag, suffix, ".rds"))
    f_pid  <- file.path(path_mod, paste0("rq2_pid_levels_",    tr, "_", tag, suffix, ".rds"))
    f_seq  <- file.path(path_mod, paste0("rq2_seq_levels_",    tr, "_", tag, suffix, ".rds"))
    f_prep <- file.path(path_mod, paste0("rq2_prepared_",      tr, "_", tag, suffix, ".rds"))
    f_excl <- file.path(path_out, paste0("rq2_excluded_minbets_", tr, "_", tag, suffix, ".csv"))
    
    skip_model <- should_skip(c(f_fit, f_pid, f_seq), cfg, "model",
                              paste0("RQ2 Stan (", tr, "/", tag, "/", m$label, ")"))
    skip_prep  <- should_skip(f_prep, cfg, "output",
                              paste0("RQ2 prepared (", tr, "/", tag, "/", m$label, ")"))
    
    if (skip_model && skip_prep) return(invisible(NULL))
    if (nrow(d0) == 0) return(invisible(NULL))
    
    d <- prepare_z(d0, tr, floor_val, dispersion)
    if (is.null(d) || nrow(d) == 0) return(invisible(NULL))
    
    # Save excluded participants
    excl_pids <- setdiff(unique(d0$pid), unique(d$pid))
    if (!should_skip(f_excl, cfg, "output",
                     paste0("RQ2 excluded (", tr, "/", tag, "/", m$label, ")"))) {
      fwrite(data.table(pid = excl_pids), f_excl)
      msg("Saved: ", f_excl)
    }
    
    if (!skip_prep) {
      saveRDS(d[, .(pid, treat, seq, stake, a_star, delta_a,
                    delta_bar, sd_star, z, n_bets)], f_prep)
      msg("Saved: ", f_prep)
    }
    
    if (!skip_model) {
      pid_levels <- sort(unique(d$pid))
      seq_levels <- sort(unique(d$seq))
      
      d[, pid_i := match(pid, pid_levels)]
      d[, sid_s := match(seq, seq_levels)]
      
      pid_const <- unique(d[, .(pid, delta_bar, sd_star)])
      setkey(pid_const, pid)
      pid_const <- pid_const[.(pid_levels)]
      
      drift_cfg      <- get_drift_cfg(tr, tag)
      include_drift  <- as.integer(isTRUE(drift_cfg$drift))
      prior_gamma_sd <- if (include_drift) {
        as.numeric(cfg$design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
      } else {
        0.3
      }
      
      if (m$label == "primary") {
        data_list <- list(
          N              = length(pid_levels),
          S              = length(seq_levels),
          T              = nrow(d),
          pid            = as.integer(d$pid_i),
          sid            = as.integer(d$sid_s),
          z              = as.vector(d$z),
          delta_bar      = as.vector(pid_const$delta_bar),
          s_star         = as.vector(pid_const$sd_star),
          e              = as.integer(e),
          include_drift  = include_drift,
          block_c        = as.numeric(d$block_c),
          prior_gamma_sd = prior_gamma_sd
        )
      } else {
        # BB alternative uses raw stakes and a_star
        a_star_pid <- d[, .(a_star_mean = mean(a_star)), by = pid]
        setkey(a_star_pid, pid)
        pid_const[, a_star_mean := a_star_pid[.(pid_levels), a_star_mean]]
        
        data_list <- list(
          N              = length(pid_levels),
          S              = length(seq_levels),
          T              = nrow(d),
          pid            = as.integer(d$pid_i),
          sid            = as.integer(d$sid_s),
          y              = as.integer(d$stake),
          a_star         = as.vector(pid_const$a_star_mean),
          e              = as.integer(e),
          include_drift  = include_drift,
          block_c        = as.numeric(d$block_c),
          prior_gamma_sd = prior_gamma_sd
        )
      }
      
      sm <- rstan::stan_model(m$stan_file)
      
      msg("RQ2 Stan: fitting tr=", tr, " tag=", tag,
          " [", m$desc, "]",
          " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T,
          " | drift=", drift_cfg$type,
          " | floor=", floor_val, " | dispersion=", dispersion)
      
      fit <- rstan::sampling(
        sm,
        data    = data_list,
        iter    = iter_val,
        warmup  = warmup_val,
        chains  = chains_val,
        seed    = seed,
        control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
      )
      
      saveRDS(fit, f_fit)
      saveRDS(pid_levels, f_pid)
      saveRDS(seq_levels, f_seq)
      msg("Saved: ", f_fit)
    }
    
    invisible(list(fit = f_fit, pid = f_pid, seq = f_seq, prep = f_prep))
  }
  
  # ---- Copy full to confirmatory (all models) ----
  copy_full_to_confirmatory <- function(tr, suffix = "") {
    for (nm in c("rq2_fit_sequences", "rq2_pid_levels", "rq2_seq_levels", "rq2_prepared")) {
      f_full <- file.path(path_mod, paste0(nm, "_", tr, "_full",         suffix, ".rds"))
      f_conf <- file.path(path_mod, paste0(nm, "_", tr, "_confirmatory", suffix, ".rds"))
      if (file.exists(f_full)) file.copy(f_full, f_conf, overwrite = TRUE)
    }
    msg("RQ2 Stan: confirmatory == full for tr=", tr, suffix, " -> copied.")
  }
  
  tau_main <- as.numeric(design$a_flags$tau[1])
  tau_nm   <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  # ---- PASS 1: full sample (primary + alternative) ----
  for (tr in tr_vec) {
    d_full <- dt[treat == tr]
    if (nrow(d_full) == 0) next
    for (m in models) fit_one(d_full, tr, "full", m)
  }
  
  # ---- PASS 2: confirmatory subset ----
  for (tr in tr_vec) {
    
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    
    f_flags <- file.path(path_mod, paste0("a_star_pid_flags_", tr, ".rds"))
    stopifnot(file.exists(f_flags))
    
    flags   <- readRDS(f_flags)
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    
    keep_pid_raw <- as.character(flags$pid_sets[[nm_keep]])
    f_pid_full   <- file.path(path_mod, paste0("rq2_pid_levels_", tr, "_full.rds"))
    
    d_conf <- dt[treat == tr & pid %in% keep_pid_raw]
    if (nrow(d_conf) == 0) next
    
    if (file.exists(f_pid_full)) {
      full_pid_prepped <- sort(as.character(readRDS(f_pid_full)))
      keep_pid_prepped <- sort(intersect(full_pid_prepped, keep_pid_raw))
      
      if (identical(full_pid_prepped, keep_pid_prepped)) {
        for (m in models) copy_full_to_confirmatory(tr, m$suffix)
      } else {
        for (m in models) fit_one(d_conf, tr, "confirmatory", m)
      }
    } else {
      for (m in models) fit_one(d_conf, tr, "confirmatory", m)
    }
    
    # ---- PASS 3: sensitivity (confirmatory only, primary model only) ----
    m_primary <- models[[1]]
    fit_one(d_conf, tr, "confirmatory", m_primary,
            floor_val = sd_floor_sens_low,
            dispersion = "sd")
    fit_one(d_conf, tr, "confirmatory",
            list(label = "primary", suffix = paste0("_floor", sd_floor_sens_low),
                 stan_file = m_primary$stan_file, desc = "Gaussian"),
            floor_val = sd_floor_sens_low)
    fit_one(d_conf, tr, "confirmatory",
            list(label = "primary", suffix = paste0("_floor", sd_floor_sens_high),
                 stan_file = m_primary$stan_file, desc = "Gaussian"),
            floor_val = sd_floor_sens_high)
    fit_one(d_conf, tr, "confirmatory",
            list(label = "primary", suffix = "_mad",
                 stan_file = m_primary$stan_file, desc = "Gaussian"),
            dispersion = "mad")
  }
  
  invisible(TRUE)
}