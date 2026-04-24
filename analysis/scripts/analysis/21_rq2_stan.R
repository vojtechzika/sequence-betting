# ============================================================
# 21_rq2_stan.R
#
# PURPOSE
#   Fits hierarchical Gaussian models for RQ2 (intensive margin
#   of betting). Fits per treatment for full sample and
#   confirmatory subset. Runs sensitivity analyses (alternative
#   floors and MAD dispersion) for confirmatory subset only.
#   Optionally fits one-inflated Beta-Binomial robustness model
#   if PPC flags boundary inflation.
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
#   stan/rq2_stakes.stan
#
# OUTPUT
#   path_mod/rq2_fit_sequences_<tr>_<tag>[_bb|_floor*|_mad].rds
#   path_mod/rq2_pid_levels_<tr>_<tag>[_bb|_floor*|_mad].rds
#   path_mod/rq2_seq_levels_<tr>_<tag>[_bb|_floor*|_mad].rds
#   path_mod/rq2_prepared_<tr>_<tag>[_bb|_floor*|_mad].rds
#   path_out/rq2_excluded_minbets_<tr>_<tag>[_bb|_floor*|_mad].csv
#
# TAGS
#   full          -- all participants in treatment
#   confirmatory  -- normative betters only
#
# SENSITIVITY (confirmatory only):
#   _floor<n> -- alternative sd_floor
#   _mad      -- MAD dispersion measure
#
# ROBUSTNESS:
#   Called with robustness = TRUE after rq2_diagnostics() has run.
#   Checks rq2_diagnostics.csv -- exits if adequate.
#
# CALL ORDER IN PIPELINE:
#   rq2_stan(cfg)                    -- primary + sensitivity
#   rq2_diagnostics(cfg)             -- PPC check
#   rq2_stan(cfg, robustness = TRUE) -- BB if needed, no-op otherwise
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq2_stan <- function(cfg, robustness = FALSE) {
  
  # ---- BB robustness early exit ----
  if (robustness) {
    f_diag <- file.path(path_out, "rq2_diagnostics.csv")
    if (!file.exists(f_diag)) {
      msg("RQ2: diagnostics CSV not found -- skipping BB robustness")
      return(invisible(TRUE))
    }
    diag <- fread(f_diag)
    if (all(diag$adequate)) {
      msg("RQ2: PPC adequate -- BB robustness not required")
      return(invisible(TRUE))
    }
    msg("RQ2: PPC inadequate -- fitting BB robustness (per preregistration)")
  }
  
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
  
  infile    <- file.path(path_src, "master_sequences.csv")
  stan_file <- here::here("stan", if (robustness) "rq2_stakes_bb.stan" else "rq2_stakes.stan")
  f_drift   <- file.path(path_mod, "drift_decisions.rds")
  
  stopifnot(file.exists(infile), file.exists(stan_file), file.exists(f_drift))
  
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
  
  sm <- rstan::stan_model(stan_file)
  
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
  
  # ---- Fit one subset ----
  fit_one <- function(d0, tr, tag, suffix = "", floor_val = sd_floor,
                      dispersion = "sd") {
    
    f_fit  <- file.path(path_mod, paste0("rq2_fit_sequences_", tr, "_", tag, suffix, ".rds"))
    f_pid  <- file.path(path_mod, paste0("rq2_pid_levels_",    tr, "_", tag, suffix, ".rds"))
    f_seq  <- file.path(path_mod, paste0("rq2_seq_levels_",    tr, "_", tag, suffix, ".rds"))
    f_prep <- file.path(path_mod, paste0("rq2_prepared_",      tr, "_", tag, suffix, ".rds"))
    f_excl <- file.path(path_out, paste0("rq2_excluded_minbets_", tr, "_", tag, suffix, ".csv"))
    
    skip_model <- should_skip(c(f_fit, f_pid, f_seq), cfg, "model",
                              paste0("RQ2 Stan (", tr, "/", tag, suffix, ")"))
    skip_prep  <- should_skip(f_prep, cfg, "output",
                              paste0("RQ2 prepared (", tr, "/", tag, suffix, ")"))
    
    if (skip_model && skip_prep) return(invisible(NULL))
    if (nrow(d0) == 0) return(invisible(NULL))
    
    d <- prepare_z(d0, tr, floor_val, dispersion)
    if (is.null(d) || nrow(d) == 0) return(invisible(NULL))
    
    # Save excluded participants
    excl_pids <- setdiff(unique(d0$pid), unique(d$pid))
    if (!should_skip(f_excl, cfg, "output",
                     paste0("RQ2 excluded (", tr, "/", tag, suffix, ")"))) {
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
      
      if (!robustness) {
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
          use_bb         = 0L,
          include_drift  = include_drift,
          block_c        = as.numeric(d$block_c),
          prior_gamma_sd = prior_gamma_sd
        )
      } else {
        # BB model uses raw stakes and a_star instead of z
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
      
      msg("RQ2 Stan: fitting tr=", tr, " tag=", tag, suffix,
          if (robustness) " [BB]" else " [Gaussian]",
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
  
  # ---- Copy full to confirmatory ----
  copy_full_to_confirmatory <- function(tr, suffix = "") {
    for (nm in c("rq2_fit_sequences", "rq2_pid_levels", "rq2_seq_levels", "rq2_prepared")) {
      f_full <- file.path(path_mod, paste0(nm, "_", tr, "_full", suffix, ".rds"))
      f_conf <- file.path(path_mod, paste0(nm, "_", tr, "_confirmatory", suffix, ".rds"))
      if (file.exists(f_full)) file.copy(f_full, f_conf, overwrite = TRUE)
    }
    msg("RQ2 Stan: confirmatory == full for tr=", tr, suffix, " -> copied.")
  }
  
  tau_main <- as.numeric(design$a_flags$tau[1])
  tau_nm   <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  # ---- PASS 1: full sample ----
  if (!robustness) {
    for (tr in tr_vec) {
      d_full <- dt[treat == tr]
      if (nrow(d_full) == 0) next
      fit_one(d_full, tr, "full")
    }
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
    
    suffix_primary <- if (robustness) "_bb" else ""
    
    if (file.exists(f_pid_full)) {
      full_pid_prepped <- sort(as.character(readRDS(f_pid_full)))
      keep_pid_prepped <- sort(intersect(full_pid_prepped, keep_pid_raw))
      
      if (identical(full_pid_prepped, keep_pid_prepped)) {
        copy_full_to_confirmatory(tr, suffix_primary)
      } else {
        fit_one(d_conf, tr, "confirmatory", suffix = suffix_primary)
      }
    } else {
      fit_one(d_conf, tr, "confirmatory", suffix = suffix_primary)
    }
    
    # ---- PASS 3: sensitivity (confirmatory only) ----
    fit_one(d_conf, tr, "confirmatory",
            suffix    = paste0("_floor", sd_floor_sens_low),
            floor_val = sd_floor_sens_low)
    fit_one(d_conf, tr, "confirmatory",
            suffix    = paste0("_floor", sd_floor_sens_high),
            floor_val = sd_floor_sens_high)
    fit_one(d_conf, tr, "confirmatory",
            suffix    = "_mad", dispersion = "mad")
  }
  
  invisible(TRUE)
}