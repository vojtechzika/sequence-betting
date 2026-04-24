# ============================================================
# 11_rq1_stan.R
#
# PURPOSE
#   Fits hierarchical Bernoulli logistic regression models for
#   RQ1 (extensive margin of betting). Fits per treatment for
#   full sample and confirmatory subset.
#   Optionally fits Beta-Binomial robustness models (per preregistration)
#   if Bernoulli PPC is inadequate (robustness = TRUE).
#
# CONFIRMATORY SUBSET FOR RQ1:
#   Participants for whom betting is EU-optimal (P(a* > 0) >= tau).
#   Controlled by cfg$design$a_flags$betting_normative and
#   cfg$design$a_flags$tau. In practice: zero normative non-betters
#   exist in this dataset, so confirmatory == full and full artifacts
#   are copied to confirmatory.
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/a_star_pid_flags_<tr>.rds  -- for confirmatory subset
#   path_mod/drift_decisions.rds        -- for drift adjustment
#   stan/rq1_bets.stan                  -- Bernoulli (use_bb=0) and BB (use_bb=1)
#
# OUTPUT
#   path_mod/rq1_fit_sequences_<tr>_<tag>[_bb].rds
#   path_mod/rq1_pid_levels_<tr>_<tag>[_bb].rds
#   path_mod/rq1_seq_levels_<tr>_<tag>[_bb].rds
#
# TAGS
#   full          -- all participants in treatment
#   confirmatory  -- normative betters only
# 
# ROBUSTNESS (Beta-Binomial):
#   Called with robustness = TRUE after rq1_diagnostics() has run.
#   Checks rq1_diagnostics.csv -- if all models adequate, exits silently.
#   If any model inadequate, fits Beta-Binomial models with _bb suffix.
#   Requires stan/rq1_bets_bb.stan to exist.
#
# CALL ORDER IN PIPELINE:
#   rq1_stan(cfg)                    -- primary Bernoulli fits
#   rq1_diagnostics(cfg)             -- PPC check
#   rq1_stan(cfg, robustness = TRUE) -- BB if needed, no-op otherwise
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq1_stan <- function(cfg, robustness = FALSE) {
  
  # ---- Beta-Binomial robustness (only if PPC inadequate) ----
  if (robustness) {
    f_diag <- file.path(path_out, "rq1_diagnostics.csv")
    if (!file.exists(f_diag)) {
      msg("RQ1: diagnostics CSV not found -- skipping Beta-Binomial robustness")
      return(invisible(TRUE))
    }
    diag <- fread(f_diag)
    if (all(diag$adequate)) {
      msg("RQ1: Bernoulli PPC adequate -- Beta-Binomial robustness not required")
      return(invisible(TRUE))
    }
    msg("RQ1: Bernoulli PPC inadequate -- fitting Beta-Binomial robustness (per preregistration)")
  }
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  model  <- cfg$model
  
  suffix    <- if (robustness) "_bb" else ""
  stan_file <- here::here("stan", "rq1_bets.stan")
  infile    <- file.path(path_src, "master_sequences.csv")
  f_drift   <- file.path(path_mod, "drift_decisions.rds")
  
  stopifnot(file.exists(infile), file.exists(stan_file), file.exists(f_drift))
  
  drift_decisions <- readRDS(f_drift)
  
  st              <- model$stan$rq1
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "seq", "block") %in% names(dt)))
  
  dt[, pid     := as.character(pid)]
  dt[, treat   := as.character(treat)]
  dt[, seq     := as.character(seq)]
  dt[, y       := as.integer(!is.na(stake) & as.numeric(stake) > 0)]
  dt[, block_c := as.numeric(block) - 2.5]
  
  sm <- rstan::stan_model(stan_file)
  
  # ---- Get drift config ----
  get_drift_cfg <- function(tr, tag) {
    key <- paste(tr, tag, "bet", sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) return(list(drift = FALSE, type = "none"))
    list(drift = isTRUE(res$drift), type = res$drift_type)
  }
  
  # ---- Fit one subset ----
  fit_one <- function(d, tr, tag) {
    
    f_fit <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_", tag, suffix, ".rds"))
    f_pid <- file.path(path_mod, paste0("rq1_pid_levels_",    tr, "_", tag, suffix, ".rds"))
    f_seq <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_", tag, suffix, ".rds"))
    
    if (should_skip(c(f_fit, f_pid, f_seq), cfg, "model",
                    paste0("RQ1 Stan (", tr, "/", tag,
                           if (robustness) "/bb" else "", ")"))) {
      return(invisible(NULL))
    }
    
    if (nrow(d) == 0) return(invisible(NULL))
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    stopifnot(length(seq_levels) == 64L)
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    drift_cfg      <- get_drift_cfg(tr, tag)
    include_drift  <- as.integer(isTRUE(drift_cfg$drift))
    prior_gamma_sd <- if (include_drift) {
      as.numeric(cfg$design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else {
      0.3
    }
    
    data_list <- list(
      N              = length(pid_levels),
      S              = length(seq_levels),
      T              = nrow(d),
      pid            = as.integer(d$pid_i),
      sid            = as.integer(d$sid_s),
      y              = as.integer(d$y),
      use_bb         = as.integer(robustness),
      include_drift  = include_drift,
      block_c        = as.numeric(d$block_c),
      prior_gamma_sd = prior_gamma_sd
    )
    
    msg("RQ1 Stan: fitting tr=", tr, " tag=", tag,
        if (robustness) " [Beta-Binomial]" else " [Bernoulli]",
        " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T,
        " | drift=", drift_cfg$type)
    
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
    
    invisible(list(fit = f_fit, pid = f_pid, seq = f_seq))
  }
  
  # ---- Copy full to confirmatory ----
  copy_full_to_confirmatory <- function(tr) {
    f_fit_full <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_full", suffix, ".rds"))
    f_pid_full <- file.path(path_mod, paste0("rq1_pid_levels_",    tr, "_full", suffix, ".rds"))
    f_seq_full <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_full", suffix, ".rds"))
    f_fit_conf <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_confirmatory", suffix, ".rds"))
    f_pid_conf <- file.path(path_mod, paste0("rq1_pid_levels_",    tr, "_confirmatory", suffix, ".rds"))
    f_seq_conf <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_confirmatory", suffix, ".rds"))
    
    if (should_skip(c(f_fit_conf, f_pid_conf, f_seq_conf), cfg, "model",
                    paste0("RQ1 Stan (", tr, "/confirmatory",
                           if (robustness) "/bb" else "", " copy)"))) {
      return(invisible(NULL))
    }
    
    stopifnot(file.exists(f_fit_full), file.exists(f_pid_full), file.exists(f_seq_full))
    stopifnot(file.copy(f_fit_full, f_fit_conf, overwrite = TRUE),
              file.copy(f_pid_full, f_pid_conf, overwrite = TRUE),
              file.copy(f_seq_full, f_seq_conf, overwrite = TRUE))
    msg("RQ1 Stan: confirmatory == full for tr=", tr,
        if (robustness) " [bb]" else "", " -> copied.")
  }
  
  tau_main <- as.numeric(design$a_flags$tau[1])
  tau_nm   <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  # ---- PASS 1: full sample ----
  for (tr in tr_vec) {
    d_full <- dt[treat == tr]
    if (nrow(d_full) == 0) next
    fit_one(d_full, tr, "full")
  }
  
  # ---- PASS 2: confirmatory subset ----
  for (tr in tr_vec) {
    
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    
    f_rds_flags <- file.path(path_mod, paste0("a_star_pid_flags_", tr, ".rds"))
    stopifnot(file.exists(f_rds_flags))
    
    flags   <- readRDS(f_rds_flags)
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    
    keep_pid <- sort(unique(as.character(flags$pid_sets[[nm_keep]])))
    full_pid <- sort(unique(dt[treat == tr, pid]))
    
    if (identical(full_pid, keep_pid)) {
      copy_full_to_confirmatory(tr)
      next
    }
    
    d_conf <- dt[treat == tr & pid %in% keep_pid]
    if (nrow(d_conf) == 0) next
    fit_one(d_conf, tr, "confirmatory")
  }
  
  invisible(TRUE)
}