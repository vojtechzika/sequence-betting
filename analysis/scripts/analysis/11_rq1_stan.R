# ============================================================
# 11_rq1_stan.R
#
# PURPOSE
#   Fits hierarchical betting models for RQ1 (extensive margin
#   of betting). Always fits both primary and alternative models.
#   Model selection is handled downstream in rq1_tables().
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
#   stan/rq1.stan                       -- use_bb=0 Bernoulli, use_bb=1 BB
#
# OUTPUT
#   path_mod/rq1_fit_sequences_<tr>_<tag>.rds        -- primary
#   path_mod/rq1_pid_levels_<tr>_<tag>.rds
#   path_mod/rq1_seq_levels_<tr>_<tag>.rds
#   path_mod/rq1_fit_sequences_<tr>_<tag>_alt.rds    -- alternative
#   path_mod/rq1_pid_levels_<tr>_<tag>_alt.rds
#   path_mod/rq1_seq_levels_<tr>_<tag>_alt.rds
#
# TAGS
#   full          -- all participants in treatment
#   confirmatory  -- normative betters only
#
# CALL ORDER IN PIPELINE:
#   rq1_stan(cfg)         -- fits primary + alternative
#   rq1_diagnostics(cfg)  -- determines selected model
#   rq1_tables(cfg)       -- produces tables from selected model
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq1_stan <- function(cfg) {
  
  # ============================================================
  # MODEL SPECIFICATION
  #   primary     : Bernoulli logistic (use_bb = 0)
  #   alternative : Beta-Binomial     (use_bb = 1)
  #   stan file   : rq1.stan (single file, use_bb flag switches model)
  # ============================================================
  models <- list(
    list(label  = "primary",
         suffix = "",
         use_bb = 0L,
         desc   = "Bernoulli"),
    list(label  = "alternative",
         suffix = "_alt",
         use_bb = 1L,
         desc   = "Beta-Binomial")
  )
  stan_file <- here::here("stan", "rq1.stan")
  # ============================================================
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  model  <- cfg$model
  
  infile  <- file.path(path_src, "master_sequences.csv")
  f_drift <- file.path(path_mod, "drift_decisions.rds")
  
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
  
  # ---- Fit one subset (all models) ----
  fit_one <- function(d, tr, tag) {
    
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
    
    for (m in models) {
      
      f_fit <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_", tag, m$suffix, ".rds"))
      f_pid <- file.path(path_mod, paste0("rq1_pid_levels_",    tr, "_", tag, m$suffix, ".rds"))
      f_seq <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_", tag, m$suffix, ".rds"))
      
      if (should_skip(c(f_fit, f_pid, f_seq), cfg, "model",
                      paste0("RQ1 Stan (", tr, "/", tag, "/", m$label, ")"))) next
      
      data_list <- list(
        N              = length(pid_levels),
        S              = length(seq_levels),
        T              = nrow(d),
        pid            = as.integer(d$pid_i),
        sid            = as.integer(d$sid_s),
        y              = as.integer(d$y),
        use_bb         = m$use_bb,
        include_drift  = include_drift,
        block_c        = as.numeric(d$block_c),
        prior_gamma_sd = prior_gamma_sd
      )
      
      msg("RQ1 Stan: fitting tr=", tr, " tag=", tag,
          " [", m$desc, "]",
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
    }
    
    invisible(NULL)
  }
  
  # ---- Copy full to confirmatory (all models) ----
  copy_full_to_confirmatory <- function(tr) {
    for (m in models) {
      for (nm in c("rq1_fit_sequences", "rq1_pid_levels", "rq1_seq_levels")) {
        f_full <- file.path(path_mod, paste0(nm, "_", tr, "_full",         m$suffix, ".rds"))
        f_conf <- file.path(path_mod, paste0(nm, "_", tr, "_confirmatory", m$suffix, ".rds"))
        if (!file.exists(f_full)) next
        if (should_skip(f_conf, cfg, "model",
                        paste0("RQ1 Stan (", tr, "/confirmatory/",
                               m$label, " copy)"))) next
        file.copy(f_full, f_conf, overwrite = TRUE)
      }
    }
    msg("RQ1 Stan: confirmatory == full for tr=", tr, " -> copied (all models).")
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