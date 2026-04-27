# ============================================================
# 31_rq3_stan.R
#
# PURPOSE
#   Fits hierarchical models for RQ3 (certainty-equivalent welfare
#   loss from suboptimal betting). Computes trial-level welfare
#   loss y_t = E_k[max(CE(a*) - CE(a), 0)] / e using posterior
#   draws of r_i and a*(r_i). Always fits all three models.
#   Model selection handled downstream in rq3_tables().
#
# CONFIRMATORY SUBSET FOR RQ3:
#   Participants for whom betting is EU-optimal (P(a* > 0) >= tau)
#   AND MPL-consistent (if cfg$run$consistent_only == TRUE).
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/mpl_r_draws_<tr>[_consistent].rds
#   path_mod/a_star_draws_<tr>.rds
#   path_mod/a_star_pid_flags_<tr>.rds
#   path_mod/drift_decisions.rds
#   stan/rq3_primary.stan      -- hurdle-Gamma (preregistered primary)
#   stan/rq3_gamma.stan        -- Gamma-only (no hurdle; TBD)
#   stan/rq3_alternative.stan  -- Gaussian (diagnostic alternative)
#
# OUTPUT
#   path_mod/rq3_fit_sequences_<tr>_<tag>.rds         -- primary
#   path_mod/rq3_fit_sequences_<tr>_<tag>_gamma.rds   -- Gamma-only
#   path_mod/rq3_fit_sequences_<tr>_<tag>_alt.rds     -- Gaussian
#   path_mod/rq3_pid_levels_<tr>_<tag>[_gamma|_alt].rds
#   path_mod/rq3_seq_levels_<tr>_<tag>[_gamma|_alt].rds
#
# TAGS
#   full          -- all participants in treatment
#   confirmatory  -- normative betters (+ consistent if consistent_only)
#
# CALL ORDER IN PIPELINE:
#   rq3_stan(cfg)         -- fits all three models
#   rq3_diagnostics(cfg)  -- determines selected model per tr/tag
#   rq3_tables(cfg)       -- produces tables from selected model
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq3_stan <- function(cfg) {
  
  # ============================================================
  # MODEL SPECIFICATION
  #   primary     : hurdle-Gamma  (preregistered; needs is_zero)
  #   gamma_only  : Gamma-only    (no hurdle; for negligible zero mass)
  #   alternative : Gaussian      (diagnostic working model)
  # NOTE: rq3_gamma.stan not yet written -- placeholder filename
  # ============================================================
  models <- list(
    list(label     = "primary",
         suffix    = "",
         stan_file = here::here("stan", "rq3_primary.stan"),
         desc      = "hurdle-Gamma",
         is_zero   = TRUE),
    list(label     = "gamma_only",
         suffix    = "_gamma",
         stan_file = here::here("stan", "rq3_gamma.stan"),
         desc      = "Gamma-only",
         is_zero   = FALSE),
    list(label     = "alternative",
         suffix    = "_alt",
         stan_file = here::here("stan", "rq3_alternative.stan"),
         desc      = "Gaussian",
         is_zero   = FALSE)
  )
  # ============================================================
  
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec          <- unique(as.character(cfg$run$treatment))
  consistent_only <- isTRUE(cfg$run$consistent_only)
  
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ---- Design constants ----
  e     <- as.numeric(design$seq$endowment)
  xmin  <- as.numeric(design$seq$xmin)
  p_win <- as.numeric(design$seq$coin_prob)
  m_map <- design$seq$treatments
  
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  stopifnot(length(xmin) == 1L, is.finite(xmin), xmin > 0)
  stopifnot(length(p_win) == 1L, is.finite(p_win), p_win >= 0, p_win <= 1)
  
  # ---- Stan settings ----
  st              <- model$stan$rq3
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  rq3_krep        <- as.integer(model$ppc$rq3_k)
  stopifnot(rq3_krep >= 10L)
  
  # ---- Drift decisions ----
  f_drift <- file.path(path_mod, "drift_decisions.rds")
  stopifnot(file.exists(f_drift))
  drift_decisions <- readRDS(f_drift)
  
  get_drift_cfg <- function(tr, tag) {
    key <- paste(tr, tag, "stake", sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) return(list(drift = FALSE, type = "none"))
    list(drift = isTRUE(res$drift), type = res$drift_type)
  }
  
  # ---- Load master ----
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "stake", "block") %in% names(dt)))
  
  dt[, pid     := as.character(pid)]
  dt[, treat   := as.character(treat)]
  dt[, seq     := as.character(seq)]
  dt[, stake   := as.numeric(stake)]
  dt[, block_c := as.numeric(block) - 2.5]
  dt[is.na(stake), stake := 0]
  
  # ---- Consistent pids ----
  consistent_pids <- if (consistent_only) {
    scored_list <- lapply(tr_vec, function(tr) {
      f <- file.path(path_out, paste0("mpl_scored_", tr, ".csv"))
      if (!file.exists(f)) return(NULL)
      d <- fread(f, encoding = "UTF-8")
      d[inconsistent == 0L, as.character(pid)]
    })
    unique(unlist(scored_list))
  } else NULL
  
  # ---- Compute welfare loss y (shared across models) ----
  compute_y <- function(d0, tr) {
    
    stopifnot(tr %in% names(m_map))
    m <- as.numeric(m_map[[tr]])
    stopifnot(is.finite(m), m > 1)
    
    r_tag        <- if (consistent_only) paste0(tr, "_consistent") else tr
    f_r          <- file.path(path_mod, paste0("mpl_r_draws_", r_tag, ".rds"))
    stopifnot(file.exists(f_r))
    r_obj        <- readRDS(f_r)
    pid_r_levels <- as.character(r_obj$pid)
    r_draws_all  <- r_obj$r_draws
    K_all        <- nrow(r_draws_all)
    stopifnot(K_all >= 10L)
    
    f_a          <- file.path(path_mod, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_a))
    a_obj        <- readRDS(f_a)
    pid_a_levels <- as.character(a_obj$pid)
    a_star_all   <- a_obj$a_star_draws
    stopifnot(nrow(a_star_all) == K_all)
    
    Krep  <- min(rq3_krep, K_all)
    set.seed(seed)
    k_idx <- sort(sample.int(K_all, Krep, replace = FALSE))
    
    pid_levels <- sort(unique(d0$pid))
    idx_r      <- match(pid_levels, pid_r_levels)
    idx_a      <- match(pid_levels, pid_a_levels)
    if (anyNA(idx_r)) stop("RQ3: pid missing from r draws for tr='", tr, "'.")
    if (anyNA(idx_a)) stop("RQ3: pid missing from a* draws for tr='", tr, "'.")
    
    r_draws <- r_draws_all[k_idx, idx_r, drop = FALSE]
    a_draws <- a_star_all[k_idx,  idx_a, drop = FALSE]
    
    d   <- copy(d0)
    d[, pid_i := match(pid, pid_levels)]
    Tn    <- nrow(d)
    y     <- numeric(Tn)
    ii    <- d$pid_i
    a_obs <- as.numeric(d$stake)
    
    for (t in seq_len(Tn)) {
      i      <- ii[t]
      r_vec  <- r_draws[, i]
      a_opt  <- a_draws[, i]
      a_act  <- rep(a_obs[t], Krep)
      ce_opt <- ce_stake_vec_r(a = a_opt, r = r_vec, m = m, e = e,
                               p_win = p_win, xmin = xmin)
      ce_act <- ce_stake_vec_r(a = a_act, r = r_vec, m = m, e = e,
                               p_win = p_win, xmin = xmin)
      y[t]   <- mean(pmax(ce_opt - ce_act, 0)) / e
    }
    
    if (!all(is.finite(y))) stop("RQ3: non-finite y for tr='", tr, "'.")
    list(y = y, d = d, pid_levels = pid_levels, Krep = Krep)
  }
  
  # ---- Fit one subset for all models ----
  fit_one <- function(d0, tr, tag) {
    
    if (nrow(d0) == 0) return(invisible(NULL))
    
    # Compute y once, shared across all three models
    y_obj      <- compute_y(d0, tr)
    y          <- y_obj$y
    d          <- y_obj$d
    pid_levels <- y_obj$pid_levels
    Krep       <- y_obj$Krep
    seq_levels <- sort(unique(d0$seq))
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
    
    drift_cfg      <- get_drift_cfg(tr, tag)
    include_drift  <- as.integer(isTRUE(drift_cfg$drift))
    prior_gamma_sd <- if (include_drift) {
      as.numeric(cfg$design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else {
      0.3
    }
    
    for (m in models) {
      
      # Skip if Stan file doesn't exist yet (placeholder)
      if (!file.exists(m$stan_file)) {
        msg("RQ3 Stan: skipping ", m$label, " -- Stan file not found: ", m$stan_file)
        next
      }
      
      f_fit <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tag, m$suffix, ".rds"))
      f_pid <- file.path(path_mod, paste0("rq3_pid_levels_",    tr, "_", tag, m$suffix, ".rds"))
      f_seq <- file.path(path_mod, paste0("rq3_seq_levels_",    tr, "_", tag, m$suffix, ".rds"))
      
      if (should_skip(c(f_fit, f_pid, f_seq), cfg, "model",
                      paste0("RQ3 Stan (", tr, "/", tag, "/", m$label, ")"))) next
      
      data_list <- list(
        N              = length(pid_levels),
        S              = length(seq_levels),
        T              = nrow(d),
        pid            = as.integer(d$pid_i),
        sid            = as.integer(d$sid_s),
        y              = as.vector(pmax(y, 0)),
        include_drift  = include_drift,
        block_c        = as.numeric(d$block_c),
        prior_gamma_sd = prior_gamma_sd
      )
      
      if (m$is_zero) {
        data_list$is_zero <- as.integer(y <= 0)
      }
      
      msg("RQ3 Stan: fitting tr=", tr, " tag=", tag,
          " [", m$desc, "]",
          " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T,
          " | Krep=", Krep, " | drift=", drift_cfg$type)
      
      sm <- rstan::stan_model(m$stan_file)
      
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
      for (nm in c("rq3_fit_sequences", "rq3_pid_levels", "rq3_seq_levels")) {
        f_full <- file.path(path_mod, paste0(nm, "_", tr, "_full",         m$suffix, ".rds"))
        f_conf <- file.path(path_mod, paste0(nm, "_", tr, "_confirmatory", m$suffix, ".rds"))
        if (!file.exists(f_full)) next
        if (should_skip(f_conf, cfg, "model",
                        paste0("RQ3 Stan (", tr, "/confirmatory/",
                               m$label, " copy)"))) next
        file.copy(f_full, f_conf, overwrite = TRUE)
      }
    }
    msg("RQ3 Stan: confirmatory == full for tr=", tr, " -> copied (all models).")
  }
  
  tau_main <- as.numeric(design$a_flags$tau[1])
  tau_nm   <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  # ---- PASS 1: full sample ----
  for (tr in tr_vec) {
    d_full <- dt[treat == tr]
    if (!is.null(consistent_pids)) d_full <- d_full[pid %in% consistent_pids]
    if (nrow(d_full) == 0) next
    fit_one(d_full, tr, "full")
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
    if (!is.null(consistent_pids)) {
      keep_pid_raw <- intersect(keep_pid_raw, consistent_pids)
    }
    
    # Use primary suffix to check full pid set
    f_pid_full       <- file.path(path_mod, paste0("rq3_pid_levels_", tr, "_full.rds"))
    full_pid_prepped <- if (file.exists(f_pid_full)) sort(as.character(readRDS(f_pid_full))) else character(0)
    keep_pid_prepped <- sort(intersect(full_pid_prepped, keep_pid_raw))
    
    if (identical(full_pid_prepped, keep_pid_prepped)) {
      copy_full_to_confirmatory(tr)
      next
    }
    
    d_conf <- dt[treat == tr & pid %in% keep_pid_raw]
    if (nrow(d_conf) == 0) next
    fit_one(d_conf, tr, "confirmatory")
  }
  
  invisible(TRUE)
}