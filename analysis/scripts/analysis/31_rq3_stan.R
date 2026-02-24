# scripts/analysis/31_rq3_stan.R
#
# RQ3: Stan (embedded diagnostics + run-level model switch)
# PREREG-COMPLIANT propagation of r and a*:
#   y_is = E_draws[ max( CE(a*_i^(k); r_i^(k)) - CE(a_is; r_i^(k)), 0 ) ] / e
# Then fit hurdle-gamma (primary) or gaussian (diagnostic).
#
# Also enforces prereg restriction:
#   keep participants for whom betting is optimal in the FN treatment with prob >= P0
#   i.e., P(a*_i,FN > 0) >= P0, computed from a_star_draws_FN.
#
# Inputs (indices outputs):
#   models/mpl_r_draws.rds           fields: pid, r_draws        (shared across treatments)
#   models/a_star_draws_<tr>.rds     fields: pid, a_star_draws   (treatment-specific)
#
# Outputs per treatment (NO dataset name in filenames):
#   models/rq3_fit_sequences_<tr>.rds
#   models/rq3_pid_levels_<tr>.rds
#   models/rq3_seq_levels_<tr>.rds

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq3_stan <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Run-level model switch
  # ----------------------------
  rq3_model <- cfg$run$rq3_model
  if (is.null(rq3_model) || !nzchar(rq3_model)) rq3_model <- "hurdle_gamma"
  stopifnot(rq3_model %in% c("hurdle_gamma", "gaussian"))
  
  # draw subsample size (approximation of E_draws)
  rq3_krep <- cfg$run$rq3_krep
  if (is.null(rq3_krep)) rq3_krep <- 200L
  rq3_krep <- as.integer(rq3_krep)
  stopifnot(length(rq3_krep) == 1L, rq3_krep >= 10L)
  
  # ----------------------------
  # Design constants
  # ----------------------------
  stopifnot(!is.null(design$seq), !is.null(design$seq$endowment))
  e <- as.numeric(design$seq$endowment)
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  
  stopifnot(!is.null(design$seq), !is.null(design$seq$xmin))
  xmin <- as.numeric(design$seq$xmin)
  stopifnot(length(xmin) == 1L, is.finite(xmin), xmin > 0)
  
  stopifnot(!is.null(design$seq), !is.null(design$seq$coin_prob))
  p_win <- as.numeric(design$seq$coin_prob)
  stopifnot(length(p_win) == 1L, is.finite(p_win), p_win >= 0, p_win <= 1)
  
  stopifnot(!is.null(design$seq), !is.null(design$seq$treatments))
  m_map <- design$seq$treatments
  stopifnot(length(m_map) >= 1L)
  
  # ----------------------------
  # Prereg restriction: "betting is optimal in FN treatment"
  # We need:
  #   FN treatment name, and probability threshold P0.
  # Defaults:
  #   FN = treatment with max multiplier
  #   P0 = max(design$exclusion$P0) if present, else 0.90
  # Can override via cfg$run$fn_treat and cfg$run$P0_main
  # ----------------------------
  fn_treat <- cfg$run$fn_treat
  if (is.null(fn_treat) || !nzchar(fn_treat)) {
    # choose max multiplier as default FN
    fn_treat <- names(m_map)[which.max(as.numeric(unlist(m_map)))]
  }
  stopifnot(fn_treat %in% names(m_map))
  
  P0_main <- cfg$run$P0_main
  if (is.null(P0_main)) {
    if (!is.null(design$exclusion) && !is.null(design$exclusion$P0)) {
      P0_main <- max(as.numeric(design$exclusion$P0))
    } else {
      P0_main <- 0.90
    }
  }
  P0_main <- as.numeric(P0_main)
  stopifnot(length(P0_main) == 1L, is.finite(P0_main), P0_main > 0, P0_main < 1)
  
  # ----------------------------
  # Stan settings
  # ----------------------------
  stopifnot(!is.null(model$stan), !is.null(model$stan$rq3), !is.null(model$stan$rq3[[ds]]))
  st <- model$stan$rq3[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  msg("RQ3 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed,
      " | rq3_model=", rq3_model,
      " | rq3_krep=", rq3_krep,
      " | fn_treat=", fn_treat,
      " | P0_main=", sprintf("%.2f", P0_main))
  
  # ----------------------------
  # Files
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- switch(
    rq3_model,
    hurdle_gamma = here::here("stan", "rq3_hurdle_gamma.stan"),
    gaussian     = here::here("stan", "rq3_gaussian.stan")
  )
  stopifnot(file.exists(stan_file))
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Safeguard: skip existing artifacts
  # ----------------------------
  to_run <- character(0)
  for (tr in tr_vec) {
    f_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq3_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq3_seq_levels_", tr, ".rds"))
    
    if (file.exists(f_fit) || file.exists(f_pid) || file.exists(f_seq)) {
      warning(
        "RQ3 artifacts already exist for dataset='", ds, "', treatment='", tr, "'.\n",
        "Stan was NOT executed.\n",
        "To rerun, delete (as applicable):\n",
        "  ", f_fit, "\n",
        "  ", f_pid, "\n",
        "  ", f_seq
      )
    } else {
      to_run <- c(to_run, tr)
    }
  }
  
  if (length(to_run) == 0L) {
    msg("RQ3: No treatments require estimation. Nothing to run.")
    return(invisible(NULL))
  }
  
  # ----------------------------
  # Load master once
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  required <- c("pid", "treat", "seq", "stake")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) stop("master_sequences.csv missing columns: ", paste(missing, collapse = ", "))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[is.na(stake), stake := 0]
  
  # ----------------------------
  # Load shared r draws
  # ----------------------------
  f_r_all <- file.path(mod_dir, "mpl_r_draws.rds")
  stopifnot(file.exists(f_r_all))
  
  r_all <- readRDS(f_r_all)
  stopifnot(is.list(r_all), all(c("pid", "r_draws") %in% names(r_all)))
  stopifnot(is.matrix(r_all$r_draws))
  
  pid_r_levels <- as.character(r_all$pid)
  r_draws_all  <- r_all$r_draws   # K_all x N_all
  K_all <- nrow(r_draws_all)
  stopifnot(K_all >= 10L)
  
  # draw subsample for E_draws approximation
  Krep <- min(rq3_krep, K_all)
  set.seed(seed)
  k_idx <- sort(sample.int(K_all, Krep, replace = FALSE))
  
  msg("RQ3 propagation: using Krep=", Krep, " draws out of K_all=", K_all)
  
  # ----------------------------
  # Load FN a* draws once to build prereg inclusion set
  #   keep pid if P(a*_FN > 0) >= P0_main
  # ----------------------------
  f_a_fn <- file.path(mod_dir, paste0("a_star_draws_", fn_treat, ".rds"))
  stopifnot(file.exists(f_a_fn))
  
  a_fn <- readRDS(f_a_fn)
  stopifnot(is.list(a_fn), all(c("pid", "a_star_draws") %in% names(a_fn)))
  stopifnot(is.matrix(a_fn$a_star_draws))
  stopifnot(nrow(a_fn$a_star_draws) == K_all)
  
  pid_fn_levels <- as.character(a_fn$pid)
  a_fn_draws    <- a_fn$a_star_draws[k_idx, , drop = FALSE]  # Krep x N_fn
  
  p_opt_fn <- colMeans(a_fn_draws > 0)  # P(a*_FN > 0) per pid
  keep_fn <- pid_fn_levels[p_opt_fn >= P0_main]
  
  msg("RQ3 prereg restriction: keep_fn=", length(keep_fn), "/", length(pid_fn_levels),
      " participants with P(a*_FN>0) >= ", sprintf("%.2f", P0_main))
  
  # ----------------------------
  # Compile Stan once
  # ----------------------------
  sm <- rstan::stan_model(stan_file)
  
  outputs <- list()
  
  for (tr in to_run) {
    
    stopifnot(tr %in% names(m_map))
    m <- as.numeric(m_map[[tr]])
    stopifnot(is.finite(m), m > 1)
    
    # treatment-specific a* draws
    f_a <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_a))
    
    a_obj <- readRDS(f_a)
    stopifnot(is.list(a_obj), all(c("pid", "a_star_draws") %in% names(a_obj)))
    stopifnot(is.matrix(a_obj$a_star_draws))
    stopifnot(nrow(a_obj$a_star_draws) == K_all)
    
    pid_a_levels <- as.character(a_obj$pid)
    a_star_all   <- a_obj$a_star_draws
    
    # filter trials: treatment + has a* + prereg keep_fn
    d <- dt[treat == tr & pid %in% pid_a_levels & pid %in% keep_fn]
    if (nrow(d) == 0) {
      warning("RQ3: No rows for treatment='", tr, "' after prereg filtering. Skipping.")
      next
    }
    
    # levels for fit
    pid_levels_fit <- sort(unique(d$pid))
    seq_levels_fit <- sort(unique(d$seq))
    d[, pid_i := match(pid, pid_levels_fit)]
    d[, sid_s := match(seq, seq_levels_fit)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
    
    # map pid levels to draw columns
    idx_a <- match(pid_levels_fit, pid_a_levels)
    idx_r <- match(pid_levels_fit, pid_r_levels)
    if (anyNA(idx_a)) stop("RQ3: some pid_levels_fit missing from a* cache for treatment='", tr, "'.")
    if (anyNA(idx_r)) stop("RQ3: some pid_levels_fit missing from r cache (mpl_r_draws.rds).")
    
    r_draws <- r_draws_all[k_idx, idx_r, drop = FALSE]     # Krep x N_fit
    a_draws <- a_star_all[k_idx, idx_a, drop = FALSE]      # Krep x N_fit
    stopifnot(nrow(r_draws) == Krep, nrow(a_draws) == Krep)
    
    # ----------------------------------------
    # Compute y_t = E_draws[ max(CE_opt - CE_act, 0) ] / e
    # ----------------------------------------
    Tn <- nrow(d)
    y <- numeric(Tn)
    
    ii <- d$pid_i
    a_obs <- d$stake
    
    for (t in seq_len(Tn)) {
      i <- ii[t]
      
      r_vec <- r_draws[, i]
      a_opt <- a_draws[, i]
      a_act <- rep(a_obs[t], Krep)
      
      ce_opt <- ce_stake_vec_r(a = a_opt, r = r_vec, m = m, e = e, p_win = p_win, xmin = xmin)
      ce_act <- ce_stake_vec_r(a = a_act, r = r_vec, m = m, e = e, p_win = p_win, xmin = xmin)
      
      y[t] <- mean(pmax(ce_opt - ce_act, 0)) / e
    }
    
    if (!all(is.finite(y))) stop("RQ3: non-finite y encountered for treatment='", tr, "'.")
    
    # diagnostics
    zero_share <- mean(y <= 0)
    q <- stats::quantile(y, probs = c(0, .25, .5, .75, .9, .95, .99, 1), na.rm = TRUE)
    
    msg("RQ3 diagnostics (", ds, " / ", tr, "):",
        " T=", Tn,
        " | zero_share=", sprintf("%.3f", zero_share),
        " | mean=", sprintf("%.4f", mean(y)),
        " | median=", sprintf("%.4f", stats::median(y)),
        " | p95=", sprintf("%.4f", q[["95%"]]),
        " | p99=", sprintf("%.4f", q[["99%"]]),
        " | max=", sprintf("%.4f", q[["100%"]])
    )
    
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
    
    msg("Fitting RQ3 Stan for treatment=", tr,
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
    
    if (length(rstan::get_sampler_params(fit, inc_warmup = FALSE)) == 0) {
      stop("RQ3: Stan produced no samples for treatment='", tr, "'. Fit will NOT be saved.")
    }
    
    f_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq3_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq3_seq_levels_", tr, ".rds"))
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels_fit, f_pid)
    saveRDS(seq_levels_fit, f_seq)
    
    msg("Saved RQ3 fit: ", f_fit)
    msg("Saved RQ3 pid levels: ", f_pid)
    msg("Saved RQ3 seq levels: ", f_seq)
    
    outputs[[tr]] <- list(
      fit_file = f_fit,
      pid_levels_file = f_pid,
      seq_levels_file = f_seq,
      model = rq3_model,
      Krep = Krep,
      fn_treat = fn_treat,
      P0_main = P0_main
    )
  }
  
  invisible(outputs)
}

# Example:
# cfg$run$rq3_model <- "hurdle_gamma"
# cfg$run$rq3_krep  <- 200L
# cfg$run$fn_treat  <- "m25"   # optional; otherwise max multiplier
# cfg$run$P0_main   <- 0.95    # optional; otherwise max(design$exclusion$P0) or 0.90
# rq3_stan(cfg)