# ============================================================
# 39_rq3_floor_robustness.R
#
# PURPOSE
#   Robustness of the RQ3 certainty-equivalent welfare-loss estimate
#   to the numerical wealth floor x_min used at the zero-wealth
#   boundary.
#
#   The losing state of an all-in stake (a = e) yields zero wealth,
#   at which CRRA utility diverges for r >= 1. The trial-level
#   welfare loss
#       y_t = E_k[ max(CE(a*) - CE(a), 0) ] / e
#   therefore depends on the floor x_min at which that state is
#   evaluated. This script recomputes y_t over a grid of floors and
#   REFITS a single, caller-specified RQ3 likelihood at each floor,
#   reporting the posterior of the grand-mean welfare loss
#   (proportion of endowment). It is the model-based companion to the
#   model-free construction check; it quantifies floor sensitivity on
#   the reported estimand itself.
#
#   The y-construction logic is identical to compute_y() in
#   31_rq3_stan.R (same Krep subsampling, seed, truncation, and
#   draw-averaging), so the baseline-floor refit reproduces the
#   production fit up to sampling noise.
#
# MODEL SELECTION
#   The fitted likelihood is a parameter (model_label). By default it
#   is read from rq3_diagnostics.csv for the given (treatment, tag),
#   so the robustness check tracks the model actually selected for the
#   headline rather than pre-committing to one. Pass model_label
#   explicitly to override.
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/mpl_r_draws_<tr>[_consistent].rds
#   path_mod/a_star_draws_<tr>.rds
#   path_mod/a_star_pid_flags_<tr>.rds
#   path_mod/drift_decisions.rds
#   path_out/rq3_diagnostics.csv          (for default model selection)
#   stan/rq3_primary.stan | rq3_gamma.stan | rq3_alternative.stan
#
# OUTPUT (all under path_out/alternatives/)
#   alt_rq3_floor_<tr>_<tag>_<model>_fit_xmin<f>.rds   -- one fit per floor
#   alt_rq3_floor_<tr>_<tag>_<model>_summary.csv       -- one row per floor:
#       posterior mean, median, and 95% CI of grand-mean mu_c,
#       plus diagnostics (n_allin, frac_zero, divergences, min_ess, max_rhat)
#
# USAGE
#   Run the pipeline through run_indices(cfg) and rq3_stan(cfg) +
#   rq3_diagnostics(cfg) (so a selected model exists), then:
#     source(here::here("scripts", "analysis", "36_rq3_floor_robustness.R"))
#     rq3_floor_robustness(cfg)                       # FN confirmatory, selected model
#     rq3_floor_robustness(cfg, model_label = "gamma_only")  # override
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq3_floor_robustness <- function(cfg,
                                 tr             = cfg$run$treatment[1],
                                 tag            = "confirmatory",
                                 floors         = c(0.01, 1.0),
                                 model_label    = NULL,
                                 ci             = c(0.025, 0.975)) {
  
  # NOTE ON FLOORS:
  #   Both the baseline floor (x_min = 0.01) and the alternative (x_min = 1.0)
  #   are refit HERE, by an identical code path, so the two rows are exactly
  #   like-for-like and the comparison carries no "were these produced the same
  #   way?" ambiguity. The baseline row will match the production headline up to
  #   MCMC noise (same model, subset, Krep, seed, drift); reporting it from this
  #   script keeps both numbers traceable to a single source.
  #
  #   x_min = 1.0 ECU is the largest floor consistent with the integer payoff
  #   grid (the smallest realizable positive payoff), hence the most
  #   conservative defensible perturbation of the baseline.
  
  design <- cfg$design
  model  <- cfg$model
  
  tr              <- as.character(tr)
  tag             <- match.arg(tag, c("confirmatory", "full"))
  consistent_only <- isTRUE(cfg$run$consistent_only)
  seed            <- as.integer(cfg$run$seed)
  
  stopifnot(length(ci) == 2L, ci[1] < ci[2], ci[1] >= 0, ci[2] <= 1)
  
  # ============================================================
  # Model registry (mirrors 31_rq3_stan.R)
  # ============================================================
  model_registry <- list(
    primary     = list(suffix = "",       stan_file = here::here("stan", "rq3_primary.stan"),
                       desc = "hurdle-Gamma", is_zero = TRUE),
    gamma_only  = list(suffix = "_gamma",  stan_file = here::here("stan", "rq3_gamma.stan"),
                       desc = "Gamma-only",   is_zero = FALSE),
    alternative = list(suffix = "_alt",    stan_file = here::here("stan", "rq3_alternative.stan"),
                       desc = "Gaussian",     is_zero = FALSE)
  )
  
  # ---- Resolve model_label: explicit > diagnostics-selected > error ----
  if (is.null(model_label)) {
    f_diag <- file.path(path_out, "rq3_diagnostics.csv")
    if (!file.exists(f_diag)) {
      stop("rq3_floor_robustness: no model_label supplied and ",
           "rq3_diagnostics.csv not found. Run rq3_diagnostics(cfg) first ",
           "or pass model_label explicitly.")
    }
    diag    <- fread(f_diag)
    tr_q    <- tr
    tag_q   <- tag
    row     <- diag[treatment == tr_q & tag == tag_q]
    if (nrow(row) != 1L || is.na(row$selected_model)) {
      stop("rq3_floor_robustness: could not resolve selected model for tr='",
           tr, "', tag='", tag, "' from rq3_diagnostics.csv. ",
           "Pass model_label explicitly.")
    }
    model_label <- as.character(row$selected_model)
    msg("rq3_floor_robustness: using diagnostics-selected model '",
        model_label, "' for tr=", tr, " tag=", tag, ".")
  }
  
  model_label <- match.arg(model_label, names(model_registry))
  mspec       <- model_registry[[model_label]]
  stopifnot(file.exists(mspec$stan_file))
  
  # ---- Design constants ----
  e     <- as.numeric(design$seq$endowment)
  p_win <- as.numeric(design$seq$coin_prob)
  m     <- as.numeric(design$seq$treatments[[tr]])
  stopifnot(is.finite(e), e > 0, is.finite(m), m > 1,
            is.finite(p_win), p_win >= 0, p_win <= 1)
  
  baseline_xmin <- as.numeric(design$seq$xmin)
  floors        <- sort(unique(as.numeric(floors)))
  stopifnot(all(is.finite(floors)), all(floors > 0), length(floors) >= 1L)
  
  # ---- Stan settings (same block as production RQ3) ----
  st              <- model$stan$rq3
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  rq3_krep        <- as.integer(model$ppc$rq3_k)
  stopifnot(rq3_krep >= 10L)
  
  # ---- Output dir ----
  path_out_alt <- file.path(path_out, "alternatives")
  dir.create(path_out_alt, showWarnings = FALSE, recursive = TRUE)
  
  fstem <- paste0("rq3_floor_robustness_", tr, "_", tag, "_", model_label)
  f_sum <- file.path(path_out, paste0(fstem, "_summary.csv"))
  
  if (should_skip(f_sum, cfg, "output",
                  paste0("RQ3 floor robustness (", tr, "/", tag, "/", model_label, ")"))) {
    return(invisible(fread(f_sum)))
  }
  
  # ============================================================
  # Build the analysis subset (mirrors 31_rq3_stan.R PASS 1/2)
  # ============================================================
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
  
  if (tag == "confirmatory") {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) {
      stop("rq3_floor_robustness: tr='", tr,
           "' is not betting-normative; no confirmatory subset.")
    }
    tau_nm  <- gsub("\\.", "", sprintf("%.2f", as.numeric(design$a_flags$tau[1])))
    flags   <- readRDS(file.path(path_mod, paste0("a_star_pid_flags_", tr, ".rds")))
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    keep_pid <- as.character(flags$pid_sets[[nm_keep]])
    
    if (consistent_only) {
      f_sc <- file.path(path_out, paste0("mpl_scored_", tr, ".csv"))
      stopifnot(file.exists(f_sc))
      cons_pid <- fread(f_sc, encoding = "UTF-8")[inconsistent == 0L, as.character(pid)]
      keep_pid <- intersect(keep_pid, cons_pid)
    }
    d0 <- dt[treat == tr & pid %in% keep_pid]
  } else {
    d0 <- dt[treat == tr]
  }
  stopifnot(nrow(d0) > 0)
  
  # ============================================================
  # Posterior draws of r and a* (same tagging as 31_rq3_stan.R)
  # ============================================================
  r_tag <- if (tag == "full") tr else if (consistent_only) paste0(tr, "_consistent") else tr
  r_obj <- readRDS(file.path(path_mod, paste0("mpl_r_draws_", r_tag, ".rds")))
  a_obj <- readRDS(file.path(path_mod, paste0("a_star_draws_", tr,
                                              if (tag == "full") "_full" else "", ".rds")))
  K_all <- nrow(r_obj$r_draws)
  stopifnot(nrow(a_obj$a_star_draws) == K_all, K_all >= 10L)
  
  Krep  <- min(rq3_krep, K_all)
  set.seed(seed)
  k_idx <- sort(sample.int(K_all, Krep, replace = FALSE))
  
  pid_levels <- sort(unique(d0$pid))
  seq_levels <- sort(unique(d0$seq))
  idx_r <- match(pid_levels, as.character(r_obj$pid))
  idx_a <- match(pid_levels, as.character(a_obj$pid))
  if (anyNA(idx_r)) stop("RQ3 floor: pid missing from r draws for tr='", tr, "'.")
  if (anyNA(idx_a)) stop("RQ3 floor: pid missing from a* draws for tr='", tr, "'.")
  
  r_draws <- r_obj$r_draws[k_idx,      idx_r, drop = FALSE]
  a_draws <- a_obj$a_star_draws[k_idx, idx_a, drop = FALSE]
  
  d <- copy(d0)
  d[, pid_i := match(pid, pid_levels)]
  d[, sid_s := match(seq, seq_levels)]
  stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
  
  Tn    <- nrow(d)
  ii    <- d$pid_i
  a_obs <- as.numeric(d$stake)
  
  # ---- Drift (read production decision; identical to 31_rq3_stan.R) ----
  drift_decisions <- readRDS(file.path(path_mod, "drift_decisions.rds"))
  dkey      <- paste(tr, tag, "stake", sep = "_")
  dcfg      <- drift_decisions[[dkey]]
  inc_drift <- as.integer(isTRUE(dcfg$drift))
  prior_gsd <- if (inc_drift) {
    as.numeric(cfg$design$drift$params[[dcfg$drift_type]]$prior_gamma_sd)
  } else 0.3
  
  # ---- Welfare loss at a given floor (verbatim compute_y logic) ----
  compute_y_floor <- function(xmin) {
    y <- numeric(Tn)
    for (t in seq_len(Tn)) {
      i      <- ii[t]
      ce_opt <- ce_stake_vec_r(a = a_draws[, i],          r = r_draws[, i],
                               m = m, e = e, p_win = p_win, xmin = xmin)
      ce_act <- ce_stake_vec_r(a = rep(a_obs[t], Krep),   r = r_draws[, i],
                               m = m, e = e, p_win = p_win, xmin = xmin)
      y[t]   <- mean(pmax(ce_opt - ce_act, 0)) / e
    }
    if (!all(is.finite(y))) stop("RQ3 floor: non-finite y at xmin=", xmin)
    y
  }
  
  # ---- Fit the chosen model at one floor; return posterior of grand mu_c ----
  sm <- rstan::stan_model(mspec$stan_file)
  
  fit_at_floor <- function(xmin) {
    
    y <- compute_y_floor(xmin)
    
    # Gamma-only is fitted on strictly positive losses (zeros are EU-optimal
    # choices and excluded); hurdle-Gamma and Gaussian use all trials.
    if (model_label == "gamma_only") {
      keep  <- y > 0
      d_fit <- d[keep]; y_fit <- y[keep]
    } else {
      d_fit <- d; y_fit <- y
    }
    stopifnot(nrow(d_fit) > 0)
    
    data_list <- list(
      N              = length(pid_levels),
      S              = length(seq_levels),
      T              = nrow(d_fit),
      pid            = as.integer(d_fit$pid_i),
      sid            = as.integer(d_fit$sid_s),
      y              = as.vector(pmax(y_fit, 0)),
      include_drift  = inc_drift,
      block_c        = as.numeric(d_fit$block_c),
      prior_gamma_sd = prior_gsd
    )
    if (isTRUE(mspec$is_zero)) data_list$is_zero <- as.integer(y <= 0)
    
    msg("RQ3 floor robustness: fit tr=", tr, " tag=", tag,
        " [", mspec$desc, "] xmin=", xmin,
        " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T,
        " | Krep=", Krep, " | drift=", if (inc_drift) dcfg$drift_type else "none")
    
    fit <- rstan::sampling(
      sm, data = data_list,
      iter = iter_val, warmup = warmup_val, chains = chains_val, seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    # Save the fit object under alternatives/ for traceability
    xtag  <- gsub("\\.", "", sprintf("%g", xmin))
    f_fit <- file.path(path_mod, paste0(fstem, "_fit_xmin", xtag, ".rds"))
    saveRDS(fit, f_fit)
    msg("Saved: ", f_fit)
    
    # Grand-mean welfare loss: average mu_c over sequences within each draw
    post <- rstan::extract(fit)
    stopifnot(!is.null(post$mu_c))
    grand <- rowMeans(post$mu_c)               # one value per posterior draw
    
    # Sampler diagnostics
    sm_summary <- rstan::summary(fit)$summary
    rhat_vec   <- sm_summary[, "Rhat"]
    ess_vec    <- sm_summary[, "n_eff"]
    n_div      <- sum(vapply(
      rstan::get_sampler_params(fit, inc_warmup = FALSE),
      function(p) sum(p[, "divergent__"]), numeric(1)
    ))
    
    data.table(
      treatment      = tr,
      tag            = tag,
      model          = model_label,
      likelihood     = mspec$desc,
      xmin           = xmin,
      is_baseline    = isTRUE(all.equal(xmin, baseline_xmin)),
      n_participants = length(pid_levels),
      n_sequences    = length(seq_levels),
      n_trials_fit   = nrow(d_fit),
      n_allin_bets   = sum(a_obs == e),
      frac_zero      = mean(y <= 0),
      grand_mean     = mean(grand),
      grand_median   = median(grand),
      grand_lo       = unname(quantile(grand, ci[1])),
      grand_hi       = unname(quantile(grand, ci[2])),
      grand_mean_pct = round(100 * mean(grand), 3),
      ci_lo_pct      = round(100 * unname(quantile(grand, ci[1])), 3),
      ci_hi_pct      = round(100 * unname(quantile(grand, ci[2])), 3),
      n_divergences  = n_div,
      max_rhat       = round(max(rhat_vec, na.rm = TRUE), 4),
      min_ess        = round(min(ess_vec,  na.rm = TRUE), 0)
    )
  }
  
  # ============================================================
  # Sweep floors, assemble summary
  # ============================================================
  out <- rbindlist(lapply(floors, fit_at_floor))
  setorder(out, xmin)
  
  fwrite(out, f_sum)
  msg("Saved: ", f_sum)
  
  span_pp <- max(out$grand_mean_pct) - min(out$grand_mean_pct)
  msg("RQ3 floor robustness [", model_label, "]: grand-mean welfare loss ",
      round(min(out$grand_mean_pct), 2), "% to ",
      round(max(out$grand_mean_pct), 2), "% (span ",
      round(span_pp, 2), " pp) across xmin in [",
      min(floors), ", ", max(floors), "].")
  
  invisible(out)
}