# ============================================================
# 22_rq2-stan-one-inflated.R
#
# PURPOSE
#   Fits the preregistered alternative for RQ2: a hierarchical
#   one-inflated Beta-Binomial model for stakes on the raw scale (1 to e),
#   intended for the case in which the Gaussian model of
#   21_rq2-stan-gaussian.R fails the check for boundary inflation from
#   all-in stakes. The model has the same hierarchical structure, samples
#   and drift decision as the Gaussian model, plus a participant-level
#   all-in probability and a Beta-Binomial concentration parameter
#   (scripts/stan/rq2-one-inflated.stan). The script reports convergence and
#   the preregistered check and takes no action on them: which of the two
#   RQ2 models is used is decided by selecting its folder under the key
#   "rq2" in selected-models.csv.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#   Each sample is fitted separately, also when the two coincide. The
#   samples (rq2-full, rq2-confirmatory) are the same as for the Gaussian
#   model and come from analysis_trials() (09_analysis-sample.R): betting
#   trials of participants with at least design$rq2$min_bets bets. The
#   participants and trials of every fit are checked against the definition
#   of the sample in diagnostics.csv. The drift decision is that of the
#   treatment, sample and outcome "stake".
#
# POSTERIOR PREDICTIVE CHECKS (rq2_ppc_stats() in _rq2-ppc.R)
#   Preregistered statistic: share of betting trials with an all-in stake (stake equal to
#   the endowment e). A replicated stake equals e with the participant's
#   all-in probability and is otherwise 1 + Beta-Binomial(e - 1, mu * phi,
#   (1 - mu) * phi), using the fitted participant, sequence and drift
#   effects. p = P(replicate >= observed) is one-sided: only an excess of
#   all-in stakes counts against the model. The check passes if
#   p >= 1 - cfg$design$rq2$ppc_cut. The ratio of the observed share to the
#   median replicated share is also reported. Two further statistics, not
#   preregistered, test the sequence dependence of the all-in share: the
#   all-in share in pure minus other sequences (two-sided) and the SD across
#   sequences of the all-in share (one-sided). Because the model contains a
#   free all-in probability for every participant, it reproduces the overall
#   all-in share by construction, so the check has little power.
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials()), drift/drift-decisions.rds
#   scripts/stan/rq2-one-inflated.stan
#
# OUTPUT (model folder path_mod/rq2-one-inflated/)
#   fit-<tr>-<sample>.rds        Stan fit
#   pid-levels-<tr>-<sample>.rds participant order used in the fit
#   seq-levels-<tr>-<sample>.rds sequence order used in the fit
#   data-<tr>-<sample>.rds       fitted data: stake, a*, block, ...
#   diagnostics.csv             convergence, sample sizes and PPC p-value
#                               of every fit, long format
#   ppc-summary.csv             PPC statistic: observed, replicate
#                               quantiles, ratio, p-value, pass
#   ppc-<tr>-<sample>.png       PPC figure
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as rq2-one-inflated.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
#   - The participant effects u_i are centred within each fit (sum to zero),
#     as are the sequence effects b_s. Without centring, the intercept and
#     the mean of u_i are identified only through their priors, which makes
#     sampling difficult on the raw stake scale. The Gaussian model centres
#     only b_s, because its outcome is standardised within participant.
#   - The model key "mpl" is resolved through selected-models.csv
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

rq2_stan_one_inflated <- function(cfg) {
  
  run <- model_begin("rq2-one-inflated")
  if (is.null(run)) return(invisible(NULL))
  
  source(here::here("scripts", "modelling", "_rq2-ppc.R"))
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  e <- as.integer(design$seq$endowment)
  stopifnot(length(e) == 1L, e >= 2L)
  
  stan_file <- here::here("scripts", "stan", "rq2-one-inflated.stan")
  stopifnot(file.exists(stan_file))
  
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  
  # ---- Sampler settings and checks ----
  st              <- cfg$model$stan$rq2_one_inflated
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  pr <- st$prior
  pr_names <- c("alpha_sd", "sigma_u_sd", "sigma_s_sd", "sigma_v_sd",
                "alpha_omega_sd", "phi_sd")
  stopifnot(!is.null(pr), all(pr_names %in% names(pr)),
            all(vapply(pr[pr_names],
                       function(x) length(x) == 1L && is.finite(x) && x > 0, logical(1))))
  
  K_ppc <- as.integer(cfg$model$ppc$rq2_k)
  cut   <- as.numeric(design$rq2$ppc_cut)
  stopifnot(length(K_ppc) == 1L, is.finite(K_ppc), K_ppc >= 50L,
            length(cut) == 1L, is.finite(cut), cut > 0, cut < 1)
  
  sm <- rstan::stan_model(stan_file)
  
  get_drift_cfg <- function(tr, sample) {
    key <- paste(tr, sample, "stake", sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) stop("No drift decision for key: ", key, call. = FALSE)
    # These Stan models implement the linear drift term only
    if (isTRUE(res$drift) && !identical(as.character(res$drift_type), "linear"))
      stop("Drift type '", res$drift_type, "' for ", key, " is not implemented in this model (linear only).", call. = FALSE)
    list(drift = isTRUE(res$drift), type = res$drift_type)
  }
  
  stem <- function(prefix, tr, tag, ext = "rds")
    paste0(prefix, "-", tr, "-", tag, ".", ext)
  
  # ---- Fit one sample ----
  fit_one <- function(tr, sample) {
    
    smp <- paste0("rq2-", sample)
    d   <- analysis_trials(tr, smp)[, .(pid, treat, seq, stake, block_c, a_star)]
    saveRDS(d, file.path(run, stem("data", tr, sample)))
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    y <- as.integer(round(d$stake))
    stopifnot(all(y >= 1L & y <= e))
    
    # one a* per participant, in the order of pid_levels
    a_star <- d[, .(a_star = a_star[1]), by = pid_i]
    setorder(a_star, pid_i)
    
    drift_cfg      <- get_drift_cfg(tr, sample)
    include_drift  <- as.integer(drift_cfg$drift)
    prior_gamma_sd <- if (include_drift) {
      as.numeric(design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else 0.3
    
    data_list <- list(
      e = e, N = length(pid_levels), S = length(seq_levels), T = nrow(d),
      pid = as.integer(d$pid_i), sid = as.integer(d$sid_s), y = y,
      a_star = as.vector(a_star$a_star),
      include_drift = include_drift, block_c = as.numeric(d$block_c),
      prior_gamma_sd       = prior_gamma_sd,
      prior_alpha_sd       = as.numeric(pr$alpha_sd),
      prior_sigma_u_sd     = as.numeric(pr$sigma_u_sd),
      prior_sigma_s_sd     = as.numeric(pr$sigma_s_sd),
      prior_sigma_v_sd     = as.numeric(pr$sigma_v_sd),
      prior_alpha_omega_sd = as.numeric(pr$alpha_omega_sd),
      prior_phi_sd         = as.numeric(pr$phi_sd)
    )
    
    msg("RQ2 one-inflated: fitting ", tr, "/", sample, " | N=", data_list$N,
        " | S=", data_list$S, " | T=", data_list$T, " | drift=", drift_cfg$type)
    
    fit <- rstan::sampling(
      sm, data = data_list, iter = iter_val, warmup = warmup_val,
      chains = chains_val, seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    saveRDS(fit,        file.path(run, stem("fit",        tr, sample)))
    saveRDS(pid_levels, file.path(run, stem("pid-levels", tr, sample)))
    saveRDS(seq_levels, file.path(run, stem("seq-levels", tr, sample)))
    
    conv <- stan_diagnostics(fit, tag = paste(tr, sample, sep = "_"), treatment = tr,
                             sample = sample, cfg = cfg, max_treedepth = treedepth_val)
    
    # size of the fitted data against the definition of the sample
    def <- sample_definition(tr, smp)
    size <- data.table(
      tag = paste(tr, sample, sep = "_"), treatment = tr, sample = sample,
      check = paste0(smp, c("_participants", "_trials")),
      value = c(length(pid_levels), nrow(d)),
      threshold = c(def$n_participants, def$n_trials)
    )
    size[, pass := value == threshold]
    
    ppc <- rq2_ppc_stats(fit, d, pid_levels, seq_levels, as.logical(include_drift), e,
                         K_ppc, seed, cut, tr, sample, run,
                         model_label = "RQ2 one-inflated model", seq_effect = FALSE)
    ppc_row <- ppc[, .(tag = paste(tr, sample, sep = "_"), treatment, sample,
                       check = paste0("ppc_", statistic, "_p"), value = p_value,
                       threshold, pass)]
    
    list(diag = rbind(conv, size, ppc_row, fill = TRUE), ppc = ppc)
  }
  
  # ---- Samples: full for each treatment, then confirmatory ----
  results <- list()
  
  for (tr in tr_vec) {
    results[[paste(tr, "full")]] <- fit_one(tr, "full")
  }
  
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    results[[paste(tr, "confirmatory")]] <- fit_one(tr, "confirmatory")
  }
  
  fwrite(rbindlist(lapply(results, `[[`, "ppc")), file.path(run, "ppc-summary.csv"))
  
  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(lapply(results, `[[`, "diag"), fill = TRUE))
  
  model_end(run, cfg, settings = st, key = "rq2")
}
