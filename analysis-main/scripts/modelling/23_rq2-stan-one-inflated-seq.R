# ============================================================
# 23_rq2-stan-one-inflated-seq.R
#
# PURPOSE
#   Fits an extension of the preregistered one-inflated model for RQ2 in
#   which the all-in probability also varies by sequence
#   (scripts/stan/rq2-one-inflated-seq.stan). The model of
#   22_rq2-stan-one-inflated.R lets sequences shift only the stakes that
#   are not all-in; here the all-in probability has a sequence effect c_s
#   with SD sigma_c, and sigma_c = 0 returns the preregistered model. The
#   extension is not preregistered. It is fitted to examine whether the
#   sequence effect operates on the decision to stake the whole endowment,
#   and the choice between the models is made by the criteria reported in
#   the paper (posterior of sigma_c, posterior predictive checks of the
#   sequence-level all-in share, and out-of-sample fit). The script reports
#   convergence and the checks and takes no action on them: the model used
#   for RQ2 is selected by the folder under the key "rq2" in
#   selected-models.csv.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#   Samples, drift decisions and size checks are those of
#   22_rq2-stan-one-inflated.R.
#
# POSTERIOR PREDICTIVE CHECKS
#   Computed by rq2_ppc_stats() (_rq2-ppc.R): overall all-in share
#   (preregistered statistic), all-in share in pure versus other sequences,
#   and the SD across sequences of the all-in share.
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials()), drift/drift-decisions.rds
#   scripts/stan/rq2-one-inflated-seq.stan, scripts/modelling/_rq2-ppc.R
#
# OUTPUT (model folder path_mod/rq2-one-inflated-seq/)
#   fit-<tr>-<sample>.rds        Stan fit (parameters incl. c, sigma_c)
#   pid-levels-<tr>-<sample>.rds, seq-levels-<tr>-<sample>.rds
#   data-<tr>-<sample>.rds       fitted data
#   diagnostics.csv             convergence, sample sizes, PPC p-values
#   ppc-summary.csv             PPC statistics: observed, replicate
#                               quantiles, ratio, p-value, pass
#   ppc-<tr>-<sample>.png       PPC figure
#   inputs.csv, run-info.txt
#
# NOTES
#   - If the model is listed in cfg$run$prior_sensitivity, a second pass
#     refits it with the broadened priors (prior_broad in settings.R) and
#     writes a complete separate folder with the suffix -broad
#   - The model folder is built as rq2-one-inflated-seq.tmp and renamed
#     when complete; an existing folder is never overwritten
#   - Estimands mu_a and mu_a_i have the same definition as in the other
#     RQ2 models, so the reporting scripts apply unchanged
#   - The sequence effects c_s are centred within each fit
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

rq2_stan_one_inflated_seq <- function(cfg) {
  for (pass in prior_passes(cfg, "rq2_one_inflated_seq")) rq2_stan_one_inflated_seq_pass(cfg, pass)
  invisible(NULL)
}

rq2_stan_one_inflated_seq_pass <- function(cfg, pass) {

  run <- model_begin(paste0("rq2-one-inflated-seq", pass$suffix))
  if (is.null(run)) return(invisible(NULL))

  source(here::here("scripts", "modelling", "_rq2-ppc.R"))

  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  e <- as.integer(design$seq$endowment)
  stopifnot(length(e) == 1L, e >= 2L)

  stan_file <- here::here("scripts", "stan", "rq2-one-inflated-seq.stan")
  stopifnot(file.exists(stan_file))

  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))

  # ---- Sampler settings and checks ----
  st              <- cfg$model$stan$rq2_one_inflated_seq
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)

  pr <- pass$prior
  pr_names <- c("alpha_sd", "sigma_u_sd", "sigma_s_sd", "sigma_v_sd", "sigma_c_sd",
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
      prior_sigma_c_sd     = as.numeric(pr$sigma_c_sd),
      prior_alpha_omega_sd = as.numeric(pr$alpha_omega_sd),
      prior_phi_sd         = as.numeric(pr$phi_sd)
    )

    msg("RQ2 one-inflated-seq: fitting ", tr, "/", sample, " | N=", data_list$N,
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
                         model_label = "RQ2 one-inflated model with sequence effect on all-in",
                         seq_effect = TRUE)
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

  settings_run <- st
  if (pass$name != "main")
    settings_run <- modifyList(settings_run, list(prior = pr, prior_pass = pass$name))
  model_end(run, cfg, settings = settings_run,
            key = if (pass$name == "main") "rq2" else NULL)
}
