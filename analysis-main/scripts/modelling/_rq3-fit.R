# ============================================================
# _rq3-fit.R
#
# PURPOSE
#   Fitting routine shared by the three RQ3 model scripts
#   (31_rq3-stan-hurdle-gamma.R, 32_rq3-stan-gamma.R,
#   33_rq3-stan-gaussian.R). Each script supplies a specification of its
#   likelihood; this routine fits the model to every sample, checks
#   convergence and the size of the fitted data, and completes the model
#   folder. The three likelihoods share the hierarchical structure
#   (participant effects, centred sequence effects, optional linear drift)
#   and differ in how the welfare loss y is modelled.
#
# FUNCTION
#   rq3_fit(cfg, spec)
#   spec is a list with
#     folder         model folder name under data/models/
#     desc           short description for messages
#     stan_file      path of the Stan file
#     settings       sampler settings (cfg$model$stan$<block>)
#     likelihood     "hurdle_gamma", "gamma" or "gaussian"
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#   Samples rq3-full and rq3-confirmatory come from analysis_trials()
#   (09_analysis-sample.R): all trials of the participants. The drift
#   decision is that of the treatment, sample and outcome "stake".
#
# LIKELIHOODS
#   hurdle_gamma  y = 0 with probability pi0 (logit-linear in participant
#                 and sequence effects); otherwise Gamma with log-linear mean
#   gamma         Gamma on the strictly positive losses only; trials with
#                 y = 0 are not used. A robustness check for samples with few zeros
#   gaussian      Normal on y; diagnostic alternative
#
# OUTPUT (model folder path_mod/<folder>/)
#   fit-<tr>-<sample>.rds, pid-levels-<tr>-<sample>.rds,
#   seq-levels-<tr>-<sample>.rds
#   data-<tr>-<sample>.rds       fitted trials with the welfare loss y
#   diagnostics.csv              convergence and sample-size checks
#   hurdle-summary.csv           (hurdle_gamma only) model-implied mean
#                                zero-loss probability and the observed share
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as <folder>.tmp and renamed when complete;
#     an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
#   - The size check compares the participants and trials of the sample with
#     its definition before any trials are dropped (the gamma likelihood
#     fits the positive losses only; their number is stored in the data file)
#   - The hurdle summary compares the model-implied mean probability of a
#     zero loss with a threshold of 0.05, below which the hurdle is negligible.
#     It is a summary, not a pass/fail check
# ============================================================

rq3_fit <- function(cfg, spec) {

  # sample functions (analysis_trials, sample_definition) are defined in 09
  if (!exists("analysis_trials", mode = "function"))
    source(here::here("scripts", "modelling", "09_analysis-sample.R"))

  run <- model_begin(spec$folder)
  if (is.null(run)) return(invisible(NULL))

  source(here::here("scripts", "modelling", "_rq3-welfare-loss.R"))

  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  stopifnot(file.exists(spec$stan_file),
            spec$likelihood %in% c("hurdle_gamma", "gamma", "gaussian"))

  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))

  # ---- Sampler settings ----
  st              <- spec$settings
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)

  sm <- rstan::stan_model(spec$stan_file)

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

  hurdle_rows <- list()

  # ---- Fit one sample ----
  fit_one <- function(tr, sample) {

    smp <- paste0("rq3-", sample)
    wl  <- rq3_welfare_loss(cfg, tr, sample)
    d   <- wl$d
    y   <- wl$y
    pid_levels <- wl$pid_levels
    seq_levels <- wl$seq_levels

    drift_cfg      <- get_drift_cfg(tr, sample)
    include_drift  <- as.integer(drift_cfg$drift)
    prior_gamma_sd <- if (include_drift) {
      as.numeric(design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else 0.3

    # size of the sample against its definition, before any trial is dropped
    def <- sample_definition(tr, smp)
    size <- data.table(
      tag = paste(tr, sample, sep = "_"), treatment = tr, sample = sample,
      check = paste0(smp, c("_participants", "_trials")),
      value = c(length(pid_levels), nrow(d)),
      threshold = c(def$n_participants, def$n_trials)
    )
    size[, pass := value == threshold]

    # trials used by the likelihood
    keep  <- if (spec$likelihood == "gamma") y > 0 else rep(TRUE, length(y))
    d_fit <- d[which(keep)]
    y_fit <- y[keep]
    stopifnot(nrow(d_fit) > 0L)

    saveRDS(list(d = d, y = y, keep = keep, krep = wl$krep),
            file.path(run, stem("data", tr, sample)))

    data_list <- list(
      N = length(pid_levels), S = length(seq_levels), T = nrow(d_fit),
      pid = as.integer(d_fit$pid_i), sid = as.integer(d_fit$sid_s),
      y = as.vector(pmax(y_fit, 0)),
      include_drift = include_drift, block_c = as.numeric(d_fit$block_c),
      prior_gamma_sd = prior_gamma_sd
    )
    if (spec$likelihood == "hurdle_gamma") data_list$is_zero <- as.integer(y <= 0)

    msg("RQ3 ", spec$desc, ": fitting ", tr, "/", sample, " | N=", data_list$N,
        " | S=", data_list$S, " | T=", data_list$T, " | Krep=", wl$krep,
        " | drift=", drift_cfg$type)

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

    # hurdle summary: model-implied mean zero-loss probability over trials
    if (spec$likelihood == "hurdle_gamma") {
      post  <- rstan::extract(fit, pars = c("a0", "u0", "b0"))
      K_use <- min(500L, length(post$a0))
      set.seed(seed + 3001L)
      k_idx <- sort(sample.int(length(post$a0), K_use))
      pi0 <- vapply(k_idx, function(k)
        mean(plogis(post$a0[k] + post$u0[k, d$pid_i] + post$b0[k, d$sid_s])),
        numeric(1L))
      hurdle_rows[[paste(tr, sample)]] <<- data.table(
        treatment = tr, sample = sample, n_trials = nrow(d),
        obs_zero_share = mean(y <= 0),
        pi0_mean = mean(pi0), pi0_q025 = unname(quantile(pi0, 0.025)),
        pi0_q975 = unname(quantile(pi0, 0.975)),
        threshold = 0.05, negligible = mean(pi0) < 0.05
      )
    }

    rbind(conv, size)
  }

  # ---- Samples: full for each treatment, then confirmatory ----
  diag <- list()
  for (tr in tr_vec) diag[[paste(tr, "full")]] <- fit_one(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    diag[[paste(tr, "confirmatory")]] <- fit_one(tr, "confirmatory")
  }

  if (length(hurdle_rows) > 0L)
    fwrite(rbindlist(hurdle_rows), file.path(run, "hurdle-summary.csv"))

  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(diag))

  model_end(run, cfg, settings = st, key = "rq3")
}
