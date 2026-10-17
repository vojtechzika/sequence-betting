# ============================================================
# 35_rq3-stan-floor-robustness.R
#
# PURPOSE
#   Fits the RQ3 model once for each value of the numerical wealth floor
#   x_min used at the zero-wealth state, for the sensitivity analysis of the
#   welfare loss. An all-in stake (a = e) can leave zero wealth, at which CRRA
#   utility diverges for r >= 1, so the certainty equivalent and hence the
#   trial-level welfare loss
#       y_t = E_k[ max(CE(a*) - CE(a_t), 0) ] / e
#   depend on the floor at which that state is evaluated. For each floor the
#   script recomputes y_t (_rq3-welfare-loss.R) and refits the likelihood of
#   the model selected for RQ3 (hurdle-Gamma, Gamma or Gaussian, decided by the
#   folder under the key "rq3" in selected-models.csv). The baseline floor
#   (cfg$design$seq$xmin) is refitted by the same code path as the other
#   floors, so that all fits are strictly comparable; it reproduces the
#   headline estimate up to sampling error. The posterior of the grand-mean
#   welfare loss is summarised in 34_rq3-floor-robustness.R.
#
# SAMPLE
#   One treatment (default: the first of cfg$run$treatment) and one sample
#   ("confirmatory" by default), with the sample, the drift decision for the
#   outcome "stake", the posterior draws of r and a* and the seed of the RQ3
#   model fit.
#
# FUNCTION
#   rq3_stan_floor_robustness(cfg, tr, sample, floors)
#     tr       treatment (default: first treatment of the run)
#     sample   "confirmatory" or "full"
#     floors   wealth floors to evaluate (default 0.01 and 1, the smallest
#              positive payoff of the integer payoff grid)
#
# INPUT
#   model "mpl": r-draws, a-star draws, drift/drift-decisions.rds
#   model "rq3": the folder name only (determines likelihood, Stan file and
#                sampler settings)
#   scripts/modelling/_rq3-welfare-loss.R, 09_analysis-sample.R
#
# OUTPUT (model folder path_mod/rq3-floor-robustness/)
#   fit-<tr>-<sample>-xmin<floor>.rds    Stan fit for each floor
#   meta-<tr>-<sample>-xmin<floor>.rds   floor, model, sample sizes, share of
#                                        zero losses, number of all-in bets
#   diagnostics.csv                      convergence of every fit, long format
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as rq3-floor-robustness.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - The Gamma likelihood is fitted to the strictly positive losses only
#   - The model is not selectable downstream; the summary is produced by
#     34_rq3-floor-robustness.R
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

rq3_stan_floor_robustness <- function(cfg,
                                      tr     = unique(as.character(cfg$run$treatment))[1],
                                      sample = "confirmatory",
                                      floors = c(0.01, 1)) {

  stopifnot(sample %in% c("confirmatory", "full"))
  if (sample == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) {
    stop("RQ3 floor robustness: treatment ", tr, " has no confirmatory sample.", call. = FALSE)
  }

  run <- model_begin("rq3-floor-robustness")
  if (is.null(run)) return(invisible(NULL))

  if (!exists("analysis_trials", mode = "function"))
    source(here::here("scripts", "modelling", "09_analysis-sample.R"))
  source(here::here("scripts", "modelling", "_rq3-welfare-loss.R"))

  # ---- Likelihood and settings of the selected RQ3 model ----
  folder     <- basename(model_path("rq3"))
  likelihood <- c("rq3-hurdle-gamma" = "hurdle_gamma", "rq3-gamma" = "gamma",
                  "rq3-gaussian" = "gaussian")[folder]
  if (is.na(likelihood)) stop("RQ3 floor robustness: unknown model folder ", folder, call. = FALSE)
  stan_file <- here::here("scripts", "stan", paste0(folder, ".stan"))
  st        <- cfg$model$stan[[gsub("-", "_", folder)]]
  stopifnot(file.exists(stan_file), !is.null(st))

  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  e      <- as.numeric(design$seq$endowment)
  floors <- sort(unique(as.numeric(floors)))
  stopifnot(all(is.finite(floors)), all(floors > 0))
  baseline_xmin <- as.numeric(design$seq$xmin)

  # ---- Drift decision (linear drift only) ----
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  dres <- drift_decisions[[paste(tr, sample, "stake", sep = "_")]]
  if (is.null(dres)) stop("No drift decision for ", tr, "/", sample, call. = FALSE)
  if (isTRUE(dres$drift) && !identical(as.character(dres$drift_type), "linear"))
    stop("Drift type '", dres$drift_type, "' is not implemented in this model (linear only).", call. = FALSE)
  include_drift  <- as.integer(isTRUE(dres$drift))
  prior_gamma_sd <- if (include_drift) {
    as.numeric(design$drift$params[[dres$drift_type]]$prior_gamma_sd)
  } else 0.3

  sm <- rstan::stan_model(stan_file)

  fit_at_floor <- function(xmin) {

    wl <- rq3_welfare_loss(cfg, tr, sample, xmin = xmin)
    d  <- wl$d; y <- wl$y
    keep  <- if (likelihood == "gamma") y > 0 else rep(TRUE, length(y))
    d_fit <- d[which(keep)]; y_fit <- y[keep]
    stopifnot(nrow(d_fit) > 0L)

    data_list <- list(
      N = length(wl$pid_levels), S = length(wl$seq_levels), T = nrow(d_fit),
      pid = as.integer(d_fit$pid_i), sid = as.integer(d_fit$sid_s),
      y = as.vector(pmax(y_fit, 0)),
      include_drift = include_drift, block_c = as.numeric(d_fit$block_c),
      prior_gamma_sd = prior_gamma_sd
    )
    if (likelihood == "hurdle_gamma") data_list$is_zero <- as.integer(y <= 0)

    msg("RQ3 floor robustness: ", folder, " ", tr, "/", sample, " | xmin=", xmin,
        " | N=", data_list$N, " | T=", data_list$T)

    fit <- rstan::sampling(
      sm, data = data_list, iter = as.integer(st$iter), warmup = as.integer(st$warmup),
      chains = as.integer(st$chains), seed = seed,
      control = list(adapt_delta = as.numeric(st$adapt_delta),
                     max_treedepth = as.integer(st$treedepth))
    )

    nm <- paste0(tr, "-", sample, "-xmin", xmin)
    saveRDS(fit, file.path(run, paste0("fit-", nm, ".rds")))
    saveRDS(list(
      treatment = tr, sample = sample, model = folder, xmin = xmin,
      is_baseline = isTRUE(all.equal(xmin, baseline_xmin)),
      n_participants = length(wl$pid_levels), n_sequences = length(wl$seq_levels),
      n_trials = nrow(d), n_trials_fit = nrow(d_fit), n_allin_bets = sum(d$stake == e),
      share_zero_loss = mean(y <= 0)
    ), file.path(run, paste0("meta-", nm, ".rds")))

    stan_diagnostics(fit, tag = paste(tr, sample, paste0("xmin", xmin), sep = "_"),
                     treatment = tr, sample = sample, cfg = cfg,
                     max_treedepth = as.integer(st$treedepth))
  }

  diag <- rbindlist(lapply(floors, fit_at_floor))

  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, diag)
  model_end(run, cfg, settings = c(st, list(floors = floors, model = folder)))
}
