# ============================================================
# 54_ex1-stan-regression.R
#
# PURPOSE
#   EX1, part II: predictors of the Gamblerish-Hot-handish Index (GHI).
#   The participant-level index chi_i is regressed on standardised
#   optimism (LOT-R), standardised log mean response time over all trials,
#   and the standardised posterior mean risk parameter r:
#       chi_i = alpha + b_opt Zopt_i + b_rt Zrt_i + b_r Zr_i + e_i,
#       e_i ~ N(0, sigma),
#   with priors alpha, b_* ~ N(0, 1) and sigma ~ HalfN(0, 1) as preregistered
#   (scripts/stan/ex1-regression.stan).
#
#   Following the preregistration, uncertainty in chi_i is propagated by
#   fitting the same regression separately for Trep posterior draws of chi_i
#   and pooling the posteriors with equal weights (the same number of draws
#   from every fit). Each fit is a separate Stan run, so convergence is
#   assessed for every fit and summarised in diagnostics.csv.
#
#   A draw in which chi_i is constant across participants (in practice zero
#   for every participant, when the similarity weights of the anchors
#   coincide) has no proper posterior under the registered model: the
#   likelihood is unbounded as sigma approaches zero. Such a draw is not
#   fitted. It enters the pooled posterior with equal weight as the
#   limiting solution (alpha equal to the constant, all slopes and sigma
#   zero), and the number of such draws is reported in diagnostics.csv.
#
# SAMPLE
#   Participants of the RQ4 fit (rq4-full) with complete covariates
#   (LOT-R, response time, r). The number of complete cases is checked
#   against the number of participants in diagnostics.csv.
#
# INPUT
#   model "ex1-chi": chi-draws-<tr>-eps<main>.rds   (51_ex1-chi-draws.R)
#   model "mpl":     r-draws-<tr>.rds, sample/analysis-trials.csv
#                    (screen time on all trials, through analysis_trials())
#   path_der/lotr-scored.csv                         (LOT-R score per participant)
#   cfg$model$stan$ex1_regression (sampler settings, per fit)
#   cfg$model$simulation$ex1_trep (number of fits)
#
# OUTPUT (model folder path_mod/ex1-regression/)
#   pooled-<tr>.rds             list(pooled [draws x parameters], t_idx, N, Trep,
#                               covariates [pid, lotr_score, rt_log_mean, r_mean])
#   fit-diagnostics-<tr>.csv    Rhat, ESS and divergences of every fit, and
#                               whether the draw was degenerate (not fitted)
#   diagnostics.csv             convergence summarised over the fitted draws,
#                               the number of degenerate draws and the number
#                               of complete cases
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as ex1-regression.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - The chi draws are those at the main tolerance (first value of
#     cfg$design$ex1$eps)
#   - A fit meets the criteria if Rhat <= rhat_max, ESS >= ess_min and there
#     are no divergent transitions (cfg$model$convergence); the check reports
#     the number of fitted draws that do not. Degenerate draws are not fitted
#     and are excluded from these summaries
#   - Fit j uses seed + j; the set of draws is sampled with cfg$run$seed
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

ex1_stan_regression <- function(cfg) {

  run <- model_begin("ex1-regression")
  if (is.null(run)) return(invisible(NULL))

  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))

  st    <- cfg$model$stan$ex1_regression
  Trep0 <- as.integer(cfg$model$simulation$ex1_trep)
  conv  <- cfg$model$convergence
  stopifnot(!is.null(st), length(Trep0) == 1L, is.finite(Trep0), Trep0 >= 10L)
  keep_per_fit <- 1000L

  eps_main <- as.numeric(cfg$design$ex1$eps)[1]
  eps_nm   <- gsub("\\.", "", sprintf("%.2f", eps_main))

  stan_file <- here::here("scripts", "stan", "ex1-regression.stan")
  stopifnot(file.exists(stan_file))
  sm <- rstan::stan_model(stan_file)

  lotr <- fread(file.path(path_der, "lotr-scored.csv"), encoding = "UTF-8")
  lotr[, `:=`(pid = as.character(pid), lotr_score = as.numeric(lotr_score))]
  lotr <- lotr[, .(pid, lotr_score)]

  z_score <- function(x) {
    s <- sd(x)
    if (!is.finite(s) || s <= 0) stop("Cannot standardise: sd <= 0", call. = FALSE)
    as.vector((x - mean(x)) / s)
  }

  par_names <- c("alpha", "b_opt", "b_rt", "b_r", "sigma")
  diag <- list()

  for (tr in tr_vec) {

    chi <- readRDS(model_path("ex1-chi", paste0("chi-draws-", tr, "-eps", eps_nm, ".rds")))
    pid_levels <- as.character(chi$pid_levels)
    stopifnot(is.matrix(chi$chi_i), ncol(chi$chi_i) == length(pid_levels))

    K_all <- nrow(chi$chi_i)
    Trep  <- min(Trep0, K_all)
    set.seed(seed)
    t_idx <- sort(sample.int(K_all, Trep))
    y_rep <- chi$chi_i[t_idx, , drop = FALSE]

    # ---- Covariates aligned to the participants of the chi draws ----
    r_obj <- readRDS(model_path("mpl", paste0("r-draws-", tr, ".rds")))
    r_dt  <- data.table(pid = as.character(r_obj$pid), r_mean = colMeans(r_obj$r_draws))

    d  <- analysis_trials(tr, "rq1-full")   # all trials, betting or not
    rt <- d[is.finite(screen_ms) & screen_ms > 0, .(rt_log_mean = mean(log(screen_ms))), by = pid]

    cov <- Reduce(function(a, b) merge(a, b, by = "pid", all.x = TRUE),
                  list(data.table(pid = pid_levels), lotr, r_dt, rt))
    cov <- cov[match(pid_levels, pid)]
    ok  <- is.finite(cov$lotr_score) & is.finite(cov$r_mean) & is.finite(cov$rt_log_mean)
    cov_ok <- cov[ok]
    y_rep  <- y_rep[, ok, drop = FALSE]
    N <- nrow(cov_ok)
    if (N < 10L) stop("EX1 regression: fewer than 10 complete cases (", tr, ").", call. = FALSE)

    Zopt <- z_score(cov_ok$lotr_score)
    Zrt  <- z_score(cov_ok$rt_log_mean)
    Zr   <- z_score(cov_ok$r_mean)

    msg("EX1 regression: ", tr, " | N=", N, " | fits=", Trep)

    fit_diag <- vector("list", Trep)
    pooled   <- vector("list", Trep)

    for (j in seq_len(Trep)) {
      yj <- as.vector(y_rep[j, ])

      # Degenerate draw: chi_i is the same for every participant, so the
      # regression has no proper posterior (the likelihood is unbounded as
      # sigma -> 0). The draw is not fitted; it contributes the limiting
      # solution to the pool with the same weight as a fitted draw.
      if (diff(range(yj)) < 1e-10) {
        fit_diag[[j]] <- data.table(
          fit = j, draw_index = t_idx[j], degenerate = TRUE,
          rhat_max = NA_real_, ess_min = NA_real_, divergences = NA_integer_)
        pooled[[j]] <- matrix(0, nrow = keep_per_fit, ncol = length(par_names),
                              dimnames = list(NULL, par_names))
        pooled[[j]][, "alpha"] <- mean(yj)
        if (j %% 20L == 0L || j == Trep) msg("EX1 regression: ", tr, " fit ", j, "/", Trep)
        next
      }

      fit <- rstan::sampling(
        sm, data = list(N = N, y = yj, Zopt = Zopt, Zrt = Zrt, Zr = Zr),
        iter = st$iter, warmup = st$warmup, chains = st$chains, seed = seed + j,
        refresh = 0, show_messages = FALSE,
        control = list(adapt_delta = st$adapt_delta, max_treedepth = st$treedepth)
      )
      summ <- rstan::summary(fit, pars = par_names)$summary
      sp   <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      fit_diag[[j]] <- data.table(
        fit = j, draw_index = t_idx[j], degenerate = FALSE,
        rhat_max = max(summ[, "Rhat"], na.rm = TRUE),
        ess_min = min(summ[, "n_eff"], na.rm = TRUE),
        divergences = as.integer(sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))))

      post <- as.matrix(fit, pars = par_names)
      sel  <- unique(round(seq(1, nrow(post), length.out = keep_per_fit)))
      pooled[[j]] <- post[sel, , drop = FALSE]

      if (j %% 20L == 0L || j == Trep) msg("EX1 regression: ", tr, " fit ", j, "/", Trep)
    }

    fd <- rbindlist(fit_diag)
    fd[, meets_criteria := NA]
    fd[degenerate == FALSE, meets_criteria := rhat_max <= conv$rhat_max &
                                              ess_min >= conv$ess_min &
                                              divergences <= conv$divergences]
    fwrite(fd, file.path(run, paste0("fit-diagnostics-", tr, ".csv")))

    saveRDS(list(pooled = do.call(rbind, pooled), t_idx = t_idx, N = N, Trep = Trep,
                 covariates = cov_ok),
            file.path(run, paste0("pooled-", tr, ".rds")))

    tag <- paste(tr, "full", sep = "_")
    diag[[tr]] <- data.table(
      tag = tag, treatment = tr, sample = "full",
      check = c("rhat_max_over_fits", "ess_min_over_fits", "divergences_total",
                "fits_not_meeting_criteria", "complete_cases", "degenerate_draws"),
      value = c(max(fd$rhat_max, na.rm = TRUE), min(fd$ess_min, na.rm = TRUE),
                sum(fd$divergences, na.rm = TRUE), sum(fd$meets_criteria == FALSE, na.rm = TRUE),
                N, sum(fd$degenerate)),
      threshold = c(conv$rhat_max, conv$ess_min, conv$divergences, 0, length(pid_levels), NA)
    )
    diag[[tr]][, pass := c(value[1] <= threshold[1], value[2] >= threshold[2],
                           value[3] <= threshold[3], value[4] == threshold[4],
                           value[5] == threshold[5], TRUE)]
  }

  write_diagnostics(run, rbindlist(diag))
  model_end(run, cfg, settings = c(st, list(ex1_trep = Trep0)), key = "ex1-regression")
}
