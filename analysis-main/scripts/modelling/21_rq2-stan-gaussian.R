# ============================================================
# 21_rq2-stan-gaussian.R
#
# PURPOSE
#   Fits the hierarchical Gaussian model for RQ2 (intensive margin of
#   betting: stakes relative to the EU-optimal stake a*) and evaluates it
#   with convergence diagnostics and the preregistered posterior
#   predictive check (PPC) of boundary inflation from all-in stakes. The
#   script reports the checks and takes no action on them: whether the
#   Gaussian model or the one-inflated alternative
#   (22_rq2-stan-one-inflated.R) is used is decided by selecting its
#   folder under the key "rq2" in selected-models.csv.
#
#   Outcome: z_it = (delta_it - delta_bar_i) / s*_i, where
#   delta_it = stake_it - a*_i (ECU) on betting trials, delta_bar_i is the
#   participant mean of delta, and s*_i = max(participant dispersion of
#   delta, floor). Participants with fewer than design$rq2$min_bets betting
#   trials are dropped (20_rq2-data.R).
#   Model: z_it ~ N(alpha + u_i + b_s [+ gamma * block_c], sigma), with
#   priors from cfg$model$stan$rq2_gaussian$prior, participant effects
#   u_i, sequence effects b_s (sum to zero) and an optional linear block
#   drift, included when the drift decision for the treatment, sample and
#   outcome "stake" says so.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#   Each sample is fitted separately, also when the two coincide. The
#   samples (rq2-full, rq2-confirmatory) come from analysis_trials()
#   (09_analysis-sample.R), and the participants and trials of every fit are
#   checked against the definition of the sample in diagnostics.csv. For the
#   confirmatory sample three sensitivity fits are added, with alternative
#   lower bounds on the participant SD (design$rq2$sd_floor_sens_low and
#   _high) and with the median absolute deviation in place of the SD.
#
# POSTERIOR PREDICTIVE CHECK
#   Statistic: share of betting trials with an all-in stake (stake equal to
#   the endowment e). Replicates draw z from the fitted model (fitted
#   participant, sequence and drift effects), map it back to stakes with the
#   participant's mean stake and dispersion, round to integers and cap at
#   e. p = P(replicate >= observed) is one-sided: only an excess of all-in
#   stakes counts against the model. The check passes if
#   p >= 1 - cfg$design$rq2$ppc_cut. The ratio of the observed share to
#   the median replicated share is also reported. The check is run for the
#   primary fits (not for the sensitivity fits).
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials()), drift/drift-decisions.rds
#   scripts/stan/rq2-gaussian.stan
#
# OUTPUT (model folder path_mod/rq2-gaussian/)
#   fit-<tr>-<tag>.rds           Stan fit
#   pid-levels-<tr>-<tag>.rds    participant order used in the fit
#   seq-levels-<tr>-<tag>.rds    sequence order used in the fit
#   data-<tr>-<tag>.rds          fitted data: stake, a*, delta, z, ...
#   excluded-minbets-<tr>-<sample>.csv   participants below min_bets
#   diagnostics.csv              convergence (all fits) and PPC p-value
#                                (primary fits), long format
#   ppc-summary.csv              PPC statistic: observed, replicate
#                                quantiles, ratio, p-value, pass
#   ppc-<tr>-<sample>.png        PPC figure
#   inputs.csv, run-info.txt
#   <tag> is full, confirmatory, or confirmatory-floor<c> and
#   confirmatory-mad for the sensitivity fits.
#
# NOTES
#   - The model folder is built as rq2-gaussian.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
#   - The model key "mpl" is resolved through selected-models.csv
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

rq2_stan_gaussian <- function(cfg) {
  
  run <- model_begin("rq2-gaussian")
  if (is.null(run)) return(invisible(NULL))
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  e                  <- as.integer(design$seq$endowment)
  sd_floor           <- as.numeric(design$rq2$sd_floor)
  sd_floor_sens_low  <- as.numeric(design$rq2$sd_floor_sens_low)
  sd_floor_sens_high <- as.numeric(design$rq2$sd_floor_sens_high)
  stopifnot(length(e) == 1L, e > 0L, length(sd_floor) == 1L, sd_floor > 0)
  
  stan_file <- here::here("scripts", "stan", "rq2-gaussian.stan")
  stopifnot(file.exists(stan_file))
  
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  
  # ---- Sampler settings and checks ----
  st              <- cfg$model$stan$rq2_gaussian
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  pr <- st$prior
  pr_names <- c("alpha_sd", "sigma_u_sd", "sigma_s_sd", "sigma_sd")
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
  
  # ---- Standardized stake deviations ----
  # delta = stake - a*; z = (delta - participant mean) / max(participant
  # dispersion, floor)
  standardise <- function(d, floor_val, dispersion) {
    d <- copy(d)
    d[, delta_a := stake - a_star]
    per_pid <- d[, {
      disp <- if (dispersion == "mad") median(abs(delta_a - median(delta_a))) else sd(delta_a)
      .(n_bets = .N, delta_bar = mean(delta_a), delta_sd = disp)
    }, by = pid]
    per_pid[, sd_star := pmax(delta_sd, floor_val)]
    d <- merge(d, per_pid[, .(pid, n_bets, delta_bar, sd_star)], by = "pid")
    d[, z := (delta_a - delta_bar) / sd_star]
    stopifnot(all(is.finite(d$z)))
    d
  }
  
  # ---- Posterior predictive check of all-in stakes for one fit ----
  ppc_allin <- function(fit, d, pid_levels, seq_levels, include_drift, tr, sample) {
    
    pid_i <- match(d$pid, pid_levels)
    sid_s <- match(d$seq, seq_levels)
    stopifnot(!anyNA(pid_i), !anyNA(sid_s))
    Tn <- nrow(d)
    
    post  <- rstan::extract(fit, pars = c("alpha", "sigma", "u", "b", "gamma_drift"))
    K_all <- length(post$alpha)
    K_use <- min(K_ppc, K_all)
    set.seed(seed + 1001L)
    k_idx <- sort(sample.int(K_all, K_use))
    
    # a_star + delta_bar is the participant's mean stake
    stake_bar  <- d$a_star + d$delta_bar
    allin_obs  <- mean(d$stake >= e)
    
    set.seed(seed + 2000L)
    allin_rep <- vapply(k_idx, function(k) {
      eta <- post$alpha[k] + post$u[k, pid_i] + post$b[k, sid_s]
      if (include_drift) eta <- eta + post$gamma_drift[k] * d$block_c
      z_new     <- rnorm(Tn, eta, post$sigma[k])
      stake_new <- stake_bar + z_new * d$sd_star
      mean(round(stake_new) >= e)
    }, numeric(1L))
    
    p_allin <- mean(allin_rep >= allin_obs)
    
    # Figure: replicated distribution, observed value (dashed), one-sided p
    p_lab <- if (p_allin < 1 / K_use) sprintf("p < %.3f", 1 / K_use) else sprintf("p = %.3f", p_allin)
    sample_lab <- paste0(c(m25 = "FN", m19 = "FP")[[tr]], " ", sample)
    g <- ggplot(data.table(v = allin_rep), aes(x = v)) +
      geom_density(fill = "#74ADD1", alpha = 0.5, colour = "#2166AC", linewidth = 0.7) +
      geom_vline(xintercept = allin_obs, colour = "#C0392B",
                 linewidth = 0.8, linetype = "dashed") +
      annotate("text", x = -Inf, y = Inf, label = p_lab,
               hjust = -0.1, vjust = 1.5, size = 3, colour = "grey30") +
      labs(title = paste0("RQ2 posterior predictive check: ", sample_lab),
           x = "Share of bets that are all-in", y = "Density") +
      theme_classic(base_size = 10) +
      theme(axis.title.y       = element_text(size = 8.5),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
            plot.title         = element_text(size = 11),
            plot.margin        = margin(4, 8, 4, 4))
    ggsave(file.path(run, stem("ppc", tr, sample, "png")), g,
           width = 5, height = 3.6, dpi = 300)
    
    data.table(
      treatment = tr, sample = sample, statistic = "allin",
      label = "Share of bets that are all-in stakes", preregistered = TRUE,
      observed = allin_obs, rep_median = median(allin_rep),
      rep_q025 = unname(quantile(allin_rep, 0.025)),
      rep_q975 = unname(quantile(allin_rep, 0.975)),
      ratio = allin_obs / median(allin_rep),
      n_rep = K_use, p_value = p_allin, threshold = 1 - cut, pass = round(p_allin, 10) >= round(1 - cut, 10)
    )
  }
  
  # ---- Fit one sample (primary, or a sensitivity variant) ----
  # variant "" is the primary fit; otherwise a suffix such as "floor1" or "mad"
  fit_one <- function(tr, sample, variant = "", floor_val = sd_floor, dispersion = "sd") {
    
    tag <- if (nzchar(variant)) paste0(sample, "-", variant) else sample
    
    smp <- paste0("rq2-", sample)
    d0  <- analysis_trials(tr, smp)
    d   <- standardise(d0[, .(pid, treat, seq, stake, block_c, a_star)], floor_val, dispersion)
    
    # participants of the sample dropped by the minimum number of bets
    # (the RQ1 sample has the same participants and all their trials)
    if (!nzchar(variant)) {
      excluded <- setdiff(unique(analysis_trials(tr, paste0("rq1-", sample))$pid),
                          unique(d0$pid))
      fwrite(data.table(pid = excluded),
             file.path(run, stem("excluded-minbets", tr, sample, "csv")))
    }
    saveRDS(d, file.path(run, stem("data", tr, tag)))
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    per_pid <- unique(d[, .(pid, delta_bar, sd_star)])
    setkey(per_pid, pid)
    per_pid <- per_pid[.(pid_levels)]
    
    drift_cfg      <- get_drift_cfg(tr, sample)
    include_drift  <- as.integer(drift_cfg$drift)
    prior_gamma_sd <- if (include_drift) {
      as.numeric(design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else 0.3
    
    data_list <- list(
      N = length(pid_levels), S = length(seq_levels), T = nrow(d),
      pid = as.integer(d$pid_i), sid = as.integer(d$sid_s), z = as.vector(d$z),
      delta_bar = as.vector(per_pid$delta_bar), s_star = as.vector(per_pid$sd_star),
      e = e, include_drift = include_drift, block_c = as.numeric(d$block_c),
      prior_gamma_sd   = prior_gamma_sd,
      prior_alpha_sd   = as.numeric(pr$alpha_sd),
      prior_sigma_u_sd = as.numeric(pr$sigma_u_sd),
      prior_sigma_s_sd = as.numeric(pr$sigma_s_sd),
      prior_sigma_sd   = as.numeric(pr$sigma_sd)
    )
    
    msg("RQ2 Gaussian: fitting ", tr, "/", tag, " | N=", data_list$N,
        " | T=", data_list$T, " | drift=", drift_cfg$type)
    
    fit <- rstan::sampling(
      sm, data = data_list, iter = iter_val, warmup = warmup_val,
      chains = chains_val, seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    saveRDS(fit,        file.path(run, stem("fit",        tr, tag)))
    saveRDS(pid_levels, file.path(run, stem("pid-levels", tr, tag)))
    saveRDS(seq_levels, file.path(run, stem("seq-levels", tr, tag)))
    
    conv <- stan_diagnostics(fit, tag = paste(tr, tag, sep = "_"), treatment = tr,
                             sample = tag, cfg = cfg, max_treedepth = treedepth_val)
    
    # size of the fitted data against the definition of the sample
    def <- sample_definition(tr, smp)
    size <- data.table(
      tag = paste(tr, tag, sep = "_"), treatment = tr, sample = tag,
      check = paste0(smp, c("_participants", "_trials")),
      value = c(length(pid_levels), nrow(d)),
      threshold = c(def$n_participants, def$n_trials)
    )
    size[, pass := value == threshold]
    conv <- rbind(conv, size)
    
    ppc <- NULL
    if (!nzchar(variant)) {
      ppc <- ppc_allin(fit, d, pid_levels, seq_levels, as.logical(include_drift), tr, sample)
      conv <- rbind(conv, ppc[, .(tag = paste(tr, tag, sep = "_"), treatment, sample,
                                  check = "ppc_allin_p", value = p_value,
                                  threshold, pass)])
    }
    list(diag = conv, ppc = ppc)
  }
  
  # ---- Samples: full for each treatment, then confirmatory with sensitivity fits ----
  results <- list()
  
  for (tr in tr_vec) {
    results[[paste(tr, "full")]] <- fit_one(tr, "full")
  }
  
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    results[[paste(tr, "confirmatory")]] <- fit_one(tr, "confirmatory")
    
    # Sensitivity (confirmatory only): alternative SD floors and MAD
    for (fl in c(sd_floor_sens_low, sd_floor_sens_high)) {
      results[[paste(tr, "floor", fl)]] <-
        fit_one(tr, "confirmatory", variant = paste0("floor", fl), floor_val = fl)
    }
    results[[paste(tr, "mad")]] <- fit_one(tr, "confirmatory", variant = "mad",
                                           dispersion = "mad")
  }
  
  ppc_all <- rbindlist(lapply(results, `[[`, "ppc"))
  fwrite(ppc_all, file.path(run, "ppc-summary.csv"))
  
  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(lapply(results, `[[`, "diag")))
  
  model_end(run, cfg, settings = c(st, list(sd_floor = sd_floor,
                                            sd_floor_sens = c(sd_floor_sens_low,
                                                              sd_floor_sens_high))),
            key = "rq2")
}
