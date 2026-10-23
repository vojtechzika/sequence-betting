# ============================================================
# 41_rq4-stan-bernoulli.R
#
# PURPOSE
#   Fits the hierarchical Bernoulli model for RQ4 (side choice conditional
#   on betting, h = 1 if Heads is chosen) and evaluates it with convergence
#   diagnostics and three posterior predictive checks (PPC). The script
#   reports the checks and takes no action on them: whether a model is used
#   is decided by selecting its folder in selected-models.csv.
#
#   Model: logit Pr(h_ist = 1) = alpha + u_i + b_s [+ drift], priors
#   alpha ~ N(0, alpha_sd), sigma_u ~ HalfN(0, sigma_u_sd),
#   sigma_s ~ HalfN(0, sigma_s_sd) from cfg$model$stan$rq4_bernoulli$prior,
#   with participant effects u_i, sequence effects b_s (sum to zero) and an
#   optional block drift (linear or categorical, as decided by the drift
#   check for the outcome "side"). The population quantities of the
#   preregistered clarification, the baseline hbar and the sequence means
#   mu_h[s], integrate over u ~ N(0, sigma_u) with cfg$model$ppc$rq4_mc_M
#   Monte Carlo draws per posterior draw.
#
# SAMPLE
#   rq4-full for each treatment in cfg$run$treatment, from analysis_trials()
#   (09_analysis-sample.R): betting trials of all participants, with no
#   restriction on participants. Side choices do not depend on r_i, so
#   neither the confirmatory flags nor the Holt-Laury consistency apply.
#   The participants and trials of every fit are checked against the
#   definition of the sample in diagnostics.csv.
#
# POSTERIOR PREDICTIVE CHECKS
#   Each statistic is compared with replicates from the fitted model.
#   p is one-sided: P(replicate >= observed). A check passes if
#   p >= 1 - cfg$design$rq4$ppc_overdisp_cut.
#     seqvar    variance across sequences of the sequence-wise Heads rate
#               (preregistered overdispersion check; a failure triggers the
#               Beta-Binomial robustness fit); replicates use the fitted
#               effects. Because the model has a free effect for every
#               sequence, this check has little power.
#     boundary  share of participants who choose the same side in all their
#               bets (not preregistered); replicates draw new participant
#               effects u_i ~ N(0, sigma_u).
#     repeat    share of consecutive bets of a participant with the same
#               side, in trial order (not preregistered); replicates use the
#               fitted effects.
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#               (through analysis_trials()), drift/drift-decisions.rds
#   scripts/stan/rq4-bernoulli.stan
#
# OUTPUT (model folder path_mod/rq4-bernoulli/)
#   fit-<tr>-full.rds, pid-levels-<tr>-full.rds, seq-levels-<tr>-full.rds
#   diagnostics.csv              convergence, sample sizes, PPC p-values
#   ppc-summary.csv              PPC statistics: observed, replicate
#                                quantiles, p-value, pass
#   ppc-<tr>-full.png            PPC figure
#   inputs.csv, run-info.txt
#
# NOTES
#   - If the model is listed in cfg$run$prior_sensitivity, a second pass
#     refits it with the broadened priors (prior_broad in settings.R) and
#     writes a complete separate folder with the suffix -broad
#   - The model folder is built as rq4-bernoulli.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
#   - Estimands (mu_h, hbar, mu_h_i, mu_h_sample, hbar_sample) are
#     evaluated at the centred block index of zero
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

rq4_stan_bernoulli <- function(cfg) {
  for (pass in prior_passes(cfg, "rq4_bernoulli")) rq4_stan_bernoulli_pass(cfg, pass)
  invisible(NULL)
}

rq4_stan_bernoulli_pass <- function(cfg, pass) {
  
  run <- model_begin(paste0("rq4-bernoulli", pass$suffix))
  if (is.null(run)) return(invisible(NULL))
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  stan_file <- here::here("scripts", "stan", "rq4-bernoulli.stan")
  stopifnot(file.exists(stan_file))
  
  # ---- Side labels ----
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads),
            length(lab_tails) == 1L, nzchar(lab_tails), lab_heads != lab_tails)
  
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  outcome_key     <- as.character(design$drift$rq4$outcome)
  stopifnot(length(outcome_key) == 1L, nzchar(outcome_key))
  
  # ---- Sampler settings and checks ----
  st              <- cfg$model$stan$rq4_bernoulli
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  pr <- pass$prior
  pr_names <- c("alpha_sd", "sigma_u_sd", "sigma_s_sd")
  stopifnot(!is.null(pr), all(pr_names %in% names(pr)),
            all(vapply(pr[pr_names], function(x) length(x) == 1L && is.finite(x) && x > 0, logical(1))))
  
  K_ppc <- as.integer(cfg$model$ppc$rq4_k)
  mc_M  <- as.integer(cfg$model$ppc$rq4_mc_M)
  cut   <- as.numeric(design$rq4$ppc_overdisp_cut)
  stopifnot(length(K_ppc) == 1L, is.finite(K_ppc), K_ppc >= 50L,
            length(mc_M) == 1L, is.finite(mc_M), mc_M >= 100L,
            length(cut) == 1L, is.finite(cut), cut > 0, cut < 1)
  
  sm <- rstan::stan_model(stan_file)
  
  get_drift_cfg <- function(tr, sample) {
    key <- paste(tr, sample, outcome_key, sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) stop("No drift decision for key: ", key, call. = FALSE)
    type <- if (isTRUE(res$drift)) as.character(res$drift_type) else "none"
    list(drift = isTRUE(res$drift), type = type,
         type_int = switch(type, none = 0L, linear = 1L, categorical = 2L,
                           stop("Unknown drift type: ", type, call. = FALSE)))
  }
  
  stem <- function(prefix, tr, sample, ext = "rds")
    paste0(prefix, "-", tr, "-", sample, ".", ext)
  
  # ---- Posterior predictive checks for one fit ----
  ppc_one <- function(fit, d, pid_levels, drift_cfg, tr, sample) {
    
    N <- length(pid_levels); S <- length(unique(d$sid_s)); Tn <- nrow(d)
    pid_i <- d$pid_i
    n_i   <- tabulate(pid_i, N)
    n_s   <- tabulate(d$sid_s, S)
    
    pars <- c("alpha", "sigma_u", "u", "beta", "gamma_drift", "gamma_block")
    post  <- rstan::extract(fit, pars = pars)
    K_all <- length(post$alpha)
    K_use <- min(K_ppc, K_all)
    set.seed(seed + 1001L)
    k_idx <- sort(sample.int(K_all, K_use))
    
    stat_seqvar   <- function(h) var(as.numeric(rowsum(h, d$sid_s)) / n_s)
    stat_boundary <- function(h) {
      s_i <- as.numeric(rowsum(h, pid_i))
      mean(s_i == 0 | s_i == n_i)
    }
    has_prev    <- which(pid_i[-1L] == pid_i[-Tn]) + 1L
    stat_repeat <- function(h) mean(h[has_prev] == h[has_prev - 1L])
    
    obs <- c(seqvar = stat_seqvar(d$h), boundary = stat_boundary(d$h),
             repeat_ = stat_repeat(d$h))
    
    blk <- as.integer(round(d$block_c + 2.5))
    lin_pred <- function(k, u_vec) {
      eta <- post$alpha[k] + u_vec[pid_i] + post$beta[k, d$sid_s]
      if (drift_cfg$type == "linear")      eta <- eta + post$gamma_drift[k] * d$block_c
      if (drift_cfg$type == "categorical") eta <- eta + post$gamma_block[k, blk]
      eta
    }
    
    set.seed(seed + 2000L)
    rep_stats <- vapply(k_idx, function(k) {
      u_new <- rnorm(N, 0, post$sigma_u[k])
      h_new <- rbinom(Tn, 1L, plogis(lin_pred(k, u_new)))
      h_fit <- rbinom(Tn, 1L, plogis(lin_pred(k, post$u[k, ])))
      c(seqvar = stat_seqvar(h_fit), boundary = stat_boundary(h_new),
        repeat_ = stat_repeat(h_fit))
    }, numeric(3L))
    
    labels <- c(seqvar   = "Variance of sequence-wise Heads rates",
                boundary = "Share of participants choosing one side in all bets",
                repeat_  = "Share of consecutive bets on the same side")
    
    tbl <- rbindlist(lapply(names(obs), function(s) {
      r <- rep_stats[s, ]
      p <- mean(r >= obs[[s]])
      data.table(
        treatment = tr, sample = sample,
        statistic = sub("_$", "", s), label = labels[[s]],
        preregistered = (s == "seqvar"),
        observed = obs[[s]], rep_median = median(r),
        rep_q025 = unname(quantile(r, 0.025)), rep_q975 = unname(quantile(r, 0.975)),
        ratio = obs[[s]] / median(r),
        n_rep = K_use, p_value = p, threshold = 1 - cut, pass = round(p, 10) >= round(1 - cut, 10)
      )
    }))
    
    # Figure: replicated distribution, observed value (dashed), one-sided p
    df_rep <- rbindlist(lapply(names(obs), function(s)
      data.table(label = labels[[s]], v = rep_stats[s, ])))
    df_rep[, label := factor(label, levels = unname(labels))]
    p_lab <- function(p) if (p < 1 / K_use) sprintf("p < %.3f", 1 / K_use) else sprintf("p = %.3f", p)
    df_obs <- data.table(label = factor(tbl$label, levels = unname(labels)),
                         v = tbl$observed, lab = vapply(tbl$p_value, p_lab, character(1)))
    g <- ggplot(df_rep, aes(x = v)) +
      geom_density(fill = "#74ADD1", alpha = 0.5, colour = "#2166AC", linewidth = 0.7) +
      geom_vline(data = df_obs, aes(xintercept = v), colour = "#C0392B",
                 linewidth = 0.8, linetype = "dashed") +
      geom_text(data = df_obs, aes(x = -Inf, y = Inf, label = lab),
                hjust = -0.1, vjust = 1.5, size = 3, colour = "grey30") +
      facet_wrap(~ label, scales = "free", nrow = 1) +
      labs(x = NULL, y = "Density") +
      theme_classic(base_size = 10) +
      theme(axis.title.y       = element_text(size = 8.5),
            strip.text         = element_text(size = 8.5),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
            plot.margin        = margin(4, 8, 4, 4))
    ggsave(file.path(run, stem("ppc", tr, sample, "png")), g,
           width = 11, height = 3.6, dpi = 300)
    
    tbl
  }
  
  # ---- Fit one treatment ----
  fit_one <- function(tr, sample) {
    
    smp <- paste0("rq4-", sample)
    d   <- analysis_trials(tr, smp)
    
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0L)
      stop("RQ4: side values outside {", lab_heads, ", ", lab_tails, "} on betting trials (",
           tr, "): ", paste(unique(head(bad$side, 10L)), collapse = ", "), call. = FALSE)
    
    d[, h := as.integer(side == lab_heads)]
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    setorder(d, pid_i, block, pos)               # trial order within participant
    
    drift_cfg      <- get_drift_cfg(tr, sample)
    include_drift  <- as.integer(drift_cfg$drift)
    prior_gamma_sd <- if (include_drift) {
      as.numeric(design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else 0.3
    
    data_list <- list(
      N = length(pid_levels), S = length(seq_levels), T = nrow(d),
      pid = as.integer(d$pid_i), sid = as.integer(d$sid_s), h = as.integer(d$h),
      include_drift = include_drift, drift_type = drift_cfg$type_int,
      block_c = as.numeric(d$block_c), prior_gamma_sd = prior_gamma_sd,
      prior_alpha_sd = as.numeric(pr$alpha_sd),
      prior_sigma_u_sd = as.numeric(pr$sigma_u_sd),
      prior_sigma_s_sd = as.numeric(pr$sigma_s_sd),
      mc_M = mc_M
    )
    
    msg("RQ4 Bernoulli: fitting ", tr, "/", sample, " | N=", data_list$N,
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
    
    ppc <- ppc_one(fit, d, pid_levels, drift_cfg, tr, sample)
    ppc_rows <- ppc[, .(tag = paste(tr, sample, sep = "_"), treatment, sample,
                        check = paste0("ppc_", statistic, "_p"),
                        value = p_value, threshold, pass)]
    
    # size of the fitted data against the definition of the sample
    def  <- sample_definition(tr, smp)
    size <- data.table(
      tag = paste(tr, sample, sep = "_"), treatment = tr, sample = sample,
      check = paste0(smp, c("_participants", "_trials")),
      value = c(length(pid_levels), nrow(d)),
      threshold = c(def$n_participants, def$n_trials)
    )
    size[, pass := value == threshold]
    
    list(diag = rbind(conv, size, ppc_rows, fill = TRUE), ppc = ppc)
  }
  
  results <- list()
  for (tr in tr_vec) results[[tr]] <- fit_one(tr, "full")
  
  fwrite(rbindlist(lapply(results, `[[`, "ppc")), file.path(run, "ppc-summary.csv"))
  
  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(lapply(results, `[[`, "diag"), fill = TRUE))
  
  settings_run <- st
  if (pass$name != "main")
    settings_run <- modifyList(settings_run, list(prior = pr, prior_pass = pass$name))
  model_end(run, cfg, settings = settings_run,
            key = if (pass$name == "main") "rq4" else NULL)
}
