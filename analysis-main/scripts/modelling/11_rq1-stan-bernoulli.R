# ============================================================
# 11_rq1-stan-bernoulli.R
#
# PURPOSE
#   Fits the hierarchical Bernoulli model for RQ1 (probability of
#   betting, y = 1 if the stake is positive) and evaluates it with
#   convergence diagnostics and three posterior predictive checks
#   (PPC). The script reports the checks and takes no action on them:
#   whether a model is used is decided by selecting its folder in
#   selected-models.csv.
#
#   Model: logit Pr(y_ist = 1) = alpha + u_i + b_s [+ gamma * block_c],
#   priors alpha ~ N(0, alpha_sd), sigma_u ~ HalfN(0, sigma_u_sd),
#   sigma_s ~ HalfN(0, sigma_s_sd) from cfg$model$stan$rq1_bernoulli$prior,
#   with participant effects u_i, sequence effects b_s (sum to zero)
#   and an optional linear block drift. The drift is included when the
#   drift decision for the treatment and sample says so.
#
#   Samples: rq1-full for each treatment in cfg$run$treatment, and
#   rq1-confirmatory for treatments with
#   cfg$design$a_flags$betting_normative = TRUE (participants kept by
#   the confirmatory flags), both from analysis_trials()
#   (09_analysis-sample.R). Each sample is fitted separately, also when
#   the two coincide. The participants and trials of every fit are checked
#   against the definition of the sample in diagnostics.csv.
#
# POSTERIOR PREDICTIVE CHECKS
#   Each statistic is compared with replicates from the fitted model.
#   p is one-sided: P(replicate >= observed). A check passes if
#   p >= 1 - cfg$design$rq1$ppc_cut.
#     seqvar    variance across sequences of the sequence-wise betting
#               frequency (preregistered); replicates use the fitted
#               effects. Because the model has a free effect for every
#               sequence, this check has little power.
#     boundary  share of participants who bet in all trials or in none
#               (not preregistered); replicates draw new participant
#               effects u_i ~ N(0, sigma_u).
#     repeat    share of consecutive bets with the same outcome within
#               participant, in trial order (not preregistered);
#               replicates use the fitted effects.
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#               (through analysis_trials()), drift/drift-decisions.rds
#   scripts/stan/rq1-bernoulli.stan
#
# OUTPUT (model folder path_mod/rq1-bernoulli/)
#   fit-<tr>-<sample>.rds        Stan fit
#   pid-levels-<tr>-<sample>.rds participant order used in the fit
#   seq-levels-<tr>-<sample>.rds sequence order used in the fit
#   diagnostics.csv              convergence and PPC p-values (long format)
#   ppc-summary.csv              PPC statistics: observed, replicate
#                                quantiles, p-value, pass
#   ppc-<tr>-<sample>.png        PPC figure
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as rq1-bernoulli.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
#   - The model key "mpl" is resolved through selected-models.csv
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

rq1_stan_bernoulli <- function(cfg) {
  
  run <- model_begin("rq1-bernoulli")
  if (is.null(run)) return(invisible(NULL))
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  stan_file <- here::here("scripts", "stan", "rq1-bernoulli.stan")
  stopifnot(file.exists(stan_file))
  
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  
  # ---- Sampler settings and checks ----
  st              <- cfg$model$stan$rq1_bernoulli
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  pr <- st$prior
  stopifnot(!is.null(pr), all(c("alpha_sd", "sigma_u_sd", "sigma_s_sd") %in% names(pr)),
            all(vapply(pr[c("alpha_sd", "sigma_u_sd", "sigma_s_sd")],
                       function(x) length(x) == 1L && is.finite(x) && x > 0, logical(1))))
  
  K_ppc <- as.integer(cfg$model$ppc$rq1_k)
  cut   <- as.numeric(design$rq1$ppc_cut)
  stopifnot(length(K_ppc) == 1L, is.finite(K_ppc), K_ppc >= 50L,
            length(cut) == 1L, is.finite(cut), cut > 0, cut < 1)
  
  # ---- Data: all trials of the sample, y = 1 if a bet was placed ----
  get_sample <- function(tr, sample) {
    d <- analysis_trials(tr, paste0("rq1-", sample))
    d[, y   := bet]
    d[, pos := as.numeric(pos)]
    d
  }
  
  sm <- rstan::stan_model(stan_file)
  
  get_drift_cfg <- function(tr, sample) {
    key <- paste(tr, sample, "bet", sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) stop("No drift decision for key: ", key, call. = FALSE)
    # These Stan models implement the linear drift term only
    if (isTRUE(res$drift) && !identical(as.character(res$drift_type), "linear"))
      stop("Drift type '", res$drift_type, "' for ", key, " is not implemented in this model (linear only).", call. = FALSE)
    list(drift = isTRUE(res$drift), type = res$drift_type)
  }
  
  stem <- function(prefix, tr, sample, ext = "rds")
    paste0(prefix, "-", tr, "-", sample, ".", ext)
  
  # ---- Posterior predictive checks for one fit ----
  ppc_one <- function(fit, d, pid_levels, seq_levels, include_drift, tr, sample) {
    
    d <- copy(d)
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    setorder(d, pid_i, block, pos)          # trial order within participant
    
    N <- length(pid_levels); S <- length(seq_levels); Tn <- nrow(d)
    
    post  <- rstan::extract(fit, pars = c("alpha", "sigma_u", "u", "b", "gamma_drift"))
    alpha <- as.numeric(post$alpha)
    K_all <- length(alpha)
    K_use <- min(K_ppc, K_all)
    set.seed(seed + 1001L)
    k_idx <- sort(sample.int(K_all, K_use))
    
    pid_i <- d$pid_i
    n_i   <- tabulate(pid_i, N)
    n_s   <- tabulate(d$sid_s, S)
    
    stat_boundary <- function(y) {
      s_i <- as.numeric(rowsum(y, pid_i))
      mean(s_i == 0 | s_i == n_i)
    }
    has_prev    <- which(pid_i[-1L] == pid_i[-Tn]) + 1L
    stat_repeat <- function(y) mean(y[has_prev] == y[has_prev - 1L])
    stat_seqvar <- function(y) var(as.numeric(rowsum(y, d$sid_s)) / n_s)
    
    obs <- c(seqvar = stat_seqvar(d$y), boundary = stat_boundary(d$y),
             repeat_ = stat_repeat(d$y))
    
    lin_pred <- function(k, u_vec) {
      eta <- alpha[k] + u_vec[pid_i] + post$b[k, d$sid_s]
      if (include_drift) eta <- eta + post$gamma_drift[k] * d$block_c
      eta
    }
    
    set.seed(seed + 2000L)
    rep_stats <- vapply(k_idx, function(k) {
      u_new <- rnorm(N, 0, post$sigma_u[k])
      y_new <- rbinom(Tn, 1L, plogis(lin_pred(k, u_new)))
      y_fit <- rbinom(Tn, 1L, plogis(lin_pred(k, post$u[k, ])))
      c(seqvar = stat_seqvar(y_fit), boundary = stat_boundary(y_new),
        repeat_ = stat_repeat(y_fit))
    }, numeric(3L))
    
    labels <- c(seqvar   = "Variance of sequence-wise betting frequencies",
                boundary = "Share of participants betting in all or no trials",
                repeat_  = "Share of consecutive bets with the same outcome")
    
    tbl <- rbindlist(lapply(names(obs), function(s) {
      r <- rep_stats[s, ]
      p <- mean(r >= obs[[s]])
      data.table(
        treatment = tr, sample = sample,
        statistic = sub("_$", "", s), label = labels[[s]],
        preregistered = (s == "seqvar"),
        observed = obs[[s]], rep_median = median(r),
        rep_q025 = unname(quantile(r, 0.025)), rep_q975 = unname(quantile(r, 0.975)),
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
  
  # ---- Fit one subset; returns diagnostics and PPC table ----
  fit_one <- function(d, tr, sample) {
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    stopifnot(length(seq_levels) == 64L)
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    drift_cfg      <- get_drift_cfg(tr, sample)
    include_drift  <- as.integer(drift_cfg$drift)
    prior_gamma_sd <- if (include_drift) {
      as.numeric(design$drift$params[[drift_cfg$type]]$prior_gamma_sd)
    } else 0.3
    
    data_list <- list(
      N = length(pid_levels), S = length(seq_levels), T = nrow(d),
      pid = as.integer(d$pid_i), sid = as.integer(d$sid_s), y = as.integer(d$y),
      include_drift = include_drift, block_c = as.numeric(d$block_c),
      prior_gamma_sd = prior_gamma_sd,
      prior_alpha_sd = as.numeric(pr$alpha_sd),
      prior_sigma_u_sd = as.numeric(pr$sigma_u_sd),
      prior_sigma_s_sd = as.numeric(pr$sigma_s_sd)
    )
    
    fit <- rstan::sampling(
      sm, data = data_list, iter = iter_val, warmup = warmup_val,
      chains = chains_val, seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    saveRDS(fit,        file.path(run, stem("fit",        tr, sample)))
    saveRDS(pid_levels, file.path(run, stem("pid-levels", tr, sample)))
    saveRDS(seq_levels, file.path(run, stem("seq-levels", tr, sample)))
    
    ppc <- ppc_one(fit, d, pid_levels, seq_levels, as.logical(include_drift), tr, sample)
    
    conv <- stan_diagnostics(fit, tag = paste(tr, sample, sep = "_"), treatment = tr,
                             sample = sample, cfg = cfg, max_treedepth = treedepth_val)
    ppc_rows <- ppc[, .(tag = paste(tr, sample, sep = "_"), treatment, sample,
                        check = paste0("ppc_", statistic, "_p"),
                        value = p_value, threshold, pass)]
    
    # size of the fitted data against the definition of the sample
    smp  <- paste0("rq1-", sample)
    def  <- sample_definition(tr, smp)
    size <- data.table(
      tag = paste(tr, sample, sep = "_"), treatment = tr, sample = sample,
      check = paste0(smp, c("_participants", "_trials")),
      value = c(length(pid_levels), nrow(d)),
      threshold = c(def$n_participants, def$n_trials)
    )
    size[, pass := value == threshold]
    
    list(diag = rbind(conv, size, ppc_rows), ppc = ppc)
  }
  
  # ---- Samples: full for each treatment, then confirmatory ----
  results <- list()
  
  for (tr in tr_vec) {
    results[[paste(tr, "full")]] <- fit_one(get_sample(tr, "full"), tr, "full")
  }
  
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    results[[paste(tr, "confirmatory")]] <-
      fit_one(get_sample(tr, "confirmatory"), tr, "confirmatory")
  }
  
  fwrite(rbindlist(lapply(results, `[[`, "ppc")), file.path(run, "ppc-summary.csv"))
  
  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(lapply(results, `[[`, "diag")))
  
  model_end(run, cfg, settings = st, key = "rq1")
}
