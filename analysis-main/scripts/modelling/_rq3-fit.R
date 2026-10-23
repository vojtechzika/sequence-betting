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
#   ppc-summary.csv, ppc-<tr>-<sample>.png
#                                (hurdle_gamma and gamma only) supplementary
#                                posterior predictive check of the loss
#                                distribution, see below
#   inputs.csv, run-info.txt
#
# POSTERIOR PREDICTIVE CHECK (hurdle_gamma and gamma)
#   Losses are replicated at the observed trials from cfg$model$ppc$rq3_k
#   posterior draws (participant and sequence effects as fitted) and the
#   replicated distribution is compared with the observed one. Statistics:
#   share of zero losses and mean loss over all trials (hurdle_gamma only),
#   and the mean, SD, 90th, 95th and 99th percentiles and maximum of the
#   positive losses. The statistics are not preregistered. The p-value is
#   twice the smaller tail probability, and the check passes if
#   p >= 1 - cfg$design$rq3$ppc_cut. For the gamma likelihood the check
#   concerns the strictly positive losses to which the model is fitted.
#
# NOTES
#   - If the model is listed in cfg$run$prior_sensitivity, a second pass
#     refits it with the broadened priors (prior_broad in settings.R) and
#     writes a complete separate folder with the suffix -broad
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
  block <- gsub("-", "_", spec$folder)   # name of the block in cfg$model$stan
  for (pass in prior_passes(cfg, block)) rq3_fit_pass(cfg, spec, pass)
  invisible(NULL)
}

rq3_fit_pass <- function(cfg, spec, pass) {

  # sample functions (analysis_trials, sample_definition) are defined in 09
  if (!exists("analysis_trials", mode = "function"))
    source(here::here("scripts", "modelling", "09_analysis-sample.R"))

  run <- model_begin(paste0(spec$folder, pass$suffix))
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
  ppc_list    <- list()

  ppc_cut <- as.numeric(design$rq3$ppc_cut)
  stopifnot(length(ppc_cut) == 1L, is.finite(ppc_cut), ppc_cut > 0, ppc_cut < 1)

  # Posterior predictive check of the loss distribution (see header)
  rq3_ppc <- function(fit, d_fit, y_fit, include_drift, tr, sample) {

    hurdle <- spec$likelihood == "hurdle_gamma"
    post   <- rstan::extract(fit, pars = c("ap", "up", "bp", "shape", "gamma_drift",
                                            if (hurdle) c("a0", "u0", "b0")))
    K_all  <- length(post$ap)
    K_use  <- min(as.integer(cfg$model$ppc$rq3_k), K_all)
    set.seed(seed + 4001L)
    k_idx  <- sort(sample.int(K_all, K_use))

    pid_i <- d_fit$pid_i
    sid_s <- d_fit$sid_s
    Tn    <- nrow(d_fit)

    labels <- c(zero_share = "Share of zero losses",
                mean_all   = "Mean loss, all trials",
                mean_pos   = "Mean positive loss",
                sd_pos     = "SD of positive losses",
                q90_pos    = "90th percentile, positive losses",
                q95_pos    = "95th percentile, positive losses",
                q99_pos    = "99th percentile, positive losses",
                max_pos    = "Largest loss")
    stat_names <- c(if (hurdle) c("zero_share", "mean_all"),
                    "mean_pos", "sd_pos", "q90_pos", "q95_pos", "q99_pos", "max_pos")

    stats_of <- function(v) {
      pos <- v[v > 0]
      q   <- unname(quantile(pos, c(0.90, 0.95, 0.99)))
      c(zero_share = mean(v <= 0), mean_all = mean(v), mean_pos = mean(pos),
        sd_pos = sd(pos), q90_pos = q[1], q95_pos = q[2], q99_pos = q[3],
        max_pos = max(pos))[stat_names]
    }

    obs <- stats_of(y_fit)

    rep <- t(vapply(k_idx, function(k) {
      eta <- post$ap[k] + post$up[k, pid_i] + post$bp[k, sid_s]
      if (include_drift) eta <- eta + post$gamma_drift[k] * d_fit$block_c
      v <- rgamma(Tn, shape = post$shape[k], rate = post$shape[k] / exp(eta))
      if (hurdle) {
        pi0 <- plogis(post$a0[k] + post$u0[k, pid_i] + post$b0[k, sid_s])
        v[runif(Tn) < pi0] <- 0
      }
      stats_of(v)
    }, numeric(length(stat_names))))

    tbl <- rbindlist(lapply(stat_names, function(st_nm) {
      r    <- rep[, st_nm]
      p    <- min(1, 2 * min(mean(r >= obs[[st_nm]]), mean(r <= obs[[st_nm]])))
      rmed <- median(r)
      data.table(
        treatment = tr, sample = sample, statistic = st_nm, label = labels[[st_nm]],
        preregistered = FALSE, observed = obs[[st_nm]], rep_median = rmed,
        rep_q025 = unname(quantile(r, 0.025)), rep_q975 = unname(quantile(r, 0.975)),
        ratio = if (rmed > 0) obs[[st_nm]] / rmed else NA_real_,
        n_rep = K_use, p_value = p, threshold = 1 - ppc_cut,
        pass = round(p, 10) >= round(1 - ppc_cut, 10)
      )
    }))

    lev  <- unname(labels[stat_names])
    long <- rbindlist(lapply(stat_names, function(st_nm)
      data.table(label = labels[[st_nm]], v = rep[, st_nm])))
    long[, label := factor(label, levels = lev)]
    ann <- data.table(label = factor(tbl$label, levels = lev), v = tbl$observed,
                      txt = ifelse(tbl$p_value < 1 / K_use, sprintf("p < %.3f", 1 / K_use),
                                   sprintf("p = %.3f", tbl$p_value)))
    sample_lab <- paste0(c(m25 = "FN", m19 = "FP")[[tr]], " ", sample)
    g <- ggplot(long, aes(x = v)) +
      geom_density(fill = "#74ADD1", alpha = 0.5, colour = "#2166AC", linewidth = 0.7) +
      geom_vline(data = ann, aes(xintercept = v), colour = "#C0392B",
                 linewidth = 0.8, linetype = "dashed") +
      geom_text(data = ann, aes(label = txt), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
                size = 3, colour = "grey30", inherit.aes = FALSE) +
      facet_wrap(~label, scales = "free", nrow = 2) +
      labs(title = paste0("RQ3 ", spec$desc, ", posterior predictive check: ", sample_lab),
           x = NULL, y = "Density") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(size = 8.5),
            axis.title.y = element_text(size = 8.5),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
            plot.title = element_text(size = 11), plot.margin = margin(4, 8, 4, 4))
    ggsave(file.path(run, stem("ppc", tr, sample, "png")), g,
           width = 11, height = 6.5, dpi = 300)

    tbl
  }

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
    # prior scales of the model (settings.R: cfg$model$stan$<model>$prior)
    pr <- pass$prior
    stopifnot(!is.null(pr))
    prior_names <- switch(spec$likelihood,
      hurdle_gamma = c("a0", "ap", "su0", "sb0", "sup", "sbp"),
      gamma        = c("ap", "sup", "sbp"),
      gaussian     = c("alpha", "sigma_u", "sigma_s", "sigma"))
    for (nm in prior_names) data_list[[paste0("prior_", nm, "_sd")]] <- as.numeric(pr[[paste0(nm, "_sd")]])

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

    ppc_rows <- NULL
    if (spec$likelihood %in% c("hurdle_gamma", "gamma")) {
      ppc_tbl <- rq3_ppc(fit, d_fit, y_fit, as.logical(include_drift), tr, sample)
      ppc_list[[paste(tr, sample)]] <<- ppc_tbl
      ppc_rows <- ppc_tbl[, .(tag = paste(tr, sample, sep = "_"), treatment, sample,
                              check = paste0("ppc_", statistic, "_p"),
                              value = p_value, threshold, pass)]
    }

    rbind(conv, size, ppc_rows, fill = TRUE)
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
  if (length(ppc_list) > 0L)
    fwrite(rbindlist(ppc_list), file.path(run, "ppc-summary.csv"))

  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(diag))

  settings_run <- st
  if (pass$name != "main")
    settings_run <- modifyList(st, list(prior = pass$prior, prior_pass = pass$name))
  model_end(run, cfg, settings = settings_run,
            key = if (pass$name == "main") "rq3" else NULL)
}
