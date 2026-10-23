# ============================================================
# 61_ex2-stan-associations.R
#
# PURPOSE
#   EX2: association of dispositional optimism (LOT-R) and response time
#   with the participant-level outcomes of RQ1-RQ4, estimated in one
#   hierarchical model per treatment with partial pooling of the
#   outcome-specific slopes (scripts/stan/ex2-associations.stan). The
#   outcomes are the posterior participant means of
#       b  betting probability          mu_b_i  (RQ1)
#       a  proportional stake deviation mu_a_i  (RQ2)
#       c  welfare loss                 mu_c_i  (RQ3)
#       h  Heads-choice probability     mu_h_i  (RQ4)
#   Uncertainty in the outcomes is propagated by a multiple-imputation
#   likelihood: for every participant the likelihood averages over Trep
#   posterior draws of the outcome, each standardised across participants
#   within the draw, inside a single fit (see the Stan file).
#
# SAMPLE
#   The full sample of each treatment. An outcome block contains the
#   participants of the corresponding RQ fit ("full"): all participants for
#   RQ1 and RQ3, participants with enough betting trials for RQ2
#   (cfg$design$rq2$min_bets) and participants with a bet for RQ4. A block
#   is dropped for the participants without a finite predictor; the
#   number of participants per block is checked against cfg$design$ex2$min_Nk.
#
# PREDICTORS
#   Zopt  LOT-R score standardised across all participants of both
#         treatments, so that it is on one scale in both fits
#   Zrt   mean log screen time over all trials of a participant (betting and
#         non-betting), standardised within the treatment over the
#         participants of RQ1
#
# INPUT
#   models "rq1", "rq2", "rq3", "rq4": fit-<tr>-full.rds, pid-levels-<tr>-full.rds
#                                      (parameters mu_b_i, mu_a_i, mu_c_i, mu_h_i)
#   model "mpl": sample/analysis-trials.csv, through analysis_trials()
#   path_der/lotr-scored.csv
#   cfg$model$stan$ex2 (sampler), cfg$model$simulation$ex2_trep (Trep),
#   cfg$design$ex2$min_Nk
#
# OUTPUT (model folder path_mod/ex2-associations/)
#   fit-<tr>.rds        Stan fit
#   data-<tr>.rds       list(blocks [pid and outcome per block], Nobs, Trep,
#                       draw_idx per block)
#   diagnostics.csv     convergence of each fit and sample sizes
#   inputs.csv, run-info.txt
#
# NOTES
#   - If the model is listed in cfg$run$prior_sensitivity, a second pass
#     refits it with the broadened priors (prior_broad in settings.R) and
#     writes a complete separate folder with the suffix -broad
#   - The model folder is built as ex2-associations.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - The outcome index of the Stan model is fixed (b=1, a=2, c=3, h=4)
#   - Draws are subsampled per block with seed + block index, then the
#     common number of draws (the minimum over blocks) is used
#   - Zopt is not re-standardised within blocks
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

ex2_stan_associations <- function(cfg) {
  for (pass in prior_passes(cfg, "ex2")) ex2_stan_associations_pass(cfg, pass)
  invisible(NULL)
}

ex2_stan_associations_pass <- function(cfg, pass) {

  run <- model_begin(paste0("ex2-associations", pass$suffix))
  if (is.null(run)) return(invisible(NULL))

  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))

  st      <- cfg$model$stan$ex2
  trep0   <- as.integer(cfg$model$simulation$ex2_trep)
  min_Nk  <- as.integer(cfg$design$ex2$min_Nk)
  stopifnot(!is.null(st), length(trep0) == 1L, trep0 >= 10L, length(min_Nk) == 1L, min_Nk >= 1L)

  stan_file <- here::here("scripts", "stan", "ex2-associations.stan")
  stopifnot(file.exists(stan_file))
  sm <- rstan::stan_model(stan_file)

  # ---- LOT-R, standardised across both treatments ----
  lotr <- fread(file.path(path_der, "lotr-scored.csv"), encoding = "UTF-8")
  lotr <- unique(lotr[, .(pid = as.character(pid), lotr_score = as.numeric(lotr_score))])
  s_l <- sd(lotr$lotr_score, na.rm = TRUE)
  if (!is.finite(s_l) || s_l <= 0) stop("LOT-R score has no variation.", call. = FALSE)
  lotr[, Zopt := (lotr_score - mean(lotr_score, na.rm = TRUE)) / s_l]

  # Outcome specification: fixed indexing used by the Stan model
  outcomes <- list(
    list(k = 1L, name = "b", key = "rq1", param = "mu_b_i"),
    list(k = 2L, name = "a", key = "rq2", param = "mu_a_i"),
    list(k = 3L, name = "c", key = "rq3", param = "mu_c_i"),
    list(k = 4L, name = "h", key = "rq4", param = "mu_h_i")
  )

  z_within_draw <- function(mat) {
    m <- rowMeans(mat)
    s <- apply(mat, 1L, sd)
    if (any(!is.finite(s) | s <= 0)) stop("A posterior draw has no variation across participants.", call. = FALSE)
    (mat - m) / s
  }

  diag <- list()

  for (tr in tr_vec) {

    # ---- Response time over all trials, standardised within treatment ----
    d_all <- analysis_trials(tr, "rq1-full")
    rt <- d_all[is.finite(screen_ms) & screen_ms > 0,
                .(rt_log_mean = mean(log(screen_ms))), by = .(pid = as.character(pid))]
    s_r <- sd(rt$rt_log_mean)
    if (!is.finite(s_r) || s_r <= 0) stop("Response time has no variation (", tr, ").", call. = FALSE)
    rt[, Zrt := (rt_log_mean - mean(rt_log_mean)) / s_r]
    pred <- merge(rt[, .(pid, Zrt)], lotr[, .(pid, Zopt)], by = "pid")

    # ---- Outcome blocks ----
    blocks <- list()
    for (ob in outcomes) {
      fit_k <- readRDS(model_path(ob$key, paste0("fit-", tr, "-full.rds")))
      pids  <- as.character(readRDS(model_path(ob$key, paste0("pid-levels-", tr, "-full.rds"))))
      y     <- rstan::extract(fit_k, pars = ob$param)[[1]]
      stopifnot(ncol(y) == length(pids))

      cov  <- pred[match(pids, pid)]
      keep <- is.finite(cov$Zopt) & is.finite(cov$Zrt)
      if (sum(keep) < min_Nk)
        stop("EX2: only ", sum(keep), " eligible participants for outcome '", ob$name,
             "' (", tr, "); min_Nk = ", min_Nk, ".", call. = FALSE)

      K_all <- nrow(y)
      Trep  <- min(trep0, K_all)
      set.seed(seed + ob$k)
      draw_idx <- sort(sample.int(K_all, Trep))
      y_rep <- z_within_draw(y[draw_idx, keep, drop = FALSE])

      blocks[[ob$name]] <- list(k = ob$k, pid = pids[keep], y = y_rep,
                                Zopt = cov$Zopt[keep], Zrt = cov$Zrt[keep],
                                draw_idx = draw_idx, n_pid_fit = length(pids))
    }

    Trep_joint <- min(vapply(blocks, function(b) nrow(b$y), integer(1L)))
    kid <- unlist(lapply(blocks, function(b) rep(b$k, length(b$pid))))
    data_list <- list(
      K = 4L, Nobs = length(kid), Trep = Trep_joint,
      y_rep = do.call(cbind, lapply(blocks, function(b) b$y[seq_len(Trep_joint), , drop = FALSE])),
      kid = as.integer(kid),
      Zopt = unlist(lapply(blocks, `[[`, "Zopt")),
      Zrt  = unlist(lapply(blocks, `[[`, "Zrt")),
      prior_alpha_sd    = as.numeric(pass$prior$alpha_sd),
      prior_sigma_sd    = as.numeric(pass$prior$sigma_sd),
      prior_beta_bar_sd = as.numeric(pass$prior$beta_bar_sd),
      prior_tau_sd      = as.numeric(pass$prior$tau_sd)
    )
    msg("EX2: ", tr, " | Nobs=", data_list$Nobs, " | Trep=", Trep_joint)

    fit <- rstan::sampling(
      sm, data = data_list, iter = st$iter, warmup = st$warmup, chains = st$chains,
      seed = seed, refresh = 0,
      control = list(adapt_delta = st$adapt_delta, max_treedepth = st$treedepth)
    )
    saveRDS(fit, file.path(run, paste0("fit-", tr, ".rds")))
    saveRDS(list(blocks = lapply(blocks, function(b) b[c("k", "pid", "draw_idx")]),
                 Nobs = data_list$Nobs, Trep = Trep_joint),
            file.path(run, paste0("data-", tr, ".rds")))

    tag <- paste(tr, "full", sep = "_")
    d1 <- stan_diagnostics(fit, tag = tag, treatment = tr, sample = "full", cfg = cfg,
                           max_treedepth = st$treedepth)
    n_blk <- vapply(blocks, function(b) length(b$pid), integer(1L))
    n_fit <- vapply(blocks, function(b) b$n_pid_fit, integer(1L))
    d2 <- data.table(tag = tag, treatment = tr, sample = "full",
                     check = paste0("participants_outcome_", names(blocks)),
                     value = n_blk, threshold = n_fit)
    d2[, pass := value == threshold]
    diag[[tr]] <- rbind(d1, d2)
  }

  write_diagnostics(run, rbindlist(diag))
  settings_run <- c(st, list(ex2_trep = trep0, min_Nk = min_Nk))
  if (pass$name != "main")
    settings_run <- modifyList(settings_run, list(prior = pass$prior, prior_pass = pass$name))
  model_end(run, cfg, settings = settings_run,
            key = if (pass$name == "main") "ex2" else NULL)
}
