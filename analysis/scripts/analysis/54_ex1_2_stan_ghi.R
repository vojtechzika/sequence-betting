# ============================================================
# 54_ex1_2_stan_ghi.R
#
# EX1.2  Predictors of the GHI
#
# Regresses participant-level GHI draws (chi_i) on:
#   Zopt  standardized LOT-R (optimism)
#   Zrt   standardized log mean RT on betting trials
#   Zr    standardized posterior mean risk parameter r
#
# Uncertainty in chi_i is propagated via a multiple-imputation
# likelihood: for each participant the log-likelihood is
# log(1/Trep * sum_t Normal(chi_i^(t) | mu_i, sigma)).
#
# INPUTS
#   path_mod/ex1_1_<tr>_participants.rds   (from 53_)
#   path_mod/mpl_r_draws_<tr>.rds
#   path_src/master_sequences.csv
#   stan/ex1_ghi.stan
#
# OUTPUTS
#   path_mod/ex1_2_fit_<tr>.rds
#   path_out/ex1_2_coeffs_<tr>.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ex1_2_stan_ghi <- function(cfg) {

  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))

  seed <- as.integer(cfg$run$seed)

  # Stan settings
  stan_cfg <- cfg$model$stan$ex1_2
  stopifnot(!is.null(stan_cfg))

  # Number of chi draws to use for MI (capped at available K)
  Trep_cfg <- as.integer(cfg$model$simulation$ex1_trep)
  stopifnot(length(Trep_cfg) == 1L, is.finite(Trep_cfg), Trep_cfg >= 10L)

  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))

  # Compile Stan model once
  stan_file <- here::here("stan", "ex1.stan")
  stopifnot(file.exists(stan_file))
  sm <- rstan::stan_model(stan_file)

  # Load master (covariates)
  f_master <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_master))

  master <- fread(f_master, encoding = "UTF-8")
  req <- c("pid", "treat", "stake", "screen_ms")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0L)
    stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))

  master[, pid       := as.character(pid)]
  master[, treat     := as.character(treat)]
  master[, stake     := as.numeric(stake)]
  master[, screen_ms := as.numeric(screen_ms)]
  master[is.na(stake), stake := 0]

  # LOT-R lives in its own scored output file
  f_lotr <- file.path(path_out, "lotr_scored.csv")
  if (!file.exists(f_lotr))
    stop("EX1.2: missing LOT-R file: ", f_lotr)
  lotr_dt <- fread(f_lotr, encoding = "UTF-8")
  lotr_dt[, pid        := as.character(pid)]
  lotr_dt[, lotr_score := as.numeric(lotr_score)]
  lotr_dt <- lotr_dt[, .(pid, lotr_score)]

  z_score <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x,   na.rm = TRUE)
    if (!is.finite(s) || s <= 0) stop("Cannot z-score: sd <= 0")
    (x - m) / s
  }

  outputs <- list()

  for (tr in tr_vec) {

    f_fit <- file.path(path_mod, paste0("ex1_2_fit_",    tr, ".rds"))
    f_csv <- file.path(path_out, paste0("ex1_2_coeffs_", tr, ".csv"))

    skip_model  <- should_skip(f_fit, cfg, "model",  paste0("EX1.2 fit (", tr, ")"))
    skip_output <- should_skip(f_csv, cfg, "output", paste0("EX1.2 coeffs (", tr, ")"))

    if (skip_model && skip_output) next

    # ---- Load chi_i draws ----
    f_chi <- file.path(path_mod, paste0("ex1_1_", tr, "_participants.rds"))
    if (!file.exists(f_chi))
      stop("EX1.2: missing EX1.1 participants file: ", f_chi)

    chi_obj    <- readRDS(f_chi)
    pid_levels <- as.character(chi_obj$pid_levels)
    chi_draws  <- chi_obj$chi_draws   # K x N matrix

    stopifnot(is.matrix(chi_draws), ncol(chi_draws) == length(pid_levels))

    K_all <- nrow(chi_draws)
    Trep  <- min(Trep_cfg, K_all)

    set.seed(seed)
    t_idx <- sort(sample.int(K_all, Trep, replace = FALSE))
    y_rep <- chi_draws[t_idx, , drop = FALSE]   # Trep x N

    # ---- Load MPL r draws ----
    f_r <- file.path(path_mod, paste0("mpl_r_draws_", tr, ".rds"))
    if (!file.exists(f_r))
      stop("EX1.2: missing MPL r file: ", f_r)

    r_obj <- readRDS(f_r)
    r_dt  <- data.table(
      pid    = as.character(r_obj$pid),
      r_mean = apply(r_obj$r_draws, 2, mean)
    )

    # ---- Build covariate table aligned to pid_levels ----
    d_tr <- master[treat == tr & pid %in% pid_levels]

    # Log RT: mean of log(screen_ms) on betting trials
    rt_pid <- d_tr[stake > 0 & is.finite(screen_ms) & screen_ms > 0,
                   .(rt_log_mean = mean(log(screen_ms))), by = pid]

    pid_cov <- Reduce(function(a, b) merge(a, b, by = "pid", all.x = TRUE),
                      list(data.table(pid = pid_levels), lotr_dt, r_dt, rt_pid))
    setkey(pid_cov, pid)
    pid_cov <- pid_cov[.(pid_levels)]   # restore original order

    # Keep only participants with complete covariates
    keep    <- is.finite(pid_cov$lotr_score) &
               is.finite(pid_cov$r_mean)     &
               is.finite(pid_cov$rt_log_mean)
    pid_cov <- pid_cov[keep]
    y_rep   <- y_rep[, keep, drop = FALSE]
    N       <- nrow(pid_cov)

    if (N < 10L) {
      warning("EX1.2: fewer than 10 complete cases for tr='", tr, "'. Skipping.")
      next
    }

    # ---- Standardize predictors ----
    Zopt <- z_score(pid_cov$lotr_score)
    Zrt  <- z_score(pid_cov$rt_log_mean)
    Zr   <- z_score(pid_cov$r_mean)

    # ---- Fit Stan model ----
    data_list <- list(
      N    = N,
      Trep = Trep,
      y_rep = y_rep,
      Zopt  = as.vector(Zopt),
      Zrt   = as.vector(Zrt),
      Zr    = as.vector(Zr)
    )

    if (!skip_model) {
      fit <- rstan::sampling(
        sm,
        data    = data_list,
        iter    = stan_cfg$iter,
        warmup  = stan_cfg$warmup,
        chains  = stan_cfg$chains,
        seed    = seed,
        control = list(adapt_delta   = stan_cfg$adapt_delta,
                       max_treedepth = stan_cfg$treedepth)
      )
      saveRDS(fit, f_fit)
      msg("Saved: ", f_fit)
    } else {
      fit <- readRDS(f_fit)
    }

    # ---- Summarize ----
    post <- rstan::extract(fit)

    summ <- function(x, nm) {
      q <- quantile(x, c(0.025, 0.975), names = FALSE)
      data.table(
        term   = nm,
        median = median(x),
        mean   = mean(x),
        q025   = q[1],
        q975   = q[2],
        p_gt0  = mean(x > 0)
      )
    }

    tbl <- rbindlist(list(
      summ(post$alpha, "alpha"),
      summ(post$b_opt, "beta_opt"),
      summ(post$b_rt,  "beta_rt"),
      summ(post$b_r,   "beta_r")
    ))

    tbl[, `:=`(
      dataset   = cfg$run$data_folder,
      treatment = tr,
      N         = N,
      Trep      = Trep
    )]

    setcolorder(tbl, c("dataset", "treatment", "N", "Trep",
                       "term", "median", "mean", "q025", "q975", "p_gt0"))

    if (!skip_output) {
      fwrite(tbl, f_csv)
      msg("Saved: ", f_csv)
    }

    outputs[[tr]] <- list(fit_file = f_fit, coeffs_csv = f_csv,
                          pid_levels = pid_cov$pid)
  }

  invisible(outputs)
}