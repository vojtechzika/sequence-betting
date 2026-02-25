# ============================================================
# scripts/analysis/54_ex1_2_stan_ghi.R
#
# EX1.2 Predictors of the GHI (per treatment run)
# - Uses chi_i draws from EX1.1 participants RDS:
#     data/clean/<ds>/models/ex1_1_participants_<tr>.rds
# - Computes covariates from master_sequences.csv:
#     lotr_z    (used as raw score input; re-standardized WITHIN analysis subsample)
#     screen_ms -> log RT (computed here on betting trials; then standardized WITHIN subsample)
#     r_median  (point estimate; standardized WITHIN subsample)
#
# Uncertainty propagation:
# - Integrate over chi-draws via MI likelihood in Stan.
#
# Outputs:
# - data/clean/<ds>/models/ex1_2_fit_<tr>.rds
# - data/clean/<ds>/output/ex1_2_coeffs_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ex1_2_stan_ghi <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ------------------------------------------------------------
  # Run controls
  # ------------------------------------------------------------
  # Trep: number of chi-draws to use (subsample for speed; full propagation if = all)
  Trep_cfg <- cfg$run$ex1_trep
  if (is.null(Trep_cfg)) Trep_cfg <- 300L
  Trep_cfg <- as.integer(Trep_cfg)
  stopifnot(length(Trep_cfg) == 1L, Trep_cfg >= 10L)
  
  # RT scope is FIXED to betting trials only (stake > 0), per your decision.
  rt_scope <- "bet"
  
  # ------------------------------------------------------------
  # Files
  # ------------------------------------------------------------
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  stan_file <- here::here("stan", "ex1_ghi.stan")
  stopifnot(file.exists(stan_file))
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # Compile once
  sm <- rstan::stan_model(stan_file)
  
  # Load master once
  master <- fread(f_master, encoding = "UTF-8")
  
  req <- c("pid", "treat", "stake", "screen_ms", "lotr_z", "r_median")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0) stop("master_sequences.csv missing: ", paste(miss, collapse = ", "))
  
  master[, pid       := as.character(pid)]
  master[, treat     := as.character(treat)]
  master[, stake     := as.numeric(stake)]
  master[, screen_ms := as.numeric(screen_ms)]
  master[, lotr_z    := as.numeric(lotr_z)]
  master[, r_median  := as.numeric(r_median)]
  master[is.na(stake), stake := 0]
  
  # Helper: z-score within analysis subsample (strict; errors if sd=0)
  z_within <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) stop("Cannot z-score within subsample: sd is non-finite or 0.")
    (x - m) / s
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ----------------------------
    # Load chi draws from EX1.1 participants RDS
    # ----------------------------
    f_chi <- file.path(mod_dir, paste0("ex1_1_participants_", tr, ".rds"))
    if (!file.exists(f_chi)) {
      stop(
        "EX1.2: Missing EX1.1 participants RDS for treatment='", tr, "': ", f_chi,
        "\nRun EX1.1 participant script and ensure it saves chi draws."
      )
    }
    
    chi_obj <- readRDS(f_chi)
    
    pid_levels <- NULL
    if (!is.null(chi_obj$pid_levels)) pid_levels <- as.character(chi_obj$pid_levels)
    if (is.null(pid_levels) && !is.null(chi_obj$pid)) pid_levels <- as.character(chi_obj$pid)
    
    chi_draws <- NULL
    if (!is.null(chi_obj$chi_draws)) chi_draws <- chi_obj$chi_draws
    if (is.null(chi_draws) && !is.null(chi_obj$chi_rep)) chi_draws <- chi_obj$chi_rep
    
    if (is.null(pid_levels) || is.null(chi_draws)) {
      stop("EX1.2: ex1_1_participants_<tr>.rds must contain pid_levels (or pid) and chi_draws (or chi_rep).")
    }
    
    stopifnot(is.matrix(chi_draws))
    stopifnot(ncol(chi_draws) == length(pid_levels))
    K_all <- nrow(chi_draws)
    stopifnot(K_all >= 10L)
    
    # Subsample draws for MI (stable given seed)
    Trep <- min(Trep_cfg, K_all)
    set.seed(seed)
    t_idx <- sort(sample.int(K_all, Trep, replace = FALSE))
    y_rep <- chi_draws[t_idx, , drop = FALSE]  # Trep x N_all (aligned to pid_levels)
    
    # ----------------------------
    # Build participant covariates from master_sequences (BETTING TRIALS ONLY)
    # ----------------------------
    d <- master[treat == tr & pid %in% pid_levels]
    if (nrow(d) == 0) stop("EX1.2: No rows in master_sequences for treatment='", tr, "' after pid filter.")
    
    # Restrict to betting trials only for RT computation (stake > 0)
    d_rt <- d[is.finite(stake) & stake > 0]
    d_rt <- d_rt[is.finite(screen_ms) & screen_ms > 0]
    
    # Aggregate RT to pid
    rt_pid <- d_rt[, .(
      rt_log_mean = mean(log(screen_ms)),
      n_rt = .N
    ), by = pid]
    
    # One row per pid for lotr_z and r_median (pid-level, repeated in master)
    pid_cov <- d[, .(
      lotr_z   = lotr_z[1],
      r_median = r_median[1]
    ), by = pid]
    
    pid_cov <- merge(pid_cov, rt_pid, by = "pid", all.x = TRUE)
    
    # Align to pid_levels (same ordering as chi_draws columns)
    pid_cov <- pid_cov[pid %in% pid_levels]
    setkey(pid_cov, pid)
    pid_cov <- pid_cov[.(pid_levels)]
    
    if (anyNA(pid_cov$pid)) stop("EX1.2: Failed to align pid covariates to pid_levels (unexpected NA pid).")
    
    # Drop pids with missing predictors (required for regression)
    keep <- is.finite(pid_cov$lotr_z) & is.finite(pid_cov$r_median) & is.finite(pid_cov$rt_log_mean)
    if (!all(keep)) {
      dropped <- pid_cov[!keep, pid]
      warning(
        "EX1.2: Dropping ", length(dropped), " participant(s) with missing covariates. Example: ",
        paste(head(dropped, 10), collapse = ", ")
      )
    }
    
    pid_keep <- pid_cov[keep, pid]
    pid_cov  <- pid_cov[keep]
    y_rep    <- y_rep[, keep, drop = FALSE]
    stopifnot(ncol(y_rep) == nrow(pid_cov))
    
    N <- nrow(pid_cov)
    if (N < 5) warning("EX1.2: Very small N=", N, " for treatment='", tr, "'. Inference will be noisy.")
    
    # ----------------------------
    # Standardize ALL predictors WITHIN THIS ANALYSIS SUBSAMPLE (as prereg)
    # ----------------------------
    Zopt <- z_within(pid_cov$lotr_z)
    Zrt  <- z_within(pid_cov$rt_log_mean)
    Zr   <- z_within(pid_cov$r_median)
    
    # ----------------------------
    # Stan data and fit
    # ----------------------------
    data_list <- list(
      N    = N,
      Trep = Trep,
      y_rep = y_rep,
      Zopt = as.vector(Zopt),
      Zrt  = as.vector(Zrt),
      Zr   = as.vector(Zr)
    )
    
    # Sampling settings (from cfg)
    stopifnot(!is.null(cfg$model$stan$ex1_2), !is.null(cfg$model$stan$ex1_2[[ds]]))
    iter_val        <- as.integer(cfg$model$stan$ex1_2[[ds]]$iter)
    warmup_val      <- as.integer(cfg$model$stan$ex1_2[[ds]]$warmup)
    chains_val      <- as.integer(cfg$model$stan$ex1_2[[ds]]$chains)
    adapt_delta_val <- as.numeric(cfg$model$stan$ex1_2[[ds]]$adapt_delta)
    treedepth_val   <- as.integer(cfg$model$stan$ex1_2[[ds]]$treedepth)
    
    stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
    stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
    stopifnot(treedepth_val > 0L)
    
    msg("EX1.2 Stan settings (dataset=", ds, ", tr=", tr, "):",
        " iter=", iter_val,
        " warmup=", warmup_val,
        " chains=", chains_val,
        " adapt_delta=", adapt_delta_val,
        " treedepth=", treedepth_val,
        " seed=", seed,
        " | Trep=", Trep,
        " | rt_scope=", rt_scope)
    
    fit <- rstan::sampling(
      sm,
      data = data_list,
      iter = iter_val,
      warmup = warmup_val,
      chains = chains_val,
      seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    if (length(rstan::get_sampler_params(fit, inc_warmup = FALSE)) == 0) {
      stop("EX1.2: Stan produced no samples for treatment='", tr, "'. Fit will NOT be saved.")
    }
    
    # ----------------------------
    # Save fit + coefficient summary
    # ----------------------------
    f_fit <- file.path(mod_dir, paste0("ex1_2_fit_", tr, ".rds"))
    saveRDS(fit, f_fit)
    msg("Saved: ", f_fit)
    
    post <- rstan::extract(fit)
    
    summ <- function(x) {
      q <- stats::quantile(x, probs = c(0.025, 0.975), names = FALSE)
      c(
        median = stats::median(x),
        mean   = mean(x),
        q025   = q[1],
        q975   = q[2],
        p_gt0  = mean(x > 0)
      )
    }
    
    coefs <- rbind(
      alpha    = summ(post$alpha),
      beta_opt = summ(post$b_opt),
      beta_rt  = summ(post$b_rt),
      beta_r   = summ(post$b_r)
    )
    
    tbl <- data.table(
      term   = rownames(coefs),
      median = coefs[, "median"],
      mean   = coefs[, "mean"],
      q025   = coefs[, "q025"],
      q975   = coefs[, "q975"],
      p_gt0  = coefs[, "p_gt0"]
    )
    
    tbl[, `:=`(
      dataset   = ds,
      treatment = tr,
      N         = N,
      Trep      = Trep,
      rt_scope  = rt_scope
    )]
    
    setcolorder(tbl, c("dataset","treatment","N","Trep","rt_scope","term","median","mean","q025","q975","p_gt0"))
    
    f_csv <- file.path(out_dir, paste0("ex1_2_coeffs_", tr, ".csv"))
    fwrite(tbl, f_csv)
    msg("Saved: ", f_csv)
    
    outputs[[tr]] <- list(
      fit_file   = f_fit,
      coeffs_csv = f_csv,
      pid_levels = pid_keep
    )
  }
  
  invisible(outputs)
}

# Example:
# ex1_2_ghi_regression(cfg)