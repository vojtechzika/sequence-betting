# ============================================================
# scripts/analysis/61_ex2_stan.R
#
# EX2: Participant-level associations with optimism and response time
# (and addendum: risk) across RQ1--RQ4.
#
# Per dataset x treatment:
# 1) Loads master_sequences.csv
# 2) Loads RQ1--RQ4 Stan fits + pid_levels
# 3) Uses participant-level posterior outcomes DIRECTLY from Stan:
#    - RQ1: mu_b_i  (baseline betting probability)
#    - RQ2: mu_a_i  (baseline mean proportional stake deviation)
#    - RQ3: mu_c_i  (participant-level expected loss / CE loss summary used in prereg)
#    - RQ4: mu_h_i  (baseline Heads probability)
#
# 4) For each outcome k, forms an outcome-specific analysis subset (per prereg):
#    participants are included if that outcome and predictors are defined.
#    No global exclusion; stacking happens across outcomes with different N_k.
#
# 5) Standardization:
#    - Outcomes: for each posterior draw t, standardize y_{t,i} within that
#      outcome subset (mean 0, sd 1).
#    - Predictors: within the same outcome subset, z-score:
#        Zopt: lotr_z
#        Zrt:  log mean RT on betting trials (stake>0)
#        Zr :  r_median  (addendum predictor)
#
# 6) Stacks all outcome blocks and fits partial-pooled regression in Stan
#    with MI likelihood over posterior draws.
#
# Outputs:
#   models/ex2_fit_<tr>.rds
#   models/ex2_data_<tr>.rds   (exact data used to fit)
#   output/ex2_coeffs_<tr>.csv (pooled + outcome-specific slopes)
# ============================================================

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ex2_stan <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Run controls
  # ----------------------------
  # Trep: number of retained posterior draws used in MI likelihood.
  # If Trep > available draws, all available draws are used.
  ex2_trep <- cfg$run$ex2_trep
  if (is.null(ex2_trep)) ex2_trep <- 500L
  ex2_trep <- as.integer(ex2_trep)
  stopifnot(length(ex2_trep) == 1L, ex2_trep >= 10L)
  
  # Minimum per-outcome participants required to include the outcome block
  min_Nk <- cfg$run$ex2_min_Nk
  if (is.null(min_Nk)) min_Nk <- 10L
  min_Nk <- as.integer(min_Nk)
  stopifnot(length(min_Nk) == 1L, min_Nk >= 3L)
  
  # RQ2 undefined rule per prereg: fewer than 3 betting trials -> undefined
  rq2_min_bets <- cfg$run$ex2_rq2_min_bets
  if (is.null(rq2_min_bets)) rq2_min_bets <- 3L
  rq2_min_bets <- as.integer(rq2_min_bets)
  stopifnot(length(rq2_min_bets) == 1L, rq2_min_bets >= 1L)
  
  # ----------------------------
  # Paths
  # ----------------------------
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  stan_file <- here::here("stan", "ex2_associations.stan")
  stopifnot(file.exists(stan_file))
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # Stan settings
  stopifnot(!is.null(cfg$model$stan), !is.null(cfg$model$stan$ex2), !is.null(cfg$model$stan$ex2[[ds]]))
  st <- cfg$model$stan$ex2[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  msg("EX2 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed,
      " | ex2_trep=", ex2_trep,
      " | min_Nk=", min_Nk,
      " | rq2_min_bets=", rq2_min_bets)
  
  # Compile once
  sm <- rstan::stan_model(stan_file)
  
  # ----------------------------
  # Load master once
  # ----------------------------
  master <- fread(f_master, encoding = "UTF-8")
  req <- c("pid", "treat", "stake", "screen_ms", "lotr_z", "r_median")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0) stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[, screen_ms := as.numeric(screen_ms)]
  master[, lotr_z := as.numeric(lotr_z)]
  master[, r_median := as.numeric(r_median)]
  master[is.na(stake), stake := 0]
  
  # z-score within a vector (strict)
  z_strict <- function(x) {
    x <- as.numeric(x)
    m <- mean(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) stop("Cannot z-score: sd is non-finite or 0.")
    (x - m) / s
  }
  
  # standardize each draw across participants (matrix Trep x Nk)
  z_within_draw <- function(mat) {
    stopifnot(is.matrix(mat))
    for (t in seq_len(nrow(mat))) {
      m <- mean(mat[t, ])
      s <- stats::sd(mat[t, ])
      if (!is.finite(s) || s <= 0) stop("Cannot draw-standardize outcome: sd is non-finite or 0 at draw ", t)
      mat[t, ] <- (mat[t, ] - m) / s
    }
    mat
  }
  
  # Build predictors for a given pid set using betting trials only for RT
  build_predictors <- function(d_trials, pid_levels) {
    
    # RT from betting trials only (stake>0), valid screen_ms
    d_rt <- d_trials[is.finite(stake) & stake > 0 & is.finite(screen_ms) & screen_ms > 0]
    rt_pid <- d_rt[, .(
      rt_log_mean = mean(log(screen_ms)),
      n_rt = .N
    ), by = pid]
    
    # pid-level covariates (repeated in master)
    pid_cov <- d_trials[, .(
      lotr_z = lotr_z[1],
      r_median = r_median[1]
    ), by = pid]
    
    pid_cov <- merge(pid_cov, rt_pid, by = "pid", all.x = TRUE)
    
    # align to pid_levels order
    setkey(pid_cov, pid)
    pid_cov <- pid_cov[.(pid_levels)]
    
    pid_cov
  }
  
  # Helper: take posterior matrix (K_all x N_all), subsample rows to Trep
  subsample_draws <- function(mat, Trep, seed_local) {
    stopifnot(is.matrix(mat))
    K_all <- nrow(mat)
    Trep_use <- min(Trep, K_all)
    set.seed(seed_local)
    idx <- sort(sample.int(K_all, Trep_use, replace = FALSE))
    mat[idx, , drop = FALSE]
  }
  
  # Helper: summarize coefficient draws
  summ <- function(x) {
    qs <- stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
    c(
      median = stats::median(x, na.rm = TRUE),
      mean   = mean(x, na.rm = TRUE),
      q025   = unname(qs[[1]]),
      q975   = unname(qs[[2]]),
      p_gt0  = mean(x > 0, na.rm = TRUE)
    )
  }
  
  # ----------------------------
  # Run per treatment
  # ----------------------------
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ---------- load fits + pid levels ----------
    f_rq1 <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_rq2 <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_rq3 <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, ".rds"))
    f_rq4 <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    
    f_p1 <- file.path(mod_dir, paste0("rq1_pid_levels_", tr, ".rds"))
    f_p2 <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_p3 <- file.path(mod_dir, paste0("rq3_pid_levels_", tr, ".rds"))
    f_p4 <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, ".rds"))
    
    stopifnot(file.exists(f_rq1), file.exists(f_rq2), file.exists(f_rq3), file.exists(f_rq4))
    stopifnot(file.exists(f_p1), file.exists(f_p2), file.exists(f_p3), file.exists(f_p4))
    
    fit1 <- readRDS(f_rq1); pid1 <- as.character(readRDS(f_p1))
    fit2 <- readRDS(f_rq2); pid2 <- as.character(readRDS(f_p2))
    fit3 <- readRDS(f_rq3); pid3 <- as.character(readRDS(f_p3))
    fit4 <- readRDS(f_rq4); pid4 <- as.character(readRDS(f_p4))
    
    post1 <- rstan::extract(fit1)
    post2 <- rstan::extract(fit2)
    post3 <- rstan::extract(fit3)
    post4 <- rstan::extract(fit4)
    
    # ---------- participant-level outcomes DIRECTLY from Stan ----------
    if (is.null(post1$mu_b_i)) stop("EX2: RQ1 fit missing 'mu_b_i' (rerun RQ1 with updated Stan).")
    if (is.null(post2$mu_a_i)) stop("EX2: RQ2 fit missing 'mu_a_i' (rerun RQ2 with updated Stan).")
    if (is.null(post3$mu_c_i)) stop("EX2: RQ3 fit missing 'mu_c_i' (rerun RQ3 with updated Stan).")
    if (is.null(post4$mu_h_i)) stop("EX2: RQ4 fit missing 'mu_h_i' (rerun RQ4 with updated Stan).")
    
    mu_b_all <- post1$mu_b_i
    mu_a_all <- post2$mu_a_i
    mu_c_all <- post3$mu_c_i
    mu_h_all <- post4$mu_h_i
    
    stopifnot(is.matrix(mu_b_all), is.matrix(mu_a_all), is.matrix(mu_c_all), is.matrix(mu_h_all))
    stopifnot(ncol(mu_b_all) == length(pid1))
    stopifnot(ncol(mu_a_all) == length(pid2))
    stopifnot(ncol(mu_c_all) == length(pid3))
    stopifnot(ncol(mu_h_all) == length(pid4))
    
    # Trials for this treatment
    dtr <- master[treat == tr]
    if (nrow(dtr) == 0) stop("EX2: no rows in master_sequences for treatment='", tr, "'.")
    
    # Betting counts for RQ2-definedness (stake>0)
    bet_n <- dtr[stake > 0, .N, by = pid]
    setkey(bet_n, pid)
    
    # Outcome definitions list
    # kid mapping: 1=b, 2=a, 3=c, 4=h (fixed)
    outcome_blocks <- list(
      list(k = 1L, name = "b", pid_levels = pid1, y_draws = mu_b_all),
      list(k = 2L, name = "a", pid_levels = pid2, y_draws = mu_a_all),
      list(k = 3L, name = "c", pid_levels = pid3, y_draws = mu_c_all),
      list(k = 4L, name = "h", pid_levels = pid4, y_draws = mu_h_all)
    )
    
    blocks_used <- list()
    
    for (ob in outcome_blocks) {
      
      k      <- ob$k
      nm     <- ob$name
      pids   <- ob$pid_levels
      y_full <- ob$y_draws
      
      # Align predictors to pid levels (order matters)
      pid_cov <- build_predictors(dtr, pids)
      
      # Apply RQ2 undefined rule: fewer than rq2_min_bets betting trials in this treatment
      # (applies only to outcome a)
      if (nm == "a") {
        n_bets <- bet_n[.(pids), N]
        n_bets[is.na(n_bets)] <- 0L
        ok_a <- n_bets >= rq2_min_bets
      } else {
        ok_a <- rep(TRUE, length(pids))
      }
      
      ok_pred <- is.finite(pid_cov$lotr_z) & is.finite(pid_cov$r_median) & is.finite(pid_cov$rt_log_mean)
      keep <- ok_a & ok_pred
      
      if (sum(keep) < min_Nk) {
        warning("EX2: dropping outcome '", nm, "' in treatment='", tr, "' (kept ", sum(keep), " < min_Nk=", min_Nk, ").")
        next
      }
      
      # Subset and subsample draws
      y_rep <- subsample_draws(y_full[, keep, drop = FALSE], Trep = ex2_trep, seed_local = seed + 1000L + k)
      
      # Outcome standardization within draw (per prereg)
      y_rep <- z_within_draw(y_rep)
      
      # Predictor standardization within this outcome subset (per prereg)
      pid_cov_k <- pid_cov[keep]
      Zopt <- z_strict(pid_cov_k$lotr_z)
      Zrt  <- z_strict(pid_cov_k$rt_log_mean)
      Zr   <- z_strict(pid_cov_k$r_median)
      
      blocks_used[[nm]] <- list(
        k = k,
        name = nm,
        pid = pids[keep],
        Trep = nrow(y_rep),
        y_rep = y_rep,
        Zopt = Zopt,
        Zrt  = Zrt,
        Zr   = Zr
      )
    }
    
    if (length(blocks_used) == 0) {
      warning("EX2: no outcome blocks usable for treatment='", tr, "'. Skipping Stan.")
      next
    }
    
    # Need a common Trep across blocks to stack into a single matrix.
    Trep_joint <- min(vapply(blocks_used, function(b) b$Trep, integer(1)))
    if (Trep_joint < 10L) stop("EX2: too few draws after alignment (Trep_joint=", Trep_joint, ") for tr='", tr, "'.")
    
    # Stack blocks
    kid <- integer(0)
    Zopt <- numeric(0)
    Zrt  <- numeric(0)
    Zr   <- numeric(0)
    y_stack <- NULL
    
    # Fixed K=4 outcomes in Stan, even if some are absent; kid indicates which ones appear.
    for (nm in names(blocks_used)) {
      b <- blocks_used[[nm]]
      Nk <- length(b$pid)
      
      kid  <- c(kid, rep(b$k, Nk))
      Zopt <- c(Zopt, b$Zopt)
      Zrt  <- c(Zrt,  b$Zrt)
      Zr   <- c(Zr,   b$Zr)
      
      yb <- b$y_rep[seq_len(Trep_joint), , drop = FALSE]  # Trep_joint x Nk
      y_stack <- if (is.null(y_stack)) yb else cbind(y_stack, yb)
    }
    
    stopifnot(ncol(y_stack) == length(kid))
    
    data_list <- list(
      K = 4L,
      Nobs = length(kid),
      Trep = Trep_joint,
      y_rep = y_stack,
      kid = as.integer(kid),
      Zopt = as.vector(Zopt),
      Zrt  = as.vector(Zrt),
      Zr   = as.vector(Zr)
    )
    
    msg("EX2: fitting (", ds, "/", tr, "):",
        " blocks=", paste(names(blocks_used), collapse = ","),
        " | Nobs=", data_list$Nobs,
        " | Trep=", data_list$Trep)
    
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
      stop("EX2: Stan produced no samples for treatment='", tr, "'. Fit will NOT be saved.")
    }
    
    # Save fit + data used
    f_fit <- file.path(mod_dir, paste0("ex2_fit_", tr, ".rds"))
    f_dat <- file.path(mod_dir, paste0("ex2_data_", tr, ".rds"))
    saveRDS(fit, f_fit)
    saveRDS(list(
      dataset = ds,
      treatment = tr,
      Trep_joint = Trep_joint,
      blocks_used = lapply(blocks_used, function(b) list(k=b$k, name=b$name, pid=b$pid)),
      data_list = data_list
    ), f_dat)
    
    msg("Saved: ", f_fit)
    msg("Saved: ", f_dat)
    
    # Coefficient summaries to CSV
    post <- rstan::extract(fit)
    
    # pooled means
    pooled <- rbind(
      beta_opt_bar = summ(post$beta_opt_bar),
      beta_rt_bar  = summ(post$beta_rt_bar),
      beta_r_bar   = summ(post$beta_r_bar)
    )
    
    # outcome-specific slopes (K=4)
    stopifnot(!is.null(post$beta_opt_k), !is.null(post$beta_rt_k), !is.null(post$beta_r_k))
    opt_k <- post$beta_opt_k
    rt_k  <- post$beta_rt_k
    r_k   <- post$beta_r_k
    
    # Build table
    rows <- list()
    
    rows[[length(rows)+1]] <- data.table(
      term = rownames(pooled),
      outcome = "pooled",
      median = pooled[, "median"],
      mean   = pooled[, "mean"],
      q025   = pooled[, "q025"],
      q975   = pooled[, "q975"],
      p_gt0  = pooled[, "p_gt0"]
    )
    
    out_names <- c("b","a","c","h")
    for (k in 1:4) {
      rows[[length(rows)+1]] <- data.table(
        term = c("beta_opt","beta_rt","beta_r"),
        outcome = out_names[k],
        median = c(median(opt_k[,k]), median(rt_k[,k]), median(r_k[,k])),
        mean   = c(mean(opt_k[,k]),   mean(rt_k[,k]),   mean(r_k[,k])),
        q025   = c(quantile(opt_k[,k],0.025), quantile(rt_k[,k],0.025), quantile(r_k[,k],0.025)),
        q975   = c(quantile(opt_k[,k],0.975), quantile(rt_k[,k],0.975), quantile(r_k[,k],0.975)),
        p_gt0  = c(mean(opt_k[,k] > 0), mean(rt_k[,k] > 0), mean(r_k[,k] > 0))
      )
    }
    
    tbl <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    tbl[, `:=`(
      dataset = ds,
      treatment = tr,
      Nobs = data_list$Nobs,
      Trep = Trep_joint,
      blocks = paste(names(blocks_used), collapse = ",")
    )]
    
    setcolorder(tbl, c("dataset","treatment","Nobs","Trep","blocks","outcome","term","median","mean","q025","q975","p_gt0"))
    
    f_csv <- file.path(out_dir, paste0("ex2_coeffs_", tr, ".csv"))
    fwrite(tbl, f_csv)
    msg("Saved: ", f_csv)
    
    outputs[[tr]] <- list(
      fit_file = f_fit,
      data_file = f_dat,
      coeffs_csv = f_csv
    )
  }
  
  invisible(outputs)
}


# Example:
# ex2_stan_associations(cfg)