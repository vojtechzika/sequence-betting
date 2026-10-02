# ============================================================
# 61_ex2_stan.R
#
# PURPOSE
#   Fits the EX2 hierarchical Stan model: participant-level
#   associations of dispositional optimism (LOT-R) and log
#   mean response time with behavioral outcomes across
#   RQ1--RQ4, with partial pooling of outcome-specific slopes.
#   Posterior uncertainty in participant-level outcomes is
#   propagated via a multiple-imputation (MI) likelihood over
#   posterior draws from the upstream RQ1--RQ4 models.
#
# INPUT
#   path_out/lotr_scored.csv           (pid, lotr_score)
#   path_src/master_sequences.csv      (pid, treat, stake, screen_ms)
#   path_mod/rq{1..4}_fit_sequences_<tr>_full.rds
#   path_mod/rq{1..4}_pid_levels_<tr>_full.rds
#   stan/ex2_associations.stan
#
# OUTPUT
#   path_mod/ex2_fit_<tr>.rds
#   path_out/ex2_coeffs_<tr>.csv
#
# CONFIGURATION KEYS USED
#   cfg$run$data_folder            dataset label (logging only)
#   cfg$run$seed                   global RNG seed
#   cfg$run$treatment              character vector of treatments
#   cfg$run$overwrite_models       skip logic
#   cfg$model$stan$ex2             Stan sampler settings
#   cfg$model$simulation$ex2_trep  MI draw count cap
#   cfg$design$ex2$min_Nk          min eligible participants per block
#   cfg$design$rq2$min_bets        min betting trials for mu_a eligibility
#
# NOTES
#   - LOT-R scores come from path_out/lotr_scored.csv, produced by
#     score_lotr() in the indices stage.  screen_ms comes from
#     master_sequences.csv (path_src).
#   - LOT-R z-scoring uses the full cross-treatment mean/SD from
#     lotr_scored.csv so Zopt is on a common scale across treatments.
#   - Outcome draws are z-scored within each posterior draw (across
#     participants) before entering the regression.
#   - K = 4 outcome indices are fixed: b=1, a=2, c=3, h=4.
#     The indexing is stable regardless of per-block filter attrition
#     so Stan parameters map consistently to outcome labels.
#   - Directories are created by 00_setup.R; no dir.create here.
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ex2_stan <- function(cfg) {
  
  # ----------------------------------------------------------
  # 0. Validate configuration
  # ----------------------------------------------------------
  stopifnot(
    is.list(cfg),
    !is.null(cfg$run),
    !is.null(cfg$model),
    !is.null(cfg$design),
    nzchar(as.character(cfg$run$data_folder)),
    is.numeric(cfg$run$seed),
    length(cfg$run$seed) == 1L
  )
  
  ds     <- as.character(cfg$run$data_folder)
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  # ----------------------------------------------------------
  # 1. Stan sampler settings
  # ----------------------------------------------------------
  st <- cfg$model$stan$ex2
  stopifnot(
    !is.null(st$iter),    !is.null(st$warmup),
    !is.null(st$chains),  !is.null(st$adapt_delta),
    !is.null(st$treedepth)
  )
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  # ----------------------------------------------------------
  # 2. Design thresholds
  # ----------------------------------------------------------
  ex2_trep     <- as.integer(cfg$model$simulation$ex2_trep)
  min_Nk       <- as.integer(cfg$design$ex2$min_Nk)
  rq2_min_bets <- as.integer(cfg$design$rq2$min_bets)
  
  stopifnot(
    length(ex2_trep) == 1L, is.finite(ex2_trep), ex2_trep >= 10L,
    length(min_Nk)   == 1L, is.finite(min_Nk),   min_Nk   >= 1L
  )
  
  # ----------------------------------------------------------
  # 3. Compile Stan model
  # ----------------------------------------------------------
  stan_file <- here::here("stan", "ex2.stan")
  if (!file.exists(stan_file))
    stop("Stan file not found: ", stan_file)
  sm <- rstan::stan_model(stan_file)
  
  # ----------------------------------------------------------
  # 4. Load covariates
  # ----------------------------------------------------------
  
  # LOT-R: produced by score_lotr() in the indices stage
  f_lotr <- file.path(path_out, "lotr_scored.csv")
  if (!file.exists(f_lotr))
    stop("lotr_scored.csv not found at: ", f_lotr,
         "\n  Run score_lotr(cfg) before ex2_stan().")
  lotr_dt <- fread(f_lotr, encoding = "UTF-8")
  lotr_dt[, pid        := as.character(pid)]
  lotr_dt[, lotr_score := as.numeric(lotr_score)]
  lotr_dt <- unique(lotr_dt[, .(pid, lotr_score)])
  
  # Z-score LOT-R once across the full sample in lotr_scored.csv
  lotr_m <- mean(lotr_dt$lotr_score, na.rm = TRUE)
  lotr_s <- sd(lotr_dt$lotr_score,   na.rm = TRUE)
  if (!is.finite(lotr_s) || lotr_s <= 0)
    stop("LOT-R score has zero or non-finite SD in lotr_scored.csv.")
  lotr_dt[, lotr_z := (lotr_score - lotr_m) / lotr_s]
  
  # Master sequences: response time and stake filter
  f_master <- file.path(path_src, "master_sequences.csv")
  if (!file.exists(f_master))
    stop("master_sequences.csv not found at: ", f_master)
  master <- fread(f_master, encoding = "UTF-8")
  
  req  <- c("pid", "treat", "stake", "screen_ms")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0L)
    stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  
  master[, pid       := as.character(pid)]
  master[, treat     := as.character(treat)]
  master[, stake     := as.numeric(stake)]
  master[, screen_ms := as.numeric(screen_ms)]
  master[is.na(stake), stake := 0]
  
  # ----------------------------------------------------------
  # 5. Helpers
  # ----------------------------------------------------------
  
  z_strict <- function(x, label = "x") {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x,   na.rm = TRUE)
    if (!is.finite(s) || s <= 0)
      stop("z_strict: SD is zero or non-finite for '", label, "'.")
    (x - m) / s
  }
  
  # Z-score each posterior draw (row) independently across participants (cols)
  z_within_draw <- function(mat) {
    for (t in seq_len(nrow(mat))) {
      m <- mean(mat[t, ])
      s <- sd(mat[t, ])
      if (!is.finite(s) || s <= 0)
        stop("z_within_draw: draw ", t, " has zero or non-finite SD.")
      mat[t, ] <- (mat[t, ] - m) / s
    }
    mat
  }
  
  # Build covariate table aligned to pids within one treatment subset.
  # Merges lotr_z (from lotr_dt) with per-pid mean log RT (from dtr).
  build_predictors <- function(pids, dtr) {
    d_rt    <- dtr[stake > 0 & is.finite(screen_ms) & screen_ms > 0]
    rt_pid  <- d_rt[, .(rt_log_mean = mean(log(screen_ms))), by = pid]
    pid_cov <- merge(
      data.table(pid = pids),
      lotr_dt[, .(pid, lotr_z)],
      by = "pid", all.x = TRUE
    )
    pid_cov <- merge(pid_cov, rt_pid, by = "pid", all.x = TRUE)
    setkey(pid_cov, pid)
    pid_cov[.(pids)]   # restore original order
  }
  
  # ----------------------------------------------------------
  # 6. Outcome specification  (K = 4, fixed indexing)
  # ----------------------------------------------------------
  outcome_specs <- list(
    list(k = 1L, name = "b", rq = "rq1", param = "mu_b_i"),
    list(k = 2L, name = "a", rq = "rq2", param = "mu_a_i"),
    list(k = 3L, name = "c", rq = "rq3", param = "mu_c_i"),
    list(k = 4L, name = "h", rq = "rq4", param = "mu_h_i")
  )
  
  # ----------------------------------------------------------
  # 7. Main loop over treatments
  # ----------------------------------------------------------
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit_out  <- file.path(path_mod, paste0("ex2_fit_",    tr, ".rds"))
    f_coef_out <- file.path(path_out, paste0("ex2_coeffs_", tr, ".csv"))
    
    if (should_skip(
      paths = c(f_fit_out, f_coef_out),
      cfg   = cfg,
      type  = "model",
      label = paste0("EX2 Stan fit (", ds, "/", tr, ")")
    )) next
    
    msg("EX2 Stan: assembling data for treatment '", tr, "'")
    
    # ---- Load upstream RQ posterior draws ----
    rq_data <- lapply(outcome_specs, function(ob) {
      f_fit <- file.path(path_mod,
                         paste0(ob$rq, "_fit_sequences_", tr, "_full.rds"))
      f_pid <- file.path(path_mod,
                         paste0(ob$rq, "_pid_levels_",    tr, "_full.rds"))
      if (!file.exists(f_fit)) stop("Missing upstream fit: ",        f_fit)
      if (!file.exists(f_pid)) stop("Missing upstream pid levels: ", f_pid)
      list(
        k    = ob$k,
        name = ob$name,
        pid  = as.character(readRDS(f_pid)),
        y    = rstan::extract(readRDS(f_fit))[[ob$param]]
      )
    })
    
    dtr   <- master[treat == tr]
    bet_n <- dtr[stake > 0, .N, by = pid]
    setkey(bet_n, pid)
    
    # ---- Build per-outcome eligible blocks ----
    blocks <- list()
    
    for (ob in rq_data) {
      
      pids   <- ob$pid
      y_full <- ob$y   # draws x participants matrix
      
      pid_cov <- build_predictors(pids, dtr)
      
      # RQ2: suppress participants below the minimum-bets threshold
      if (ob$name == "a") {
        n_bets <- bet_n[.(pids), N]
        n_bets[is.na(n_bets)] <- 0L
        ok_rq2 <- n_bets >= rq2_min_bets
      } else {
        ok_rq2 <- rep(TRUE, length(pids))
      }
      
      ok_pred <- is.finite(pid_cov$lotr_z) & is.finite(pid_cov$rt_log_mean)
      keep    <- ok_rq2 & ok_pred
      
      if (sum(keep) < min_Nk) {
        msg("  Outcome '", ob$name, "': only ", sum(keep),
            " eligible participants (min_Nk = ", min_Nk, "). Skipping block.")
        next
      }
      
      # Subsample posterior draws for the MI likelihood
      K_all    <- nrow(y_full)
      Trep     <- min(ex2_trep, K_all)
      set.seed(seed + ob$k)
      draw_idx <- sort(sample.int(K_all, Trep))
      
      y_rep <- y_full[draw_idx, keep, drop = FALSE]
      y_rep <- z_within_draw(y_rep)
      
      pid_cov_k <- pid_cov[keep]
      
      blocks[[ob$name]] <- list(
        k    = ob$k,
        pid  = pids[keep],
        y    = y_rep,
        Zopt = z_strict(pid_cov_k$lotr_z,     "lotr_z"),
        Zrt  = z_strict(pid_cov_k$rt_log_mean, "rt_log_mean")
      )
    }
    
    if (length(blocks) == 0L) {
      msg("  No eligible blocks for treatment '", tr, "'. Skipping.")
      next
    }
    
    # ---- Harmonize Trep across retained blocks ----
    Trep_joint <- min(vapply(blocks, function(b) nrow(b$y), integer(1L)))
    
    # ---- Stack blocks into the Stan data list ----
    # K = 4 is fixed; outcome indices must match the Stan model regardless
    # of which blocks survive the eligibility filter.
    kid     <- integer(0)
    Zopt    <- numeric(0)
    Zrt     <- numeric(0)
    y_stack <- NULL
    
    for (b in blocks) {
      Nk      <- length(b$pid)
      kid     <- c(kid,  rep(b$k, Nk))
      Zopt    <- c(Zopt, b$Zopt)
      Zrt     <- c(Zrt,  b$Zrt)
      yb      <- b$y[seq_len(Trep_joint), , drop = FALSE]
      y_stack <- if (is.null(y_stack)) yb else cbind(y_stack, yb)
    }
    
    data_list <- list(
      K     = 4L,
      Nobs  = length(kid),
      Trep  = Trep_joint,
      y_rep = y_stack,
      kid   = as.integer(kid),
      Zopt  = Zopt,
      Zrt   = Zrt
    )
    
    msg("EX2 Stan: sampling | tr=", tr,
        " | Nobs=", data_list$Nobs,
        " | Trep=", data_list$Trep,
        " | blocks=", paste(names(blocks), collapse = ","))
    
    fit <- rstan::sampling(
      sm,
      data    = data_list,
      iter    = iter_val,
      warmup  = warmup_val,
      chains  = chains_val,
      seed    = seed,
      control = list(
        adapt_delta   = adapt_delta_val,
        max_treedepth = treedepth_val
      )
    )
    
    saveRDS(fit, f_fit_out)
    msg("Saved: ", f_fit_out)
    
    # ---- Pooled coefficient summary ----
    post <- rstan::extract(fit)
    
    summ <- function(x) {
      q <- quantile(x, c(0.025, 0.975), names = FALSE)
      c(median = median(x), mean = mean(x),
        q025 = q[1L], q975 = q[2L], p_gt0 = mean(x > 0))
    }
    
    pooled_mat <- rbind(
      beta_opt_bar = summ(post$beta_opt_bar),
      beta_rt_bar  = summ(post$beta_rt_bar)
    )
    
    tbl <- data.table(
      dataset   = ds,
      treatment = tr,
      term      = rownames(pooled_mat),
      median    = pooled_mat[, "median"],
      mean      = pooled_mat[, "mean"],
      q025      = pooled_mat[, "q025"],
      q975      = pooled_mat[, "q975"],
      p_gt0     = pooled_mat[, "p_gt0"]
    )
    
    fwrite(tbl, f_coef_out)
    msg("Saved: ", f_coef_out)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}