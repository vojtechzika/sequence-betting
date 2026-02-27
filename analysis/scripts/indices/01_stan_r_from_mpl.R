# ============================================================
# scripts/indices/01_stan_r_from_mpl.R
#
# Estimates individual CRRA risk preferences (r_i) and choice
# precision (lambda_i) from the Holt–Laury task using a
# hierarchical Bayesian model in Stan.
#
# The model is estimated:
#   (1) pooled across all participants, and
#   (2) separately by treatment (e.g., m25, m19),
# using identical likelihoods and priors.
#
# Outputs (models/):
#   - mpl_fit.rds                 (pooled)
#   - mpl_r_draws.rds             (pooled posterior draws)
#   - mpl_fit_<tr>.rds            (by treatment)
#   - mpl_r_draws_<tr>.rds        (by treatment)
#
# Outputs (data/clean/<ds>/):
#   - mpl_scored.csv              (pooled summary)
#   - mpl_scored_<tr>.csv         (treatment-specific summary)
#
# Treatment labels are taken from design$seq$treatments.
# No hardcoding of multipliers or group names.
# ============================================================

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_r_from_mpl <- function(cfg) {
  
  stopifnot(
    is.list(cfg),
    !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model),
    !is.null(cfg$run$dataset),
    !is.null(cfg$run$seed),
    !is.null(cfg$design$mpl),
    !is.null(cfg$model$stan$mpl)
  )
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  
  design <- cfg$design
  model  <- cfg$model
  
  # ----------------------------
  # Files
  # ----------------------------
  stan_file <- here::here("stan", "holt-laurie-model.stan")
  f_mpl     <- file.path(path_clean_ds(ds), "mpl.csv")
  f_master  <- file.path(path_clean_ds(ds), "master_sequences.csv")
  
  stopifnot(file.exists(stan_file), file.exists(f_mpl), file.exists(f_master))
  
  mpl     <- fread(f_mpl)
  master  <- fread(f_master)
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  
  pid_treat <- unique(master[, .(pid, treat)])
  mpl[, pid := as.character(pid)]
  mpl <- merge(mpl, pid_treat, by = "pid", all.x = TRUE)
  
  if (anyNA(mpl$treat)) {
    stop("Some MPL participants have no treatment match.")
  }
  
  # ----------------------------
  # Design parameters
  # ----------------------------
  K_cfg  <- design$mpl$K
  A_high <- design$mpl$A_high
  A_low  <- design$mpl$A_low
  B_high <- design$mpl$B_high
  B_low  <- design$mpl$B_low
  
  if (!is.null(design$mpl$p_schedule)) {
    p_sched <- design$mpl$p_schedule(K_cfg)
  } else {
    p_sched <- (1:K_cfg) / K_cfg
  }
  
  # ----------------------------
  # Detect HL choice columns
  # ----------------------------
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c","",choice_cols)))]
  stopifnot(length(choice_cols) == K_cfg)
  
  # ----------------------------
  # Compile Stan once
  # ----------------------------
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  sm <- rstan::stan_model(stan_file)
  
  st <- model$stan$mpl[[ds]]
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ============================================================
  # Internal fitter
  # ============================================================
  
  fit_one <- function(mpl_subset, tag) {
    
    fit_file  <- file.path(mod_dir,
                           if (nzchar(tag)) paste0("mpl_fit_", tag, ".rds")
                           else "mpl_fit.rds")
    
    draw_file <- file.path(mod_dir,
                           if (nzchar(tag)) paste0("mpl_r_draws_", tag, ".rds")
                           else "mpl_r_draws.rds")
    
    scored_csv <- file.path(path_clean_ds(ds),
                            if (nzchar(tag)) paste0("mpl_scored_", tag, ".csv")
                            else "mpl_scored.csv")
    
    if (should_skip(c(fit_file, draw_file), cfg, "model",
                    paste0("MPL Stan (", ds,
                           if (nzchar(tag)) paste0("/", tag), ")"))) return()
    
    if (should_skip(scored_csv, cfg, "output",
                    paste0("MPL scored (", ds,
                           if (nzchar(tag)) paste0("/", tag), ")"))) return()
    
    long <- melt(
      mpl_subset,
      id.vars = c("pid","treat"),
      measure.vars = choice_cols,
      variable.name = "row",
      value.name = "choice"
    )
    
    long[, row := as.integer(sub("^c","",row))]
    long[, choice := toupper(trimws(choice))]
    long[, y := fifelse(choice=="A",1L,0L)]
    
    pid_levels <- sort(unique(long$pid))
    long[, uid := match(pid, pid_levels)]
    
    Tobs <- nrow(long)
    
    data_list <- list(
      N   = length(pid_levels),
      T   = Tobs,
      uid = long$uid,
      p   = p_sched[long$row],
      A1  = rep(A_high, Tobs),
      A2  = rep(A_low,  Tobs),
      B1  = rep(B_high, Tobs),
      B2  = rep(B_low,  Tobs),
      y   = long$y
    )
    
    fit <- rstan::sampling(
      sm,
      data = data_list,
      iter = st$iter,
      warmup = st$warmup,
      chains = st$chains,
      seed = seed,
      control = list(
        adapt_delta = st$adapt_delta,
        max_treedepth = st$treedepth
      )
    )
    
    saveRDS(fit, fit_file)
    
    post <- rstan::extract(fit)
    
    saveRDS(
      list(pid = pid_levels, r_draws = post$r),
      draw_file
    )
    
    scored <- data.table(
      pid = pid_levels,
      treat = if (nzchar(tag)) tag else "pooled",
      r_mean   = apply(post$r,2,mean),
      r_median = apply(post$r,2,median),
      r_q025   = apply(post$r,2,quantile,0.025),
      r_q975   = apply(post$r,2,quantile,0.975),
      lambda_mean   = apply(post$lambda,2,mean),
      lambda_median = apply(post$lambda,2,median)
    )
    
    fwrite(scored, scored_csv)
  }
  
  # ============================================================
  # 1) POOLED
  # ============================================================
  
  fit_one(mpl, "")
  
  # ============================================================
  # 2) PER TREATMENT
  # ============================================================
  
  for (tr in names(design$seq$treatments)) {
    mpl_tr <- mpl[treat == tr]
    fit_one(mpl_tr, tr)
  }
}