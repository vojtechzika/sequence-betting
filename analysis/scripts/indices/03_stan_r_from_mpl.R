# ============================================================
# scripts/indices/01_stan_r_from_mpl.R
#
# Stan HL -> CRRA r_i (+ lambda_i), pooled + by treatment
# Also appends r_mean, r_median, inconsistent into master_sequences.csv
# ============================================================

source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_r_from_mpl <- function(cfg) {
  
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
  
  mpl    <- fread(f_mpl, encoding = "UTF-8")
  master <- fread(f_master, encoding = "UTF-8")
  
  stopifnot("pid" %in% names(mpl), "pid" %in% names(master), "treat" %in% names(master))
  
  mpl[, pid := as.character(pid)]
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  
  # Attach treat to mpl via master (master built earlier without mpl columns)
  pid_treat <- unique(master[, .(pid, treat)])
  mpl <- merge(mpl, pid_treat, by = "pid", all.x = TRUE)
  
  if (anyNA(mpl$treat)) {
    bad <- mpl[is.na(treat), unique(pid)]
    stop("Some MPL participants have no treatment match in master_sequences.csv. Example pid(s): ",
         paste(head(bad, 5), collapse = ", "))
  }
  
  # ----------------------------
  # Design parameters
  # ----------------------------
  K_cfg  <- as.integer(design$mpl$K)
  A_high <- as.numeric(design$mpl$A_high)
  A_low  <- as.numeric(design$mpl$A_low)
  B_high <- as.numeric(design$mpl$B_high)
  B_low  <- as.numeric(design$mpl$B_low)
  
  stopifnot(K_cfg > 0L, is.finite(A_high), is.finite(A_low), is.finite(B_high), is.finite(B_low))
  
  p_sched <- (1:K_cfg) / K_cfg
  
  stopifnot(
    length(p_sched) == K_cfg,
    all(is.finite(p_sched)),
    all(p_sched > 0),
    all(p_sched <= 1)
  )
  
  # ----------------------------
  # Detect HL choice columns
  # ----------------------------
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  stopifnot(length(choice_cols) == K_cfg)
  
  # ----------------------------
  # Stan sampling settings
  # ----------------------------
  stopifnot(!is.null(model$stan$mpl[[ds]]))
  st <- model$stan$mpl[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  # ----------------------------
  # Compile Stan once
  # ----------------------------
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  sm <- rstan::stan_model(stan_file)
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ============================================================
  # Helper: inconsistency flag from raw HL choices
  # - y=1 for A (safe), y=0 for B (risky)
  # - consistent pattern: non-increasing in y over rows; at most one switch
  # ============================================================
  hl_inconsistent_flag <- function(long_dt) {
    # long_dt has pid, row (1..K), y (0/1)
    setorder(long_dt, pid, row)
    long_dt[, {
      yy <- y
      # count switches between consecutive rows
      n_switch <- sum(yy[-1] != yy[-length(yy)])
      # reversal exists if ever goes 0 -> 1 later (risky then safe)
      reversal <- any(diff(yy) > 0)
      .(inconsistent = as.integer(n_switch > 1 || reversal))
    }, by = pid]
  }
  
  # ============================================================
  # Generate pooled and per-treatment outputs
  # ============================================================
  fit_one <- function(mpl_subset, tag) {
    
    fit_file  <- file.path(mod_dir, if (nzchar(tag)) paste0("mpl_fit_", tag, ".rds") else "mpl_fit.rds")
    draw_file <- file.path(mod_dir, if (nzchar(tag)) paste0("mpl_r_draws_", tag, ".rds") else "mpl_r_draws.rds")
    scored_csv <- file.path(path_out_ds(ds), if (nzchar(tag)) paste0("mpl_scored_", tag, ".csv") else "mpl_scored.csv")
    
    skip_model <- should_skip(
      paths = c(fit_file, draw_file),
      cfg   = cfg,
      type  = "model",
      label = paste0("MPL Stan (", ds, if (nzchar(tag)) paste0("/", tag) else "", ")")
    )
    
    skip_scored <- should_skip(
      paths = scored_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("MPL scored (", ds, if (nzchar(tag)) paste0("/", tag) else "", ")")
    )
    
    if (skip_model && skip_scored) return(invisible(NULL))
    
    # Long format
    long <- melt(
      mpl_subset,
      id.vars = c("pid", "treat"),
      measure.vars = choice_cols,
      variable.name = "row",
      value.name = "choice"
    )
    
    long[, row := as.integer(sub("^c", "", row))]
    long[, choice := toupper(trimws(as.character(choice)))]
    long[, y := fifelse(choice == "A", 1L, 0L)]
    
    pid_levels <- sort(unique(long$pid))
    long[, uid := match(pid, pid_levels)]
    Tobs <- nrow(long)
    
    # inconsistency (deterministic, from raw choices)
    inc_tbl <- hl_inconsistent_flag(long[, .(pid, row, y)])
    setkey(inc_tbl, pid)
    
    # Fit or load
    fit <- NULL
    post <- NULL
    
    if (!skip_model) {
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
      
      msg("MPL Stan: fitting ", if (nzchar(tag)) tag else "pooled",
          " | N=", data_list$N, " | T=", data_list$T)
      
      fit <- rstan::sampling(
        sm,
        data = data_list,
        iter = iter_val,
        warmup = warmup_val,
        chains = chains_val,
        seed = seed,
        control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
      )
      
      saveRDS(fit, fit_file)
      msg("Saved: ", fit_file)
      
      post <- rstan::extract(fit)
      saveRDS(list(pid = pid_levels, r_draws = post$r), draw_file)
      msg("Saved: ", draw_file)
      
    } else {
      # model skipped but outputs may not be: load fit
      stopifnot(file.exists(fit_file))
      fit <- readRDS(fit_file)
      post <- rstan::extract(fit)
      stopifnot(!is.null(post$r))
    }
    
    # Write scored csv if needed
    if (!skip_scored) {
      scored <- data.table(
        pid = pid_levels,
        treat = if (nzchar(tag)) tag else "pooled",
        r_mean   = apply(post$r, 2, mean),
        r_median = apply(post$r, 2, median)
      )
      scored <- merge(scored, inc_tbl, by = "pid", all.x = TRUE)
      scored[is.na(inconsistent), inconsistent := 0L]
      
      fwrite(scored, scored_csv)
      msg("Saved: ", scored_csv)
    }
    
    invisible(list(fit_file = fit_file, draw_file = draw_file, scored_csv = scored_csv))
  }
  
  # ============================================================
  # 1) POOLED
  # ============================================================
  pooled_res <- fit_one(mpl, "")
  
  # ============================================================
  # 2) PER TREATMENT
  # ============================================================
  for (tr in names(design$seq$treatments)) {
    mpl_tr <- mpl[treat == tr]
    if (nrow(mpl_tr) == 0) {
      warning("MPL Stan: no MPL rows for treat='", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    fit_one(mpl_tr, tr)
  }
  
  # ============================================================
  # Append r_mean, r_median, inconsistent to master_sequences.csv (from pooled scored)
  # ============================================================
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  f_mpl    <- file.path(path_clean_ds(ds), "mpl_scored.csv")
  
  stopifnot(file.exists(f_master), file.exists(f_mpl))
  
  master <- fread(f_master)
  mpl_scored <- fread(f_mpl)
  
  master[, pid := as.character(pid)]
  mpl_scored[, pid := as.character(pid)]
  
  # Drop previous versions
  master[, c("inconsistent", "r_inconsistent", "mpl_inconsistent") := NULL]
  
  master <- merge(
    master,
    mpl_scored[, .(pid, mpl_inconsistent = as.integer(inconsistent))],
    by = "pid",
    all.x = TRUE
  )
  
  fwrite(master, f_master)
  msg("Updated master_sequences with mpl_inconsistent")
  
  invisible(list(pooled = pooled_res))
}