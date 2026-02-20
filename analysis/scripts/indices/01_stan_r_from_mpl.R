# scripts/indices/01_stan_r_from_mpl.R

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_r_from_mpl <- function(run, design, model) {
  
  stopifnot(
    is.list(run), is.list(design), is.list(model),
    "dataset" %in% names(run),
    "seed" %in% names(run),
    !is.null(design$mpl),
    !is.null(model$stan$mpl)
  )
  
  ds   <- run$dataset
  seed <- run$seed
  
  # ----------------------------
  # Design parameters (from config)
  # ----------------------------
  K_cfg   <- as.integer(design$mpl$K)
  A_high  <- as.numeric(design$mpl$A_high)
  A_low   <- as.numeric(design$mpl$A_low)
  B_high  <- as.numeric(design$mpl$B_high)
  B_low   <- as.numeric(design$mpl$B_low)
  
  stopifnot(
    is.finite(K_cfg), K_cfg > 0,
    is.finite(A_high), is.finite(A_low), is.finite(B_high), is.finite(B_low)
  )
  
  # Probability schedule (default: (1:K)/K if not provided)
  if (!is.null(design$mpl$p_schedule) && is.function(design$mpl$p_schedule)) {
    p_high_by_row <- design$mpl$p_schedule(K_cfg)
  } else {
    p_high_by_row <- (1:K_cfg) / K_cfg
  }
  stopifnot(length(p_high_by_row) == K_cfg, all(is.finite(p_high_by_row)),
            all(p_high_by_row >= 0), all(p_high_by_row <= 1))
  
  # ----------------------------
  # Files
  # ----------------------------
  stan_file <- here::here("stan", "holt-laurie-model.stan")
  infile    <- file.path(path_clean_ds(ds), "mpl.csv")
  
  stopifnot(file.exists(stan_file), file.exists(infile))
  
  mpl <- fread(infile, encoding = "UTF-8")
  
  out_dir <- path_out_ds(ds)
  mod_dir <- path_mod_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  fit_file  <- file.path(mod_dir, "mpl_fit.rds")
  draw_file <- file.path(mod_dir, "mpl_r_draws.rds")
  
  if (file.exists(fit_file) || file.exists(draw_file)) {
    warning(
      "MPL Stan results already exist for dataset='", ds, "'.\n",
      "Stan was NOT executed.\n",
      "To rerun, delete:\n",
      "  ", fit_file, "\n",
      "  ", draw_file
    )
    return(invisible(NULL))
  }
  
  # ----------------------------
  # Detect choice columns and enforce K
  # ----------------------------
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  K <- length(choice_cols)
  
  stopifnot(K == K_cfg)
  stopifnot(all(c("pid", choice_cols) %in% names(mpl)))
  
  # ----------------------------
  # Wide -> long
  # ----------------------------
  long <- melt(
    mpl,
    id.vars = "pid",
    measure.vars = choice_cols,
    variable.name = "row",
    value.name = "choice"
  )
  long[, row := as.integer(sub("^c", "", row))]
  long[, choice := toupper(trimws(as.character(choice)))]
  
  bad_choice <- long[!choice %in% c("A", "B")]
  if (nrow(bad_choice) > 0) {
    stop(
      "Found non A/B choices in HL data. First rows:\n",
      paste(capture.output(print(head(bad_choice, 10))), collapse = "\n")
    )
  }
  
  long[, y := fifelse(choice == "A", 1L, 0L)]
  setorder(long, pid, row)
  
  counts <- long[, .N, by = pid]
  stopifnot(all(counts$N == K))
  
  # ----------------------------
  # Diagnostics
  # ----------------------------
  diag <- long[, .(
    nA = sum(choice == "A"),
    nB = sum(choice == "B"),
    switch_row = {
      bb <- row[choice == "B"]
      if (length(bb) == 0) NA_integer_ else min(bb)
    },
    inconsistent = {
      z <- fifelse(choice == "A", 0L, 1L)  # A=0, B=1
      any(diff(z) == -1L)                  # A after B
    }
  ), by = pid]
  
  # ----------------------------
  # Stan data
  # ----------------------------
  pid_levels <- sort(unique(long$pid))
  long[, uid := match(pid, pid_levels)]
  
  Tobs <- nrow(long)
  
  data_list <- list(
    N   = length(pid_levels),
    T   = Tobs,
    uid = long$uid,
    p   = p_high_by_row[long$row],
    A1  = rep(A_high, Tobs),
    A2  = rep(A_low,  Tobs),
    B1  = rep(B_high, Tobs),
    B2  = rep(B_low,  Tobs),
    y   = long$y
  )
  
  # Normalize line endings (hash stability)
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # ----------------------------
  # Sampling settings (from config)
  # ----------------------------
  if (!ds %in% names(model$stan$mpl)) {
    stop("model$stan$mpl has no entry for dataset='", ds, "'.")
  }
  st <- model$stan$mpl[[ds]]
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(
    iter_val > 0, warmup_val >= 0, chains_val > 0,
    is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1,
    treedepth_val > 0
  )
  
  msg("MPL Stan settings:",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed)
  
  sm <- rstan::stan_model(stan_file)
  
  fit <- rstan::sampling(
    sm,
    data = data_list,
    iter = iter_val,
    warmup = warmup_val,
    chains = chains_val,
    seed = seed,
    control = list(
      adapt_delta = adapt_delta_val,
      max_treedepth = treedepth_val
    )
  )
  
  # ----------------------------
  # Outputs
  # ----------------------------
  logfile <- file.path(out_dir, "mpl_run_log.txt")
  cat(
    paste0(
      "dataset=", ds, "\n",
      "timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "K=", K, "\n",
      "iter=", iter_val, " warmup=", warmup_val, " chains=", chains_val, "\n",
      "adapt_delta=", adapt_delta_val, " treedepth=", treedepth_val, "\n",
      "seed=", seed, "\n",
      "A_high=", A_high, " A_low=", A_low, " B_high=", B_high, " B_low=", B_low, "\n"
    ),
    file = logfile
  )
  
  saveRDS(fit, file.path(mod_dir, "mpl_fit.rds"))
  
  post <- rstan::extract(fit)
  r_draws <- post$r  # iterations x N
  
  saveRDS(
    list(pid = pid_levels, r_draws = r_draws),
    file.path(mod_dir, "mpl_r_draws.rds")
  )
  
  scored <- data.table(
    pid = pid_levels,
    r_mean   = apply(post$r, 2, mean),
    r_median = apply(post$r, 2, median),
    r_q025   = apply(post$r, 2, quantile, 0.025),
    r_q975   = apply(post$r, 2, quantile, 0.975),
    lambda_mean   = apply(post$lambda, 2, mean),
    lambda_median = apply(post$lambda, 2, median)
  )
  scored <- merge(scored, diag, by = "pid", all.x = TRUE)
  
  outfile <- file.path(path_clean_ds(ds), "mpl_scored.csv")
  fwrite(scored, outfile)
  msg("MPL scored table saved:", outfile, "| rows:", nrow(scored))
  
  invisible(list(
    fit = fit,
    scored = scored,
    r_draws = r_draws,
    pid = pid_levels,
    K = K,
    hl_params = list(
      K = K,
      A_high = A_high, A_low = A_low,
      B_high = B_high, B_low = B_low
    ),
    out_dir = out_dir,
    mod_dir = mod_dir
  ))
}