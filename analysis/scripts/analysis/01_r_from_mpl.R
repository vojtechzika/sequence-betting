# scripts/analysis/01_r_from_mpl.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

r_from_mpl <- function(dataset = "pilot") {
  
  # =================================================
  # DESIGN PARAMETERS (edit here if HL design changes)
  # =================================================
  
  HL_NUM_CHOICES <- 10L  # expected number of HL rows (choices c1..cK)
  
  # Lottery payoffs (currency units, must match oTree settings)
  A_high <- 20.0
  A_low  <- 16.0
  B_high <- 38.5
  B_low  <- 1.0
  
  # Probabilities by row: 1/K, 2/K, ..., 1 (standard HL)
  # If you ever change the probability schedule, edit this function.
  p_high_by_row_fun <- function(K) (1:K) / K
  
  # Stan model file (relative to analysis/ project root)
  stan_file <- here::here("stan", "holt-laurie-model.stan")
  
  # =================================================
  # INPUT
  # =================================================
  
  infile <- file.path(path_clean_ds(dataset), "mpl.csv")
  stopifnot(file.exists(infile))
  stopifnot(file.exists(stan_file))
  
  mpl <- data.table::fread(infile, encoding = "UTF-8")
  
  # =================================================
  # DETECT CHOICE COLUMNS AND ENFORCE DESIGN
  # =================================================
  
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  K <- length(choice_cols)
  
  stopifnot(K == HL_NUM_CHOICES)
  stopifnot(all(c("pid", choice_cols) %in% names(mpl)))
  
  # =================================================
  # WIDE -> LONG
  # =================================================
  
  long <- data.table::melt(
    mpl,
    id.vars = "pid",
    measure.vars = choice_cols,
    variable.name = "row",
    value.name = "choice"
  )
  long[, row := as.integer(sub("^c", "", row))]
  long[, choice := toupper(trimws(as.character(choice)))]
  
  # Enforced choices: require all A/B, no missing / junk values
  bad_choice <- long[!choice %in% c("A", "B")]
  if (nrow(bad_choice) > 0) {
    stop(
      "Found non A/B choices in HL data. First rows:\n",
      paste(capture.output(print(head(bad_choice, 10))), collapse = "\n")
    )
  }
  
  # y = 1 if choose A, 0 if choose B
  long[, y := data.table::fifelse(choice == "A", 1L, 0L)]
  data.table::setorder(long, pid, row)
  
  # Enforced choices: each pid must have exactly K rows
  counts <- long[, .N, by = pid]
  stopifnot(all(counts$N == K))
  
  # =================================================
  # DIAGNOSTICS (A safe, B risky)
  # =================================================
  
  diag <- long[, .(
    nA = sum(choice == "A"),
    nB = sum(choice == "B"),
    switch_row = {
      bb <- row[choice == "B"]
      if (length(bb) == 0) NA_integer_ else min(bb)
    },
    inconsistent = {
      z <- data.table::fifelse(choice == "A", 0L, 1L)  # A=0, B=1
      any(diff(z) == -1)                               # A after B
    }
  ), by = pid]
  
  # =================================================
  # STAN DATA
  # =================================================
  
  p_high_by_row <- p_high_by_row_fun(K)
  
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
  
  # Normalize line endings (prevents hash mismatch churn)
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # =================================================
  # SAMPLING SETTINGS (dataset-specific)
  # =================================================
  
  if (dataset == "pilot") {
    
    iter_val        <- 1000
    warmup_val      <- 500
    chains_val      <- 2
    adapt_delta_val <- 0.9
    treedepth_val   <- 10
    
  } else if (dataset == "main") {
    
    iter_val        <- 4000
    warmup_val      <- 2000
    chains_val      <- 4
    adapt_delta_val <- 0.995
    treedepth_val   <- 15
    
  } else {
    stop("dataset must be 'pilot' or 'main'")
  }
  
  msg("Running Stan with settings:")
  msg("iter:", iter_val,
      "| warmup:", warmup_val,
      "| chains:", chains_val,
      "| adapt_delta:", adapt_delta_val,
      "| max_treedepth:", treedepth_val)
  
  # Compile once per session
  sm <- rstan::stan_model(stan_file)
  
  fit <- rstan::sampling(
    sm,
    data = data_list,
    iter = iter_val,
    warmup = warmup_val,
    chains = chains_val,
    seed = 12345,
    control = list(
      adapt_delta = adapt_delta_val,
      max_treedepth = treedepth_val
    )
  )
  
  # =================================================
  # OUTPUTS
  # =================================================
  
  saveRDS(fit, file.path(path_clean_ds(dataset), "mpl_fit.rds"))
  
  post <- rstan::extract(fit)
  
  # Posterior draws for uncertainty propagation
  r_draws <- post$r  # iterations x N
  saveRDS(
    list(pid = pid_levels, r_draws = r_draws),
    file.path(path_clean_ds(dataset), "mpl_r_draws.rds")
  )
  
  # Convenience summary table
  scored <- data.table::data.table(
    pid = pid_levels,
    r_mean   = apply(post$r, 2, mean),
    r_median = apply(post$r, 2, median),
    r_q025   = apply(post$r, 2, quantile, 0.025),
    r_q975   = apply(post$r, 2, quantile, 0.975),
    lambda_mean   = apply(post$lambda, 2, mean),
    lambda_median = apply(post$lambda, 2, median)
  )
  scored <- merge(scored, diag, by = "pid", all.x = TRUE)
  
  outfile <- file.path(path_clean_ds(dataset), "mpl_scored.csv")
  data.table::fwrite(scored, outfile)
  
  msg("MPL scored table saved:", outfile, "| rows:", nrow(scored))
  
  invisible(list(
    fit = fit,
    scored = scored,
    r_draws = r_draws,
    pid = pid_levels,
    K = K,
    hl_params = list(
      HL_NUM_CHOICES = HL_NUM_CHOICES,
      A_high = A_high,
      A_low  = A_low,
      B_high = B_high,
      B_low  = B_low
    )
  ))
}