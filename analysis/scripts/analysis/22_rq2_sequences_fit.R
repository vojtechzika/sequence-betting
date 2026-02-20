# scripts/analysis/22_rq2_sequences_fit.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq2_fit_sequences <- function(dataset = "pilot",
                              treat_fn = "m25",
                              seed = 12345) {
  
  f_dat <- file.path(path_mod_ds(dataset), paste0("rq2_data_", treat_fn, ".rds"))
  stopifnot(file.exists(f_dat))
  obj <- readRDS(f_dat)
  
  d <- obj$d
  pid_levels <- obj$pid_levels
  seq_levels <- obj$seq_levels
  
  stan_file <- here::here("stan", "rq2_stakes_z.stan")
  stopifnot(file.exists(stan_file))
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  data_list <- list(
    N = length(pid_levels),
    S = length(seq_levels),
    T = nrow(d),
    pid = d$pid_i,
    sid = d$sid_s,
    z = d$z
  )
  
  if (dataset == "pilot") {
    iter_val <- 4000; warmup_val <- 2000; chains_val <- 4
    adapt_delta_val <- 0.99; treedepth_val <- 15
  } else if (dataset == "main") {
    iter_val <- 4000; warmup_val <- 2000; chains_val <- 4
    adapt_delta_val <- 0.99; treedepth_val <- 15
  } else stop("dataset must be 'pilot' or 'main'")
  
  sm <- rstan::stan_model(stan_file)
  fit <- rstan::sampling(
    sm,
    data = data_list,
    iter = iter_val,
    warmup = warmup_val,
    chains = chains_val,
    seed = seed,
    control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
  )
  
  f_fit <- file.path(path_mod_ds(dataset), paste0("rq2_fit_", treat_fn, ".rds"))
  saveRDS(fit, f_fit)
  
  saveRDS(pid_levels, file.path(path_mod_ds(dataset), paste0("rq2_pid_levels_", treat_fn, ".rds")))
  saveRDS(seq_levels, file.path(path_mod_ds(dataset), paste0("rq2_seq_levels_", treat_fn, ".rds")))
  
  msg("Saved RQ2 fit: ", f_fit)
  invisible(list(fit = fit, pid_levels = pid_levels, seq_levels = seq_levels))
}

# Example:
# rq2_fit_sequences("pilot")
# rq2_fit_sequences("main")