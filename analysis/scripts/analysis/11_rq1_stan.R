# ============================================================
# scripts/analysis/11_rq1_stan.R
#   - fits RQ1 model per treatment (NO pooled)
#   - writes ONLY Stan artifacts (fit + levels)
#   - has safeguard: if artifacts exist, skips that treatment
# ============================================================

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq1_stan <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  run    <- cfg$run
  model  <- cfg$model
  
  ds   <- as.character(run$dataset)
  seed <- as.integer(run$seed)
  
  overwrite_models <- isTRUE(cfg$run$overwrite_models)
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Files
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- here::here("stan", "rq1_bets.stan")
  stopifnot(file.exists(stan_file))
  
  # Normalize line endings (prevents hash-mismatch churn)
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # ----------------------------
  # Stan sampling settings (from model cfg)
  # ----------------------------
  stopifnot(!is.null(model$stan), !is.null(model$stan$rq1), !is.null(model$stan$rq1[[ds]]))
  st <- model$stan$rq1[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  msg("RQ1 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed)
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Load + schema checks
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "stake", "seq")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(
      "master_sequences.csv missing columns: ", paste(missing, collapse = ", "),
      "\nExpected at least: ", paste(required, collapse = ", ")
    )
  }
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  
  # -------------------------------------------------
  # Pre-check artifacts BEFORE compiling Stan
  # -------------------------------------------------
  to_run <- character(0)
  
  for (tr in tr_vec) {
    
    f_fit        <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq1_pid_levels_", tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq1_seq_levels_", tr, ".rds"))
    
    if (should_skip(
      paths = c(f_fit, f_pid_levels, f_seq_levels),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ1 Stan (", ds, "/", tr, ")")
    )) {
      next
    }
    
    to_run <- c(to_run, tr)
  }
  
  # Compile once per session (ONLY if needed)
  sm <- rstan::stan_model(stan_file)
  
  # ----------------------------
  # Run per treatment (NO pooled)
  # ----------------------------
  outputs <- list()
  
  for (tr in to_run) {
    
    f_fit        <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq1_pid_levels_",  tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq1_seq_levels_", tr, ".rds"))
    
    d <- dt[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ1 Stan: No rows after filtering treat == '", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    
    d[, y := as.integer(!is.na(stake) & stake > 0)]
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    
    stopifnot(length(seq_levels) == 64L)
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    
    data_list <- list(
      N   = length(pid_levels),
      S   = length(seq_levels),
      T   = nrow(d),
      pid = d$pid_i,
      sid = d$sid_s,
      y   = d$y
    )
    
    msg("RQ1 Stan: fitting treatment=", tr, " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T)
    
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
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels, f_pid_levels)
    saveRDS(seq_levels, f_seq_levels)
    
    msg("Saved:", f_fit)
    msg("Saved:", f_pid_levels)
    msg("Saved:", f_seq_levels)
    
    outputs[[tr]] <- list(
      fit_file = f_fit,
      pid_levels_file = f_pid_levels,
      seq_levels_file = f_seq_levels,
      N = data_list$N,
      S = data_list$S,
      T = data_list$T
    )
  }
  
  invisible(outputs)
}