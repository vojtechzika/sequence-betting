# ============================================================
# scripts/analysis/41_rq4_stan.R
#   RQ4 (side choice on betting trials): Stan ONLY (per treatment)
#
# Per dataset x treatment:
# 1) Loads master_sequences.csv
# 2) Keeps betting trials only: stake > 0
# 3) Uses design$seq$side_labels: heads / tails / nobet
# 4) Constructs y = 1(side == heads), 0(side == tails)
# 5) Fits stan/rq4_side.stan (per treatment; no pooling)
# 6) Saves ONLY Stan artifacts:
#      models/rq4_fit_sequences_<tr>.rds
#      models/rq4_pid_levels_<tr>.rds
#      models/rq4_seq_levels_<tr>.rds
# ============================================================


library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq4_stan <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Side labels from design
  # ----------------------------
  stopifnot(!is.null(design$seq), !is.null(design$seq$side_labels))
  stopifnot(is.list(design$seq$side_labels))
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  lab_nobet <- as.character(design$seq$side_labels$nobet)
  
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads))
  stopifnot(length(lab_tails) == 1L, nzchar(lab_tails))
  stopifnot(length(lab_nobet) == 1L, nzchar(lab_nobet))
  stopifnot(length(unique(c(lab_heads, lab_tails, lab_nobet))) == 3L)
  
  msg("RQ4 side labels:",
      " heads='", lab_heads,
      "' tails='", lab_tails,
      "' nobet='", lab_nobet, "'")
  
  # ----------------------------
  # Paths
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- here::here("stan", "rq4_side.stan")
  stopifnot(file.exists(stan_file))
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Stan sampling settings
  # ----------------------------
  st <- model$stan$rq4[[ds]]
  stopifnot(!is.null(st))
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  msg("RQ4 Stan settings (dataset=", ds, "):",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " seed=", seed)
  
  # ----------------------------
  # Load master + schema
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "seq", "stake", "side")
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
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  # Compile once
  sm <- rstan::stan_model(stan_file)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, "_full.rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, "_full.rds"))
    
    if (should_skip(
      paths = c(f_fit, f_pid, f_seq),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ4 Stan (", ds, "/", tr, ")")
    )) next
    
    d <- dt[treat == tr]
    if (nrow(d) == 0) next
    
    # betting trials only
    d <- d[is.finite(stake) & stake > 0]
    if (nrow(d) == 0) next
    
    # side must be heads/tails after stake>0
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0) {
      stop(
        "RQ4: Found side values outside {heads, tails} after stake>0 (tr='", tr, "').\n",
        "Expected: {'", lab_heads, "','", lab_tails, "'}\n",
        "Examples: ", paste(unique(head(bad$side, 10)), collapse = ", ")
      )
    }
    
    # y = 1(Heads), 0(Tails)
    d[, y := as.integer(side == lab_heads)]
    
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
    
    data_list <- list(
      N   = length(pid_levels),
      S   = length(seq_levels),
      T   = nrow(d),
      pid = as.integer(d$pid_i),
      sid = as.integer(d$sid_s),
      h   = as.integer(d$y)
    )
    
    msg("RQ4 Stan: fitting tr=", tr,
        " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T)
    
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
      stop("RQ4: Stan produced no samples for tr='", tr, "'. Fit will NOT be saved.")
    }
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels, f_pid)
    saveRDS(seq_levels, f_seq)
    
    msg("Saved: ", f_fit)
    msg("Saved: ", f_pid)
    msg("Saved: ", f_seq)
    
    outputs[[tr]] <- list(
      fit_file = f_fit,
      pid_levels_file = f_pid,
      seq_levels_file = f_seq,
      N = data_list$N,
      S = data_list$S,
      T = data_list$T
    )
  }
  
  invisible(outputs)
}