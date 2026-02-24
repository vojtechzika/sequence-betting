# ============================================================
# scripts/analysis/41_rq4_stan.R
#   RQ4 (side choice on betting trials): Stan ONLY (per treatment)
#
# What this does (per dataset x treatment):
# 1) Loads master_sequences.csv
# 2) Keeps betting trials only: stake > 0
# 3) Uses cfg$design$seq$side_labels to interpret side values:
#      heads = "H", tails = "O", nobet = "NB" (default in your design)
#    After stake>0, side must be heads or tails.
# 4) Constructs y = 1(side == heads), 0(side == tails)
# 5) Fits Stan model stan/rq4_side.stan (per treatment; NO pooling)
# 6) Saves ONLY Stan artifacts:
#      models/rq4_fit_sequences_<tr>.rds
#      models/rq4_pid_levels_<tr>.rds
#      models/rq4_seq_levels_<tr>.rds
#
# Safeguard:
# - If any artifact exists for a treatment, that treatment is skipped.
# ============================================================

library(data.table)
library(rstan)


### diagnostic helper to detect the need for beta-binominal robustness check (per prereg):
rq4_beta_bi_reminder <- function(fit, d_fit, seq_levels, seed, cut = 0.95) {
  stopifnot(!is.null(fit), is.data.table(d_fit))
  stopifnot(all(c("sid_s", "y") %in% names(d_fit)))
  stopifnot(length(seq_levels) >= 2L)
  
  post <- rstan::extract(fit)
  if (is.null(post$mu_h)) {
    msg("RQ4 PPC: fit has no mu_h; cannot run overdispersion diagnostic.")
    return(invisible(NULL))
  }
  
  mu_draws <- post$mu_h   # iters x S
  iters <- nrow(mu_draws)
  S <- ncol(mu_draws)
  stopifnot(S == length(seq_levels))
  
  # observed counts by sequence (betting trials only)
  obs <- d_fit[, .(
    n_s = .N,
    h_s = sum(y, na.rm = TRUE)
  ), by = sid_s]
  
  # ensure all sequences present in 1..S
  obs_full <- data.table(sid_s = 1:S)
  obs_full <- merge(obs_full, obs, by = "sid_s", all.x = TRUE)
  obs_full[is.na(n_s), `:=`(n_s = 0L, h_s = 0L)]
  
  # if any sequences have 0 trials, variance of p_hat becomes unstable;
  # we exclude those sequences from the dispersion check
  keep <- obs_full$n_s > 0
  if (sum(keep) < 3L) {
    msg("RQ4 PPC: <3 sequences with n_s>0; skipping overdispersion diagnostic.")
    return(invisible(NULL))
  }
  
  n_s <- as.integer(obs_full$n_s[keep])
  p_obs <- obs_full$h_s[keep] / n_s
  var_obs <- stats::var(p_obs)
  
  # simulate replicated dispersion under the fitted model
  # use mu_h draws; for each draw, Binomial replicate counts per sequence
  set.seed(seed)
  var_rep <- numeric(iters)
  
  # subset mu_draws to sequences with n_s>0
  mu_sub <- mu_draws[, keep, drop = FALSE]
  
  for (k in seq_len(iters)) {
    h_rep <- stats::rbinom(n = length(n_s), size = n_s, prob = mu_sub[k, ])
    p_rep <- h_rep / n_s
    var_rep[k] <- stats::var(p_rep)
  }
  
  p_under <- mean(var_rep < var_obs)  # model underpredicts dispersion
  ratio <- var_obs / stats::median(var_rep)
  
  msg(
    "RQ4 PPC (overdispersion check): ",
    "var_obs=", sprintf("%.6f", var_obs),
    " | median(var_rep)=", sprintf("%.6f", stats::median(var_rep)),
    " | var_obs/median(var_rep)=", sprintf("%.2f", ratio),
    " | P(var_rep < var_obs)=", sprintf("%.3f", p_under)
  )
  
  if (is.finite(p_under) && p_under >= cut) {
    msg(
      "RQ4 prereg reminder: Bernoulli-logit PPC suggests underdispersion (sequence-wise Heads rates more variable than Binomial).\n",
      "-> Run Beta–Binomial robustness model (aggregated counts by sequence) as preregistered."
    )
  }
  
  invisible(list(var_obs = var_obs, var_rep = var_rep, p_under = p_under, ratio = ratio))
}

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq4_stan <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Side labels from cfg
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
  # Files
  # ----------------------------
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  stan_file <- here::here("stan", "rq4_side.stan")
  stopifnot(file.exists(stan_file))
  
  # Normalize line endings (prevents hash-mismatch churn)
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  # ----------------------------
  # Stan sampling settings
  # ----------------------------
  stopifnot(!is.null(model$stan), !is.null(model$stan$rq4), !is.null(model$stan$rq4[[ds]]))
  st <- model$stan$rq4[[ds]]
  
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
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Load + schema checks
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
  
  # -------------------------------------------------
  # Pre-check artifacts BEFORE compiling Stan
  # -------------------------------------------------
  to_run <- character(0)
  
  for (tr in tr_vec) {
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    
    if (file.exists(f_fit) || file.exists(f_pid) || file.exists(f_seq)) {
      warning(
        "RQ4 Stan artifacts already exist for dataset='", ds, "', treatment='", tr, "'.\n",
        "Stan was NOT executed.\n",
        "To rerun, delete (as applicable):\n",
        "  ", f_fit, "\n",
        "  ", f_pid, "\n",
        "  ", f_seq
      )
    } else {
      to_run <- c(to_run, tr)
    }
  }
  
  if (length(to_run) == 0L) {
    msg("RQ4 Stan: No treatments require estimation. Nothing to run.")
    return(invisible(NULL))
  }
  
  # Compile once per session (ONLY if needed)
  sm <- rstan::stan_model(stan_file)
  
  outputs <- list()
  
  for (tr in to_run) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    
    d <- dt[treat == tr]
    if (nrow(d) == 0) {
      warning("RQ4 Stan: No rows after filtering treat == '", tr, "' (dataset='", ds, "'). Skipping.")
      next
    }
    
    # ----------------------------
    # Keep betting trials only
    # ----------------------------
    d <- d[is.finite(stake) & stake > 0]
    if (nrow(d) == 0) {
      warning("RQ4 Stan: No betting trials (stake > 0) for treatment='", tr, "'. Skipping.")
      next
    }
    
    # After stake>0, side must be heads or tails
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0) {
      stop(
        "RQ4: Found side values outside {heads, tails} after stake>0 (treatment='", tr, "').\n",
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
    
    msg("RQ4 Stan: fitting treatment=", tr,
        " | N=", data_list$N, " | S=", data_list$S, " | T=", data_list$T)
    
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
    
    # --- prereg PPC reminder for Beta–Binomial robustness (if needed) ---
    # configurable cutoff (optional)
    cut_ppc <- cfg$design$rhos$rq4_ppc_overdisp_cut
    if (is.null(cut_ppc)) cut_ppc <- 0.95
    cut_ppc <- as.numeric(cut_ppc)
    
    rq4_beta_bi_reminder(
      fit = fit,
      d_fit = d[, .(sid_s, y)],         # betting trials only already
      seq_levels = seq_levels,
      seed = seed,
      cut = cut_ppc
    )
    
    # Do NOT save empty fits
    if (length(rstan::get_sampler_params(fit, inc_warmup = FALSE)) == 0) {
      stop("RQ4: Stan produced no samples for treatment='", tr, "'. Fit will NOT be saved.")
    }
    
    saveRDS(fit, f_fit)
    saveRDS(pid_levels, f_pid)
    saveRDS(seq_levels, f_seq)
    
    msg("Saved:", f_fit)
    msg("Saved:", f_pid)
    msg("Saved:", f_seq)
    
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

# Example:
# rq4_stan(cfg)