# ============================================================
# scripts/analysis/42_rq4_diagnostics.R
#   RQ4 diagnostics (PPC overdispersion check) -- NO refits
#
# Per dataset x treatment:
# - Loads rq4 fit + seq_levels
# - Reconstructs betting-trial sequence-wise Heads rates from master
# - PPC: compares observed variance of sequence-wise p_hat to Binomial reps
# - If P(var_rep < var_obs) >= cut: prints prereg reminder to run Beta–Binomial robustness
# - Saves: output/rq4_ppc_overdispersion_<tr>.csv
#
# Requires in cfg:
# - cfg$design$rq4$ppc_overdisp_cut (optional; default 0.95)
# - cfg$model$ppc$rq4_k[[ds]]       (optional; default uses all posterior draws)
# ============================================================

library(data.table)
library(rstan)

rq4_diagnostics <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # side labels
  stopifnot(!is.null(design$seq), !is.null(design$seq$side_labels))
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, length(lab_tails) == 1L)
  
  # cutoff (prereg reminder trigger)
  cut <- design$rq4$ppc_overdisp_cut
  if (is.null(cut)) cut <- 0.95
  cut <- as.numeric(cut)
  stopifnot(is.finite(cut), cut > 0, cut < 1)
  
  # PPC draw subsampling (optional)
  krep <- NULL
  if (!is.null(model$ppc) && !is.null(model$ppc$rq4_k) && !is.null(model$ppc$rq4_k[[ds]])) {
    krep <- as.integer(model$ppc$rq4_k[[ds]])
    stopifnot(length(krep) == 1L, krep >= 10L)
  }
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid","treat","seq","stake","side") %in% names(dt)))
  
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq4_seq_levels_", tr, ".rds"))
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) next
    
    fit <- readRDS(f_fit)
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    stopifnot(length(pid_levels) >= 1L, length(seq_levels) >= 2L)
    
    # reconstruct betting trials used (restricted to pid_levels/seq_levels)
    d <- dt[treat == tr & pid %in% pid_levels & seq %in% seq_levels]
    d <- d[is.finite(stake) & stake > 0]
    if (nrow(d) == 0) next
    
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0) {
      stop("RQ4 diagnostics: invalid side values after stake>0 in tr='", tr, "'.")
    }
    
    d[, y := as.integer(side == lab_heads)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$sid_s))
    
    # observed sequence-wise heads rates
    obs <- d[, .(n_s = .N, h_s = sum(y)), by = sid_s]
    
    S <- length(seq_levels)
    obs_full <- data.table(sid_s = 1:S)
    obs_full <- merge(obs_full, obs, by = "sid_s", all.x = TRUE)
    obs_full[is.na(n_s), `:=`(n_s = 0L, h_s = 0L)]
    
    keep <- obs_full$n_s > 0
    if (sum(keep) < 3L) {
      msg("RQ4 diagnostics: <3 sequences with n_s>0 for tr=", tr, " -> skipping PPC overdispersion check.")
      next
    }
    
    n_s <- as.integer(obs_full$n_s[keep])
    p_obs <- obs_full$h_s[keep] / n_s
    var_obs <- stats::var(p_obs)
    
    post <- rstan::extract(fit)
    if (is.null(post$mu_h)) stop("RQ4 fit missing generated quantity mu_h: ", f_fit)
    
    mu_draws <- post$mu_h  # iters x S
    stopifnot(is.matrix(mu_draws), ncol(mu_draws) == S)
    iters_all <- nrow(mu_draws)
    
    # subsample posterior draws if requested
    if (is.null(krep)) {
      idx <- seq_len(iters_all)
    } else {
      Krep <- min(krep, iters_all)
      set.seed(seed)
      idx <- sort(sample.int(iters_all, Krep, replace = FALSE))
    }
    
    mu_sub <- mu_draws[idx, keep, drop = FALSE]
    iters <- nrow(mu_sub)
    
    # replicate dispersion under Binomial model
    set.seed(seed)
    var_rep <- numeric(iters)
    
    for (k in seq_len(iters)) {
      h_rep <- stats::rbinom(n = length(n_s), size = n_s, prob = mu_sub[k, ])
      p_rep <- h_rep / n_s
      var_rep[k] <- stats::var(p_rep)
    }
    
    p_under <- mean(var_rep < var_obs)
    ratio <- var_obs / stats::median(var_rep)
    
    msg(
      "RQ4 PPC (overdispersion check): tr=", tr,
      " | var_obs=", sprintf("%.6f", var_obs),
      " | median(var_rep)=", sprintf("%.6f", stats::median(var_rep)),
      " | var_obs/median(var_rep)=", sprintf("%.2f", ratio),
      " | P(var_rep < var_obs)=", sprintf("%.3f", p_under),
      " | cut=", sprintf("%.2f", cut)
    )
    
    if (is.finite(p_under) && p_under >= cut) {
      msg(
        "RQ4 PPC: Posterior predictive check suggests underpredicted dispersion in sequence-wise Heads rates (relative to Binomial).\n",
        "-> Run Beta--Binomial robustness model (aggregated counts by sequence) as preregistered."
      )
    } else {
      msg(
        "RQ4 PPC: No strong indication of underpredicted dispersion at cut=", sprintf("%.2f", cut), ".\n",
        "-> Beta--Binomial robustness remains optional."
      )
    }
    
    # save diagnostics row (+ keep count)
    f_out <- file.path(out_dir, paste0("rq4_ppc_overdispersion_", tr, ".csv"))
    diag_tbl <- data.table(
      treat = tr,
      n_sequences_total = S,
      n_sequences_used = sum(keep),
      var_obs = var_obs,
      var_rep_median = stats::median(var_rep),
      ratio_var = ratio,
      p_under = p_under,
      cut = cut,
      iters_used = iters
    )
    
    if (!should_skip(
      paths = f_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ4 PPC overdispersion (", ds, "/", tr, ")")
    )) {
      fwrite(diag_tbl, f_out)
      msg("Saved: ", f_out)
    }
    
    outputs[[tr]] <- diag_tbl
  }
  
  invisible(outputs)
}