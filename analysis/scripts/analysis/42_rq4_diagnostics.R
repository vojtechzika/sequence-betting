# ============================================================
# 42_rq4_diagnostics.R
#   RQ4 diagnostics: convergence checks + PPC overdispersion test
#
# PURPOSE
#   1) Convergence: Rhat, ESS, divergences for each treatment fit
#   2) PPC overdispersion check: compares observed variance of
#      sequence-wise Heads rates to posterior predictive replicates
#      (y_rep drawn in Stan generated quantities, not re-simulated here)
#   3) If P(var_rep >= var_obs) < (1 - ppc_overdisp_cut): flags that
#      the Beta-Binomial robustness fit is required per preregistration
#   4) Saves a PPC density figure per treatment
#   5) Saves a single diagnostics CSV: output/rq4_diagnostics.csv
#
# STRUCTURE
#   Single "full" tag per treatment -- no confirmatory subset for RQ4
#   (see 41_rq4_stan.R for rationale).
#
# CONFIG KEYS USED
#   cfg$run$treatment
#   cfg$run$seed
#   cfg$run$overwrite_outputs
#   cfg$design$rq4$ppc_overdisp_cut  -- trigger threshold (default 0.95)
#   cfg$design$drift$rq4$outcome     -- "side" (for drift key lookup)
#   cfg$model$ppc$rq4_k              -- number of y_rep draws to use (scalar)
#
# INPUTS
#   path_src/master_sequences.csv
#   path_mod/rq4_fit_sequences_<tr>_full.rds   (contains y_rep, mu_h)
#   path_mod/rq4_pid_levels_<tr>_full.rds
#   path_mod/rq4_seq_levels_<tr>_full.rds
#   path_mod/drift_decisions.rds
#
# OUTPUTS
#   path_out/rq4_diagnostics.csv
#     columns: treatment, drift, N, S, T, K_all, K_used,
#              rhat_max, ess_min, divergences,
#              var_obs, var_rep_median, var_rep_q025, var_rep_q975,
#              ratio_var, p_overdisp, cut, bb_required
#   path_fig/rq4_ppc_<tr>.png
#
# CALL ORDER
#   rq4_stan(cfg)        -- 41_rq4_stan.R
#   rq4_diagnostics(cfg) -- this script
#   rq4_tables(cfg)      -- 43_rq4_tables.R
# ============================================================

library(data.table)
library(rstan)
library(ggplot2)

rq4_diagnostics <- function(cfg) {
  
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Config
  # ----------------------------
  
  # Overdispersion trigger cutoff (preregistered): Beta-Binomial robustness
  # required when P(var_rep >= var_obs) < (1 - cut), i.e. the Binomial model
  # underpredicts observed dispersion with posterior probability >= cut.
  cut <- as.numeric(design$rq4$ppc_overdisp_cut %||% 0.95)
  stopifnot(is.finite(cut), cut > 0, cut < 1)
  
  # Number of y_rep posterior draws to use for PPC (scalar from 02_models.R)
  K_ppc <- as.integer(model$ppc$rq4_k)
  stopifnot(length(K_ppc) == 1L, is.finite(K_ppc), K_ppc >= 50L)
  
  # Side labels
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(nzchar(lab_heads), nzchar(lab_tails))
  
  # Drift outcome key
  outcome_key <- as.character(design$drift$rq4$outcome)   # "side"
  stopifnot(nzchar(outcome_key))
  
  # ----------------------------
  # Paths
  # ----------------------------
  infile  <- file.path(path_src, "master_sequences.csv")
  f_drift <- file.path(path_mod, "drift_decisions.rds")
  f_out   <- file.path(path_out, "rq4_diagnostics.csv")
  
  stopifnot(file.exists(infile))
  stopifnot(file.exists(f_drift))
  
  if (should_skip(f_out, cfg, "output", "RQ4 diagnostics")) {
    return(invisible(NULL))
  }
  
  dir.create(path_fig, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Load shared data
  # ----------------------------
  drift_decisions <- readRDS(f_drift)
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "stake", "side") %in% names(dt)))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  # ----------------------------
  # Loop over treatments
  # ----------------------------
  all_rows <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(path_mod, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_pid <- file.path(path_mod, paste0("rq4_pid_levels_",    tr, "_full.rds"))
    f_seq <- file.path(path_mod, paste0("rq4_seq_levels_",    tr, "_full.rds"))
    
    if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) {
      warning("RQ4 diagnostics: missing artifacts for tr='", tr, "'. Skipping.")
      next
    }
    
    fit        <- readRDS(f_fit)
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    S          <- length(seq_levels)
    stopifnot(length(pid_levels) >= 1L, S >= 2L)
    
    # ---- Drift type (for logging) ----
    drift_key  <- paste(tr, "full", outcome_key, sep = "_")
    drift_res  <- drift_decisions[[drift_key]]
    drift_type <- if (!is.null(drift_res) && isTRUE(drift_res$drift)) {
      as.character(drift_res$drift_type)
    } else {
      "none"
    }
    
    # ---- Convergence diagnostics ----
    summ     <- summary(fit)$summary
    rhat_max <- max(summ[, "Rhat"],  na.rm = TRUE)
    ess_min  <- min(summ[, "n_eff"], na.rm = TRUE)
    sp       <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divs     <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    
    msg("RQ4 diagnostics (", tr, "):",
        " Rhat_max=",     round(rhat_max, 4),
        " ESS_min=",      round(ess_min,  0),
        " divergences=",  divs,
        " drift=",        drift_type)
    
    if (rhat_max > 1.01) {
      warning("RQ4 diagnostics: Rhat_max=", round(rhat_max, 4),
              " > 1.01 for tr='", tr, "'. Results may be unreliable.")
    }
    if (divs > 0L) {
      warning("RQ4 diagnostics: ", divs, " divergent transitions for tr='",
              tr, "'.")
    }
    
    # ---- Observed sequence-wise Heads rates ----
    # Restrict to the exact betting trials used in the fit (same pid/seq sets)
    d <- dt[treat == tr & pid %in% pid_levels & seq %in% seq_levels]
    d <- d[is.finite(stake) & stake > 0]
    
    if (nrow(d) == 0L) {
      warning("RQ4 diagnostics: no betting trials found for tr='", tr, "'. Skipping.")
      next
    }
    
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0L) {
      stop("RQ4 diagnostics: invalid side values after stake > 0 for tr='", tr, "'.")
    }
    
    d[, h     := as.integer(side == lab_heads)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$sid_s))
    
    # Sequence-wise observed rates (sequences with no betting trials get NA -> excluded)
    obs <- d[, .(n_s = .N, h_s = sum(h)), by = sid_s]
    obs <- obs[n_s > 0L]
    
    if (nrow(obs) < 3L) {
      warning("RQ4 diagnostics: fewer than 3 sequences with betting trials",
              " for tr='", tr, "'. Skipping PPC.")
      next
    }
    
    p_obs   <- obs$h_s / obs$n_s
    var_obs <- stats::var(p_obs)
    keep_s  <- obs$sid_s                  # integer indices into seq_levels
    
    # ---- PPC overdispersion using y_rep from Stan ----
    # y_rep[T] was drawn in generated quantities of rq4.stan.
    # Each draw k gives one replicated dataset; we compute sequence-wise
    # replicated rates from those draws and compare their variance to var_obs.
    # This correctly integrates over all sources of posterior uncertainty
    # (alpha, u_i, beta_s, drift) without re-simulation in R.
    post   <- rstan::extract(fit)
    
    if (is.null(post$y_rep)) {
      stop("RQ4 diagnostics: y_rep not found in fit for tr='", tr,
           "'. Refit with updated rq4.stan that includes y_rep.")
    }
    
    y_rep_all <- post$y_rep          # K_all x T matrix
    K_all     <- nrow(y_rep_all)
    Tn        <- ncol(y_rep_all)
    stopifnot(Tn == nrow(d))
    
    K_use <- min(K_ppc, K_all)
    set.seed(seed)
    k_idx <- sort(sample.int(K_all, K_use, replace = FALSE))
    y_rep_sub <- y_rep_all[k_idx, , drop = FALSE]
    
    # Map trial index -> sid_s for the obs-subset sequences only
    sid_vec   <- d$sid_s             # length T; integer in 1..S
    
    var_rep <- numeric(K_use)
    for (k in seq_len(K_use)) {
      h_rep_k  <- as.integer(y_rep_sub[k, ])
      # sequence-wise replicated rate (only over sequences in keep_s)
      p_rep_s <- vapply(keep_s, function(s) {
        idx <- which(sid_vec == s)
        if (length(idx) == 0L) NA_real_ else mean(h_rep_k[idx])
      }, numeric(1L))
      var_rep[k] <- stats::var(p_rep_s, na.rm = TRUE)
    }
    
    # P(var_rep >= var_obs): probability that Binomial model replicates
    # at least as much dispersion as observed.
    # Low values -> model underpredicts dispersion -> Beta-Binomial required.
    p_overdisp  <- mean(var_rep >= var_obs)
    ratio_var   <- var_obs / stats::median(var_rep)
    bb_required <- is.finite(p_overdisp) && p_overdisp < (1 - cut)
    
    msg("RQ4 PPC overdispersion (", tr, "):",
        " var_obs=",          sprintf("%.6f", var_obs),
        " median(var_rep)=",  sprintf("%.6f", stats::median(var_rep)),
        " ratio=",            sprintf("%.2f",  ratio_var),
        " P(var_rep>=var_obs)=", sprintf("%.3f", p_overdisp),
        " cut=",              sprintf("%.2f",  cut),
        " BB_required=",      bb_required)
    
    if (bb_required) {
      msg("RQ4 PPC: Overdispersion detected for tr='", tr, "'.\n",
          "  -> Beta-Binomial robustness fit required per preregistration.\n",
          "  -> Run rq4_side_bb.stan on aggregated sequence-level counts.")
    } else {
      msg("RQ4 PPC: No overdispersion detected for tr='", tr,
          "' at cut=", sprintf("%.2f", cut), ".")
    }
    
    # ---- PPC density figure ----
    f_fig <- file.path(path_fig, paste0("rq4_ppc_", tr, ".png"))
    if (!should_skip(f_fig, cfg, "output",
                     paste0("RQ4 PPC figure (", tr, ")"))) {
      
      df_dens <- data.table(var_rep = var_rep)
      
      p_fig <- ggplot(df_dens, aes(x = var_rep)) +
        geom_density(fill  = "#74ADD1", alpha = 0.5,
                     colour = "#2166AC", linewidth = 0.7) +
        geom_vline(xintercept = var_obs,
                   colour = "#C0392B", linewidth = 0.8, linetype = "dashed") +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf(
                   "P(var_rep \u2265 var_obs) = %.3f\nDrift: %s\nBB required: %s",
                   p_overdisp, drift_type, bb_required
                 ),
                 hjust = 1.1, vjust = 1.5, size = 3.2, colour = "grey30") +
        labs(
          title = paste0("RQ4 PPC overdispersion check (", tr, ")"),
          x     = "Var(sequence-wise Heads rates) — replicated",
          y     = "Density"
        ) +
        theme_classic(base_size = 11)
      
      ggsave(f_fig, p_fig, width = 6, height = 4, dpi = 300)
      msg("Saved: ", f_fig)
    }
    
    # ---- Collect row ----
    all_rows[[tr]] <- data.table(
      treatment      = tr,
      drift          = drift_type,
      N              = length(pid_levels),
      S              = S,
      T              = Tn,
      K_all          = K_all,
      K_used         = K_use,
      rhat_max       = round(rhat_max, 4),
      ess_min        = round(ess_min,  0),
      divergences    = divs,
      var_obs        = var_obs,
      var_rep_median = stats::median(var_rep),
      var_rep_q025   = stats::quantile(var_rep, 0.025),
      var_rep_q975   = stats::quantile(var_rep, 0.975),
      ratio_var      = ratio_var,
      p_overdisp     = p_overdisp,
      cut            = cut,
      bb_required    = bb_required
    )
  }
  
  # ----------------------------
  # Save combined diagnostics table
  # ----------------------------
  if (length(all_rows) > 0L) {
    fwrite(rbindlist(all_rows, fill = TRUE), f_out)
    msg("Saved: ", f_out)
  } else {
    warning("RQ4 diagnostics: no rows produced; diagnostics CSV not written.")
  }
  
  invisible(all_rows)
}