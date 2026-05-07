# ============================================================
# scripts/analysis/72_ex3_rq2_implied_prob.R
#
# PURPOSE
#   EX3 / RQ2: Implied subjective win probabilities from
#   positive-stake betting trials. Run for every treatment
#   in cfg$run$treatment (both FN and FP).
#
# INPUT
#   path_src/master_sequences.csv
#   path_mod/mpl_r_draws_<tr>.rds  (preferred)
#   path_mod/mpl_r_draws.rds        (fallback, pooled)
#
# OUTPUT (per treatment <tr>)
#   path_out/ex3_rq2_trials_<tr>.csv
#   path_out/ex3_rq2_participants_incl_<tr>.csv
#   path_out/ex3_rq2_participants_excl_<tr>.csv
#   path_mod/ex3_rq2_<tr>.rds
#
# Additional output
#   path_out/ex3_rq2_treatment_comparison.csv
#
# NOTES
#   - Implied p is the FOC solution under CRRA expected utility.
#   - Trials where the FOC has no interior solution on [0,1]
#     are flagged; summaries reported both including and
#     excluding flagged trials (prereg §EX3 RQ2).
#   - Paths and dataset tag come from cfg / 00_setup.R globals;
#     no values are hardcoded.
# ============================================================

ex3_rq2_implied_p <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  
  # ----------------------------------------------------------
  # Resolve design constants from cfg
  # ----------------------------------------------------------
  
  endowment <- as.numeric(cfg$design$seq$endowment)
  stopifnot(length(endowment) == 1L, is.finite(endowment), endowment > 0)
  
  tr_mult <- vapply(cfg$design$seq$treatments, as.numeric, numeric(1L))
  tr_run  <- unique(as.character(cfg$run$treatment))
  stopifnot(all(tr_run %in% names(tr_mult)))
  
  # ----------------------------------------------------------
  # Load master data (ETL output lives in path_src)
  # ----------------------------------------------------------
  
  f_master <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  master <- fread(f_master, encoding = "UTF-8")
  
  req_cols <- c("pid", "treat", "seq", "stake")
  miss     <- setdiff(req_cols, names(master))
  if (length(miss)) stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  
  # ----------------------------------------------------------
  # Helper: posterior summaries for a probability vector
  # ----------------------------------------------------------
  
  summ_prob <- function(x) {
    if (!length(x) || all(!is.finite(x))) {
      return(c(median = NA_real_, mean = NA_real_,
               q025 = NA_real_, q975 = NA_real_, p_gt_05 = NA_real_))
    }
    q <- quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    c(
      median  = median(x, na.rm = TRUE),
      mean    = mean(x,   na.rm = TRUE),
      q025    = q[1L],
      q975    = q[2L],
      p_gt_05 = mean(x > 0.5, na.rm = TRUE)
    )
  }
  
  # ----------------------------------------------------------
  # Empty-output scaffolds (reused when a treatment has no
  # positive-stake trials)
  # ----------------------------------------------------------
  
  empty_trials <- function()
    data.table(
      data_folder   = character(), treatment = character(),
      pid           = character(), sequence  = character(),
      stake         = numeric(),   n_draws   = integer(),
      n_flagged     = integer(),   share_flagged = numeric(),
      median        = numeric(),   mean      = numeric(),
      q025          = numeric(),   q975      = numeric(),
      p_gt_05       = numeric()
    )
  
  empty_pid <- function()
    data.table(
      data_folder   = character(), treatment     = character(),
      pid           = character(), n_trials      = integer(),
      n_draws_total = integer(),   n_flagged     = integer(),
      share_flagged = numeric(),   median        = numeric(),
      mean          = numeric(),   q025          = numeric(),
      q975          = numeric(),   p_gt_05       = numeric()
    )
  
  # ----------------------------------------------------------
  # Main loop over treatments
  # ----------------------------------------------------------
  
  outputs    <- list()
  treat_rows <- list()
  
  for (tr in tr_run) {
    
    msg("EX3 RQ2 implied probabilities: treatment ", tr)
    
    m <- as.numeric(tr_mult[[tr]])
    stopifnot(is.finite(m), m > 1)
    
    f_trials <- file.path(path_out, paste0("ex3_rq2_trials_",             tr, ".csv"))
    f_pincl  <- file.path(path_out, paste0("ex3_rq2_participants_incl_",  tr, ".csv"))
    f_pexcl  <- file.path(path_out, paste0("ex3_rq2_participants_excl_",  tr, ".csv"))
    f_rds    <- file.path(path_mod, paste0("ex3_rq2_",                    tr, ".rds"))
    
    skip_trials <- should_skip(f_trials, cfg, "output",
                               paste0("EX3 RQ2 trials (",         tr, ")"))
    skip_pincl  <- should_skip(f_pincl,  cfg, "output",
                               paste0("EX3 RQ2 participants incl (", tr, ")"))
    skip_pexcl  <- should_skip(f_pexcl,  cfg, "output",
                               paste0("EX3 RQ2 participants excl (", tr, ")"))
    skip_rds    <- should_skip(f_rds,    cfg, "model",
                               paste0("EX3 RQ2 RDS (",            tr, ")"))
    
    # ---- subset positive-stake trials for this treatment ----
    
    d <- master[treat == tr & is.finite(stake) & stake > 0]
    
    # ---- load treatment-specific (preferred) or pooled r draws ----
    
    f_r_tr   <- file.path(path_mod, paste0("mpl_r_draws_", tr, ".rds"))
    f_r_base <- file.path(path_mod, "mpl_r_draws.rds")
    f_r      <- if (file.exists(f_r_tr)) f_r_tr else f_r_base
    
    if (!file.exists(f_r))
      stop("Missing MPL r file for treatment ", tr,
           ". Looked for:\n  ", f_r_tr, "\n  ", f_r_base)
    
    r_obj <- readRDS(f_r)
    stopifnot(is.list(r_obj), !is.null(r_obj$pid), !is.null(r_obj$r_draws))
    
    r_pid   <- as.character(r_obj$pid)
    r_draws <- r_obj$r_draws
    stopifnot(is.matrix(r_draws), ncol(r_draws) == length(r_pid))
    
    # ---- handle no positive-stake trials -----------------------
    
    if (nrow(d) == 0L) {
      
      msg("  No positive-stake trials for treatment ", tr, "; writing empty outputs")
      
      trials_tbl <- empty_trials()
      pincl_tbl  <- empty_pid()
      pexcl_tbl  <- empty_pid()
      
      if (!skip_trials) { fwrite(trials_tbl, f_trials); msg("  Saved: ", f_trials) }
      if (!skip_pincl)  { fwrite(pincl_tbl,  f_pincl);  msg("  Saved: ", f_pincl)  }
      if (!skip_pexcl)  { fwrite(pexcl_tbl,  f_pexcl);  msg("  Saved: ", f_pexcl)  }
      if (!skip_rds) {
        saveRDS(list(
          data_folder = cfg$run$data_folder, treatment = tr,
          multiplier  = m, endowment = endowment,
          trials = trials_tbl, participants_incl = pincl_tbl,
          participants_excl = pexcl_tbl, p_hat_draws = NULL, flagged = NULL
        ), f_rds)
        msg("  Saved: ", f_rds)
      }
      
      outputs[[tr]] <- list(trials = trials_tbl,
                            participants_incl = pincl_tbl,
                            participants_excl = pexcl_tbl)
      next
    }
    
    # ---- restrict to participants with r draws -----------------
    
    d <- d[pid %in% r_pid]
    if (nrow(d) == 0L)
      stop("No positive-stake trials with matching MPL r draws for treatment ", tr, ".")
    
    # ================================================================
    # Vectorised implied-p computation (draws × trials)
    # ================================================================
    
    R     <- nrow(r_draws)
    J     <- nrow(d)
    r_mat <- r_draws[, match(d$pid, r_pid), drop = FALSE]  # R × J
    
    gain <- endowment + d$stake * (m - 1)   # J-length
    loss <- endowment - d$stake
    
    gain_mat  <- matrix(gain, nrow = R, ncol = J, byrow = TRUE)
    loss_mat  <- matrix(loss, nrow = R, ncol = J, byrow = TRUE)
    ratio_mat <- loss_mat / gain_mat
    
    # FOC solution: p* = 1 / (1 + (m-1) * (loss/gain)^r)
    raw <- 1 / (1 + (m - 1) * (ratio_mat ^ r_mat))
    
    # Flag draws where the FOC has no interior solution
    invalid <- !is.finite(r_mat) |
      gain_mat <= 0 | loss_mat <= 0 |
      !is.finite(raw)
    
    p_hat            <- pmin(pmax(raw, 0), 1)
    p_hat[loss_mat <= 0 & !is.finite(raw)] <- 1
    p_hat[gain_mat <= 0 & !is.finite(raw)] <- 0
    flagged          <- invalid | raw < 0 | raw > 1
    flagged[is.na(p_hat)] <- TRUE
    
    # ================================================================
    # Trial-level summaries
    # ================================================================
    
    trial_sum <- t(vapply(seq_len(J), function(j)
      summ_prob(p_hat[, j]), numeric(5L)))
    
    trials_tbl <- data.table(
      data_folder   = cfg$run$data_folder,
      treatment     = tr,
      pid           = as.character(d$pid),
      sequence      = as.character(d$seq),
      stake         = as.numeric(d$stake),
      n_draws       = R,
      n_flagged     = colSums(flagged, na.rm = TRUE),
      share_flagged = colMeans(flagged, na.rm = TRUE),
      median        = trial_sum[, "median"],
      mean          = trial_sum[, "mean"],
      q025          = trial_sum[, "q025"],
      q975          = trial_sum[, "q975"],
      p_gt_05       = trial_sum[, "p_gt_05"]
    )
    
    # ================================================================
    # Participant-level summaries (incl. and excl. flagged draws)
    # ================================================================
    
    pid_levels <- unique(as.character(d$pid))
    
    make_pid_tbl <- function(exclude_flagged) {
      rows <- vector("list", length(pid_levels))
      for (i in seq_along(pid_levels)) {
        pid_i    <- pid_levels[i]
        j_idx    <- which(d$pid == pid_i)
        vals_vec <- as.vector(p_hat[, j_idx, drop = FALSE])
        flg_vec  <- as.vector(flagged[, j_idx, drop = FALSE])
        vals_use <- if (exclude_flagged) vals_vec[!flg_vec] else vals_vec
        s        <- summ_prob(vals_use)
        rows[[i]] <- data.table(
          data_folder   = cfg$run$data_folder,
          treatment     = tr,
          pid           = pid_i,
          n_trials      = length(j_idx),
          n_draws_total = length(vals_vec),
          n_flagged     = sum(flg_vec, na.rm = TRUE),
          share_flagged = mean(flg_vec, na.rm = TRUE),
          median  = s["median"], mean    = s["mean"],
          q025    = s["q025"],   q975    = s["q975"],
          p_gt_05 = s["p_gt_05"]
        )
      }
      rbindlist(rows, fill = TRUE)
    }
    
    pincl_tbl <- make_pid_tbl(exclude_flagged = FALSE)
    pexcl_tbl <- make_pid_tbl(exclude_flagged = TRUE)
    
    # ================================================================
    # Treatment-level summary row (for comparison table)
    # ================================================================
    
    s_tr <- summ_prob(as.vector(p_hat))
    treat_rows[[tr]] <- data.table(
      data_folder = cfg$run$data_folder,
      treatment   = tr,
      median  = s_tr["median"], mean    = s_tr["mean"],
      q025    = s_tr["q025"],   q975    = s_tr["q975"],
      p_gt_05 = s_tr["p_gt_05"]
    )
    
    # ================================================================
    # Save outputs
    # ================================================================
    
    if (!skip_trials) { fwrite(trials_tbl, f_trials); msg("  Saved: ", f_trials) }
    if (!skip_pincl)  { fwrite(pincl_tbl,  f_pincl);  msg("  Saved: ", f_pincl)  }
    if (!skip_pexcl)  { fwrite(pexcl_tbl,  f_pexcl);  msg("  Saved: ", f_pexcl)  }
    if (!skip_rds) {
      saveRDS(list(
        data_folder       = cfg$run$data_folder,
        treatment         = tr,
        multiplier        = m,
        endowment         = endowment,
        trial_data        = d,
        trials            = trials_tbl,
        participants_incl = pincl_tbl,
        participants_excl = pexcl_tbl,
        p_hat_draws       = p_hat,
        flagged           = flagged
      ), f_rds)
      msg("  Saved: ", f_rds)
    }
    
    outputs[[tr]] <- list(trials = trials_tbl,
                          participants_incl = pincl_tbl,
                          participants_excl = pexcl_tbl)
  }
  
  # ----------------------------------------------------------
  # Treatment comparison table
  # ----------------------------------------------------------
  
  f_treat     <- file.path(path_out, "ex3_rq2_treatment_comparison.csv")
  skip_treat  <- should_skip(f_treat, cfg, "output",
                             "EX3 RQ2 treatment comparison")
  
  if (!skip_treat) {
    fwrite(rbindlist(treat_rows, fill = TRUE), f_treat)
    msg("Saved: ", f_treat)
  }
  
  invisible(outputs)
}