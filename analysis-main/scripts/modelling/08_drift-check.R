# ============================================================
# 08_drift-check.R
#
# PURPOSE
#   Tests for drift across the four experimental blocks in betting
#   rate, stake size and side choice, with within-participant
#   permutation tests, and records whether a drift term enters the
#   Stan models. For each sample and outcome the block means are
#   compared with replicates in which the block labels are permuted
#   within participant:
#     delta_41  difference between block 4 and block 1 (two-sided p)
#     range_b   range of the four block means (one-sided p)
#   Drift is flagged if either p-value is below alpha (0.05); the type is
#   "linear" if the delta test triggers, "categorical" if only the range
#   test triggers, and "none" otherwise. The number of permutations
#   (2000) and alpha are fixed in the script; the drift priors are in
#   cfg$design$drift$params.
#
#   The script works on the MPL model selected under the key "mpl" in
#   selected-models.csv and writes into a subfolder of that model, because
#   the confirmatory sample depends on its a* classification.
#
# SAMPLES (within each treatment)
#   full          all participants
#   consistent    Holt-Laury-consistent participants (when
#                 cfg$run$consistent_only is TRUE)
#   confirmatory  participants kept by the a* classification (treatments
#                 with cfg$design$a_flags$betting_normative = TRUE)
#
# OUTCOMES
#   bet    betting indicator, 1(stake > 0), all trials
#   stake  stake, trials with a bet
#   side   Heads indicator, trials with a bet
#
# INPUT
#   path_der/master-sequences.csv
#   path_der/switching-point-r.csv       -- inconsistency flag per pid
#   model "mpl": a-star/a-star-flags-<tr>-<sample>.csv  (via a_star_pid_set)
#
# OUTPUT (subfolder drift/ of the selected MPL model folder)
#   drift-decisions.rds   named list, one entry per <tr>_<sample>_<outcome>
#                         (read by the Stan model scripts)
#   drift-decisions.csv   the same results as a table
#   diagnostics.csv       integrity checks: trials per test and blocks present
#   inputs.csv, run-info.txt
#
# NOTES
#   - The subfolder is built as drift.tmp and renamed when complete; an
#     existing drift folder is never overwritten (delete it to rerun)
#   - A sample or outcome without observations stops the script
# ============================================================

block_centered <- function(block_int) {
  stopifnot(all(block_int %in% 1:4))
  as.numeric(block_int) - 2.5
}

bet_indicator <- function(stake) {
  as.integer(as.numeric(stake) > 0)
}

compute_block_means <- function(y, block) {
  out <- tapply(y, block, mean)
  m   <- rep(NA_real_, 4)
  names(m) <- paste0("b", 1:4)
  if (!is.null(out)) {
    for (b in names(out)) {
      bb <- as.integer(b)
      if (!is.na(bb) && bb >= 1 && bb <= 4) m[bb] <- as.numeric(out[[b]])
    }
  }
  m
}

safe_range <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x) - min(x)
}

perm_test_within_pid <- function(dt, B = 2000L, seed = 1L) {
  stopifnot(all(c("pid", "y", "block") %in% names(dt)))
  stopifnot(is.integer(dt$block))
  stopifnot(all(dt$block %in% 1:4))
  
  dt <- data.table::copy(dt)
  dt[, y := as.numeric(y)]
  
  m_obs     <- compute_block_means(dt$y, dt$block)
  delta_obs <- as.numeric(m_obs[4] - m_obs[1])
  range_obs <- as.numeric(safe_range(m_obs))
  
  dt[, row_id := .I]
  idx_by_pid <- split(dt$row_id, dt$pid)
  
  set.seed(seed)
  delta_rep <- numeric(B)
  range_rep <- numeric(B)
  
  for (b in seq_len(B)) {
    block_perm <- dt$block
    for (rows in idx_by_pid) {
      if (length(rows) <= 1L) next
      block_perm[rows] <- sample(block_perm[rows], size = length(rows), replace = FALSE)
    }
    m_rep        <- compute_block_means(dt$y, block_perm)
    delta_rep[b] <- as.numeric(m_rep[4] - m_rep[1])
    range_rep[b] <- as.numeric(safe_range(m_rep))
  }
  
  list(
    block_means       = m_obs,
    delta_41          = delta_obs,
    range_b           = range_obs,
    p_delta_two_sided = as.numeric(mean(abs(delta_rep) >= abs(delta_obs))),
    p_range_one_sided = as.numeric(mean(range_rep >= range_obs))
  )
}

drift_check <- function(cfg) {
  
  mpl_dir <- model_path("mpl")
  run     <- model_begin(file.path(basename(mpl_dir), "drift"))
  if (is.null(run)) return(invisible(NULL))
  
  tr_vec          <- unique(as.character(cfg$run$treatment))
  consistent_only <- isTRUE(cfg$run$consistent_only)
  B               <- 2000L
  alpha           <- 0.05
  seed            <- as.integer(cfg$run$seed)
  heads_label     <- as.character(cfg$design$seq$side_labels$heads)
  tau_main        <- as.numeric(cfg$design$a_flags$tau[1])
  a_sample        <- hl_sample(cfg)
  
  # ---- Trials and the Holt-Laury consistency flag ----
  f_master <- file.path(path_der, "master-sequences.csv")
  f_sw     <- file.path(path_der, "switching-point-r.csv")
  stopifnot(file.exists(f_master), file.exists(f_sw))
  
  dt0 <- fread(f_master, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "block", "side") %in% names(dt0)))
  dt0[, pid   := as.character(pid)]
  dt0[, treat := as.character(treat)]
  dt0[, stake := as.numeric(stake)]
  dt0[, block := as.integer(block)]
  stopifnot(all(dt0$block %in% 1:4))
  dt0[, y_bet   := bet_indicator(stake)]
  dt0[, block_c := block_centered(block)]
  
  sw <- fread(f_sw, encoding = "UTF-8", colClasses = c(pid = "character"))
  stopifnot(all(c("pid", "inconsistent") %in% names(sw)))
  consistent_pids <- sw[inconsistent == 0L, pid]
  
  # ---- One permutation test: sample x outcome ----
  compute_one <- function(d, tr, tag, outcome) {
    
    d_out <- switch(outcome,
                    bet   = d[, .(pid, block, y = y_bet)],
                    stake = d[y_bet == 1L, .(pid, block, y = stake)],
                    side  = d[y_bet == 1L & !is.na(side),
                              .(pid, block, y = as.integer(side == heads_label))],
                    stop("Unknown outcome: ", outcome))
    
    if (nrow(d_out) == 0L) {
      stop("No observations for ", outcome, " drift (", tr, "/", tag, ").", call. = FALSE)
    }
    d_out[, block := as.integer(block)]
    
    res        <- perm_test_within_pid(d_out, B = B, seed = seed)
    drift_flag <- isTRUE(res$p_delta_two_sided < alpha) || isTRUE(res$p_range_one_sided < alpha)
    drift_type <- if (isTRUE(res$p_delta_two_sided < alpha)) "linear" else
      if (isTRUE(res$p_range_one_sided < alpha)) "categorical" else "none"
    
    list(
      full = list(
        treatment = tr, tag = tag, outcome = outcome,
        block_means = res$block_means, delta_41 = res$delta_41,
        range_b = res$range_b, p_delta_two_sided = res$p_delta_two_sided,
        p_range_one_sided = res$p_range_one_sided,
        drift = drift_flag, drift_type = drift_type
      ),
      row = data.table(
        treatment = tr, tag = tag, outcome = outcome, n = nrow(d_out),
        n_blocks = uniqueN(d_out$block),
        b1 = res$block_means[1], b2 = res$block_means[2],
        b3 = res$block_means[3], b4 = res$block_means[4],
        delta_41 = res$delta_41, range_b = res$range_b,
        p_delta_two_sided = res$p_delta_two_sided,
        p_range_one_sided = res$p_range_one_sided,
        drift = drift_flag, drift_type = drift_type
      )
    )
  }
  
  all_results <- list()
  all_rows    <- list()
  
  run_all_outcomes <- function(d, tr, tag) {
    for (outcome in c("bet", "stake", "side")) {
      res <- compute_one(d, tr, tag, outcome)
      key <- paste(tr, tag, outcome, sep = "_")
      all_results[[key]] <<- res$full
      all_rows[[key]]    <<- res$row
    }
  }
  
  for (tr in tr_vec) {
    
    run_all_outcomes(dt0[treat == tr], tr, "full")
    
    if (consistent_only) {
      run_all_outcomes(dt0[treat == tr & pid %in% consistent_pids], tr, "consistent")
    }
    
    if (!isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
    
    keep_pid <- a_star_pid_set(tr, a_sample, tau_main, "keep")
    run_all_outcomes(dt0[treat == tr & pid %in% keep_pid], tr, "confirmatory")
  }
  
  res_tbl <- rbindlist(all_rows)
  saveRDS(all_results, file.path(run, "drift-decisions.rds"))
  fwrite(res_tbl[, !"n_blocks"], file.path(run, "drift-decisions.csv"))
  
  # ---- Diagnostics: every test has observations in all four blocks ----
  write_diagnostics(run, rbindlist(list(
    res_tbl[, .(tag = paste(treatment, tag, sep = "_"), treatment, sample = tag,
                check = paste0("n_trials_", outcome), value = n,
                threshold = 1, pass = n >= 1)],
    res_tbl[, .(tag = paste(treatment, tag, sep = "_"), treatment, sample = tag,
                check = paste0("blocks_present_", outcome), value = n_blocks,
                threshold = 4, pass = n_blocks == 4L)]
  )))
  
  model_end(run, cfg, settings = list(
    permutations = B, alpha = alpha, seed = seed, tau = tau_main,
    consistent_only = consistent_only
  ))
}
