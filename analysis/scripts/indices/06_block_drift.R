# ============================================================
# 06_block_drift.R
#
# PURPOSE
#   Runs within-participant permutation tests for drift in betting
#   rate, stake size, and side choice across the four experimental
#   blocks. Results inform drift adjustment decisions for Stan models.
#
# INPUT
#   path_src/master_sequences.csv
#   path_out/mpl_scored_<tr>.csv       -- for consistent subset
#   path_mod/a_star_pid_flags_<tr>.rds -- for conf subset
#
# OUTPUT
#   path_mod/drift_decisions.rds  -- full results list
#   path_out/drift_decisions.csv  -- combined summary table
#
# OUTCOMES TESTED
#   bet   -- betting rate by block: y = 1(stake > 0)
#   stake -- mean stake by block (conditional on betting)
#   side  -- Heads rate by block (conditional on betting)
#
# NOTES
#   - Drift type: "linear" if delta test triggers, "categorical"
#     if only range triggers, "none" otherwise
#   - Results saved to drift_decisions.rds and read at runtime
#     by Stan scripts (rq1_stan, rq2_stan, rq4_stan)
#   - Drift parametrization (prior SD, block centers) is in
#     cfg$design$drift$params
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

block_drift_check <- function(cfg) {
  
  tr_vec      <- unique(as.character(cfg$run$treatment))
  B           <- 2000L
  alpha       <- 0.05
  seed        <- as.integer(cfg$run$seed)
  heads_label <- as.character(cfg$design$seq$side_labels$heads)
  tau_main    <- as.numeric(cfg$design$a_flags$tau[1])
  tau_nm      <- gsub("\\.", "", sprintf("%.2f", tau_main))
  
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt0 <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "block", "side") %in% names(dt0)))
  
  dt0[, pid   := as.character(pid)]
  dt0[, treat := as.character(treat)]
  dt0[, stake := as.numeric(stake)]
  dt0[, block := as.integer(block)]
  stopifnot(all(dt0$block %in% 1:4))
  
  dt0[, y_bet   := bet_indicator(stake)]
  dt0[, block_c := block_centered(block)]
  
  f_rds <- file.path(path_mod, "drift_decisions.rds")
  f_csv <- file.path(path_out, "drift_decisions.csv")
  
  if (should_skip(
    paths = c(f_rds, f_csv),
    cfg   = cfg,
    type  = "model",
    label = "Drift decisions"
  )) return(invisible(NULL))
  
  compute_one <- function(d, tr, tag, outcome) {
    
    d_out <- switch(outcome,
                    bet = d[, .(pid, block, y = y_bet)],
                    stake = {
                      d_bet <- d[y_bet == 1L, .(pid, block, y = stake)]
                      if (nrow(d_bet) == 0) {
                        warning("No betting trials for stake drift (tr=", tr, ", tag=", tag, ")")
                        return(NULL)
                      }
                      d_bet
                    },
                    side = {
                      d_bet <- d[y_bet == 1L & !is.na(side),
                                 .(pid, block, y = as.integer(side == heads_label))]
                      if (nrow(d_bet) == 0) {
                        warning("No betting trials for side drift (tr=", tr, ", tag=", tag, ")")
                        return(NULL)
                      }
                      d_bet
                    },
                    stop("Unknown outcome: ", outcome)
    )
    
    d_out[, block := as.integer(block)]
    if (nrow(d_out) == 0) return(NULL)
    
    res        <- perm_test_within_pid(d_out, B = B, seed = seed)
    drift_flag <- isTRUE(res$p_delta_two_sided < alpha) || isTRUE(res$p_range_one_sided < alpha)
    drift_type <- if (isTRUE(res$p_delta_two_sided < alpha)) "linear" else
      if (isTRUE(res$p_range_one_sided  < alpha)) "categorical" else "none"
    
    bm <- sprintf("b1=%.3f b2=%.3f b3=%.3f b4=%.3f",
                  res$block_means[1], res$block_means[2],
                  res$block_means[3], res$block_means[4])
    
    msg(if (drift_flag) "Drift detected" else "No drift detected",
        " (", tr, "/", tag, "/", outcome, "):  ",
        "Δ41=", sprintf("%.3f", res$delta_41),
        ", pΔ=", sprintf("%.4f", res$p_delta_two_sided),
        " | range=", sprintf("%.3f", res$range_b),
        ", pR=", sprintf("%.4f", res$p_range_one_sided),
        " | ", bm, " | recommend=", drift_type)
    
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
      if (!is.null(res)) {
        key <- paste(tr, tag, outcome, sep = "_")
        all_results[[key]] <<- res$full
        all_rows[[key]]    <<- res$row
      }
    }
  }
  
  for (tr in tr_vec) {
    
    run_all_outcomes(dt0[treat == tr], tr, "full")
    
    f_scored <- file.path(path_out, paste0("mpl_scored_", tr, ".csv"))
    if (file.exists(f_scored)) {
      scored <- fread(f_scored, encoding = "UTF-8")
      scored[, pid := as.character(pid)]
      run_all_outcomes(
        dt0[treat == tr & pid %in% scored[inconsistent == 0L, pid]],
        tr, "consistent"
      )
    } else {
      warning("Scored file not found for tr=", tr, "; skipping consistent drift check.")
    }
    
    if (!isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
    
    f_rds_flags <- file.path(path_mod, paste0("a_star_pid_flags_", tr, ".rds"))
    stopifnot(file.exists(f_rds_flags))
    flags   <- readRDS(f_rds_flags)
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    run_all_outcomes(
      dt0[treat == tr & pid %in% sort(unique(as.character(flags$pid_sets[[nm_keep]])))],
      tr, "conf"
    )
  }
  
  saveRDS(all_results, f_rds)
  fwrite(rbindlist(all_rows), f_csv)
  msg("Saved: ", f_rds)
  msg("Saved: ", f_csv)
  
  invisible(all_results)
}