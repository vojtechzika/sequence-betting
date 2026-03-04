# ------------------------------------------------------------
# 03_drift_check.R
# Central drift diagnostics (dataset x treatment x tag)
#
# Output: drift_decision_<tr>_<tag>.rds saved into path_mod_ds(ds)
#
# Drift diagnostics are model-agnostic and based on betting rate by block:
#   y = 1(stake > 0)
#
# Two drift statistics (both within-participant permutation tested):
#   1) delta_41  = mean(y|block=4) - mean(y|block=1)        [monotone drift]
#   2) range_b   = max_b mean(y|block=b) - min_b mean(y|b)  [non-monotone drift]
#
# Decision: drift = (p_delta_two_sided < alpha) OR (p_range_one_sided < alpha)
#
# Drift-type recommendation (prereg-aligned):
#   - "linear"       if delta test triggers (monotone)
#   - "categorical"  if only range triggers (non-monotone)
#   - "none"         otherwise
#
# Requirements:
# - master_sequences.csv contains: pid, treat, stake, block
# - a_star_pid_flags_<tr>.rds exists for conf tags when betting_normative[[tr]] is TRUE
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# --- helpers -------------------------------------------------

block_centered <- function(block_int) {
  # maps {1,2,3,4} -> {-1.5,-0.5,0.5,1.5}
  stopifnot(all(block_int %in% 1:4))
  as.numeric(block_int) - 2.5
}

bet_indicator <- function(stake) {
  # stake is always numeric per your invariant (0 = no bet)
  as.integer(as.numeric(stake) > 0)
}

compute_block_means <- function(y, block) {
  out <- tapply(y, block, mean)
  m <- rep(NA_real_, 4)
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
  # dt must contain columns: pid, y, block
  stopifnot(all(c("pid", "y", "block") %in% names(dt)))
  stopifnot(is.integer(dt$y), is.integer(dt$block))
  stopifnot(all(dt$block %in% 1:4))
  
  dt <- data.table::copy(dt)
  
  # observed stats
  m_obs <- compute_block_means(dt$y, dt$block)
  delta_obs <- as.numeric(m_obs[4] - m_obs[1])
  range_obs <- as.numeric(safe_range(m_obs))
  
  # pre-split indices by pid for speed
  dt[, row_id := .I]
  idx_by_pid <- split(dt$row_id, dt$pid)
  
  # permutation distribution
  set.seed(seed)
  delta_rep <- numeric(B)
  range_rep <- numeric(B)
  
  for (b in seq_len(B)) {
    block_perm <- dt$block
    for (rows in idx_by_pid) {
      if (length(rows) <= 1L) next
      block_perm[rows] <- sample(block_perm[rows], size = length(rows), replace = FALSE)
    }
    
    m_rep <- compute_block_means(dt$y, block_perm)
    delta_rep[b] <- as.numeric(m_rep[4] - m_rep[1])
    range_rep[b] <- as.numeric(safe_range(m_rep))
  }
  
  # p-values
  p_delta_two <- mean(abs(delta_rep) >= abs(delta_obs))
  # range is non-negative; one-sided is appropriate
  p_range_one <- mean(range_rep >= range_obs)
  
  list(
    block_means = m_obs,
    delta_41 = delta_obs,
    range_b = range_obs,
    delta_rep = delta_rep,
    range_rep = range_rep,
    p_delta_two_sided = as.numeric(p_delta_two),
    p_range_one_sided = as.numeric(p_range_one)
  )
}

# --- main ----------------------------------------------------

drift_check <- function(cfg) {
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  B <- 2000L
  alpha <- 0.05
  seed <- as.integer(cfg$run$seed)
  
  # explicit main threshold
  tau_main <- 0.90
  tau_nm <- gsub("\\.", "", sprintf("%.2f", tau_main)) # "090"
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt0 <- data.table::fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "block") %in% names(dt0)))
  
  dt0[, pid := as.character(pid)]
  dt0[, treat := as.character(treat)]
  dt0[, stake := as.numeric(stake)]
  dt0[, block := as.integer(block)]
  stopifnot(all(dt0$block %in% 1:4))
  
  dt0[, y := bet_indicator(stake)]
  dt0[, block_c := block_centered(block)]
  
  mod_dir <- path_mod_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  compute_one <- function(d, tr, tag) {
    
    f_out <- file.path(mod_dir, paste0("drift_decision_", tr, "_", tag, ".rds"))
    
    if (should_skip(
      paths = f_out,
      cfg   = cfg,
      type  = "model",
      label = paste0("Drift decision (", ds, "/", tr, "/", tag, ")")
    )) return(invisible(NULL))
    
    if (nrow(d) == 0) return(invisible(NULL))
    
    res <- perm_test_within_pid(d[, .(pid, y, block)], B = B, seed = seed)
    
    drift_flag <- isTRUE(res$p_delta_two_sided < alpha) || isTRUE(res$p_range_one_sided < alpha)
    
    # prereg-aligned recommendation:
    drift_type <- "none"
    if (isTRUE(res$p_delta_two_sided < alpha)) {
      drift_type <- "linear"
    } else if (isTRUE(res$p_range_one_sided < alpha)) {
      drift_type <- "categorical"
    }
    
    out <- list(
      dataset = ds,
      treatment = tr,
      tag = tag,
      method = "within_pid_permutation_on_betting_rate",
      B = B,
      alpha = alpha,
      block_means = res$block_means,
      delta_41 = res$delta_41,
      range_b = res$range_b,
      p_delta_two_sided = res$p_delta_two_sided,
      p_range_one_sided = res$p_range_one_sided,
      drift = drift_flag,
      drift_type_recommended = drift_type,
      block_centering = c("1" = -1.5, "2" = -0.5, "3" = 0.5, "4" = 1.5),
      decision_rule = "drift = (p_delta_two_sided < alpha) OR (p_range_one_sided < alpha)"
    )
    
    saveRDS(out, f_out)
    
    # ------------------------------------------------------------
    # Console message
    # ------------------------------------------------------------
    bm <- sprintf("b1=%.3f b2=%.3f b3=%.3f b4=%.3f",
                  out$block_means[1],
                  out$block_means[2],
                  out$block_means[3],
                  out$block_means[4])
    
    if (drift_flag) {
      msg(
        "Drift detected (", ds, "/", tr, "/", tag, "): ",
        "Δ41=", sprintf("%.3f", out$delta_41),
        ", pΔ=", sprintf("%.4f", out$p_delta_two_sided),
        " | range=", sprintf("%.3f", out$range_b),
        ", pR=", sprintf("%.4f", out$p_range_one_sided),
        " | ", bm,
        " | recommend=", drift_type
      )
    } else {
      msg(
        "No drift detected (", ds, "/", tr, "/", tag, "): ",
        "Δ41=", sprintf("%.3f", out$delta_41),
        ", pΔ=", sprintf("%.4f", out$p_delta_two_sided),
        " | range=", sprintf("%.3f", out$range_b),
        ", pR=", sprintf("%.4f", out$p_range_one_sided),
        " | ", bm,
        " | recommend=", drift_type
      )
    }
    
    msg("Saved: ", f_out)
    
    invisible(out)
  }
  
  # FULL + CONF (per treatment)
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # FULL
    d_full <- dt0[treat == tr]
    outputs[[paste0(tr, "_full")]] <- compute_one(d_full, tr, "full")
    
    # CONF only if betting_normative TRUE for that treatment
    if (!isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
    
    f_rds_flags <- file.path(mod_dir, paste0("a_star_pid_flags_", tr, ".rds"))
    stopifnot(file.exists(f_rds_flags))
    
    flags <- readRDS(f_rds_flags)
    stopifnot(is.list(flags), !is.null(flags$pid_sets))
    
    nm_keep <- paste0("pid_keep_tau", tau_nm)
    stopifnot(nm_keep %in% names(flags$pid_sets))
    
    keep_pid <- sort(unique(as.character(flags$pid_sets[[nm_keep]])))
    
    d_conf <- dt0[treat == tr & pid %in% keep_pid]
    outputs[[paste0(tr, "_conf")]] <- compute_one(d_conf, tr, "conf")
  }
  
  invisible(outputs)
}

# If you source this file, run:
# drift_check(cfg)