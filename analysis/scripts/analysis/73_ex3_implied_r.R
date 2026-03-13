# ============================================================
# scripts/analysis/73_rq2_implied_r.R
#
# Implied risk aversion from observed stakes under three belief
# assumptions:
#   1) objective p = 0.5
#   2) RQ4 baseline p = mean(hbar)
#   3) self-reported p = perceived_win_probability / 100
#
# Outputs per treatment:
#   data/clean/<ds>/output/ex3_implied_r_participants_<tr>.csv
#   data/clean/<ds>/output/ex3_implied_r_summary_<tr>.csv
#
# Notes:
# - No trial-level table is written.
# - Participant table compares MPL r with participant-level implied r
#   summaries under all three belief assumptions.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex3_implied_r <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$treatment))
  
  ds <- as.character(cfg$run$dataset)
  
  mod_dir   <- path_mod_ds(ds)
  out_dir   <- path_out_ds(ds)
  clean_dir <- path_clean_ds(ds)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(clean_dir, "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  master <- fread(f_master, encoding = "UTF-8")
  
  req <- c("pid", "treat", "seq", "stake")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0) {
    stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  }
  
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  
  if ("perceived_win_probability" %in% names(master)) {
    master[, perceived_win_probability := as.numeric(perceived_win_probability)]
  }
  
  endowment <- as.numeric(cfg$design$seq$endowment)
  stopifnot(length(endowment) == 1L, is.finite(endowment), endowment > 0)
  
  tr_cfg  <- cfg$design$seq$treatments
  tr_mult <- vapply(tr_cfg, as.numeric, numeric(1))
  tr_run  <- unique(as.character(cfg$run$treatment))
  
  stopifnot(all(tr_run %in% names(tr_mult)))
  
  summ <- function(x) {
    if (!length(x) || all(!is.finite(x))) {
      return(c(
        median = NA_real_,
        mean   = NA_real_,
        q025   = NA_real_,
        q975   = NA_real_
      ))
    }
    q <- quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    c(
      median = median(x, na.rm = TRUE),
      mean   = mean(x, na.rm = TRUE),
      q025   = q[1],
      q975   = q[2]
    )
  }
  
  # implied r from FOC:
  # r = log( p(m-1)/(1-p) ) / log( (W + a(m-1)) / (W-a) )
  calc_r <- function(p, m, gain, loss) {
    
    denom <- log(gain / loss)
    
    if (length(p) == 1L) {
      p <- rep(p, length(denom))
    }
    
    out <- rep(NA_real_, length(denom))
    
    ok <- is.finite(p) &
      p > 0 & p < 1 &
      is.finite(gain) & gain > 0 &
      is.finite(loss) & loss > 0 &
      is.finite(denom) & denom != 0
    
    if (any(ok)) {
      num <- log((p[ok] * (m - 1)) / (1 - p[ok]))
      out[ok] <- num / denom[ok]
    }
    
    out
  }
  
  outputs <- list()
  
  for (tr in tr_run) {
    
    msg("Running implied r for treatment: ", tr)
    
    m <- as.numeric(tr_mult[[tr]])
    stopifnot(is.finite(m), m > 1)
    
    f_part <- file.path(out_dir, paste0("ex3_implied_r_participants_", tr, ".csv"))
    f_sum  <- file.path(out_dir, paste0("ex3_implied_r_summary_", tr, ".csv"))
    
    skip_all <- should_skip(
      paths = c(f_part, f_sum),
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 implied r (", ds, "/", tr, ")")
    )
    if (skip_all) next
    
    d <- master[treat == tr & is.finite(stake) & stake > 0]
    if (nrow(d) == 0L) next
    
    # --------------------------------------------------------
    # MPL r draws: prefer treatment-specific file
    # --------------------------------------------------------
    
    f_r_tr   <- file.path(mod_dir, paste0("mpl_r_draws_", tr, ".rds"))
    f_r_base <- file.path(mod_dir, "mpl_r_draws.rds")
    
    f_r <- if (file.exists(f_r_tr)) f_r_tr else f_r_base
    if (!file.exists(f_r)) {
      stop("Missing MPL r file for treatment ", tr,
           ". Looked for: ", f_r_tr, " and ", f_r_base)
    }
    
    r_obj <- readRDS(f_r)
    stopifnot(is.list(r_obj), !is.null(r_obj$pid), !is.null(r_obj$r_draws))
    
    r_pid   <- as.character(r_obj$pid)
    r_draws <- r_obj$r_draws
    
    stopifnot(is.matrix(r_draws), ncol(r_draws) == length(r_pid))
    
    # only participants with MPL draws
    d <- d[pid %in% r_pid]
    if (nrow(d) == 0L) {
      stop("No positive-stake trials have matching MPL r draws for treatment ", tr, ".")
    }
    
    # --------------------------------------------------------
    # RQ4 baseline p = mean(hbar)
    # prefer full fit if available
    # --------------------------------------------------------
    
    p_rq4 <- NA_real_
    
    f_rq4_full <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_rq4_base <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    f_rq4 <- if (file.exists(f_rq4_full)) f_rq4_full else f_rq4_base
    
    if (file.exists(f_rq4)) {
      post_rq4 <- rstan::extract(readRDS(f_rq4))
      if (!is.null(post_rq4$hbar)) {
        p_rq4 <- mean(post_rq4$hbar, na.rm = TRUE)
      }
    }
    
    # --------------------------------------------------------
    # self-report p
    # --------------------------------------------------------
    
    if ("perceived_win_probability" %in% names(master)) {
      self_tbl <- unique(
        master[treat == tr, .(pid, perceived_win_probability)]
      )
      self_tbl[, p_self := perceived_win_probability / 100]
      d <- merge(d, self_tbl[, .(pid, p_self)], by = "pid", all.x = TRUE)
    } else {
      d[, p_self := NA_real_]
    }
    
    # --------------------------------------------------------
    # implied r calculations
    # --------------------------------------------------------
    
    gain <- endowment + d$stake * (m - 1)
    loss <- endowment - d$stake
    
    r_obj05 <- calc_r(0.5,   m, gain, loss)
    r_rq4   <- calc_r(p_rq4, m, gain, loss)
    r_self  <- calc_r(d$p_self, m, gain, loss)
    
    d_imp <- data.table(
      pid = d$pid,
      seq = d$seq,
      stake = d$stake,
      p_self = d$p_self,
      r_obj05 = r_obj05,
      r_rq4   = r_rq4,
      r_self  = r_self
    )
    
    # --------------------------------------------------------
    # participant-level comparison table
    # --------------------------------------------------------
    
    pid_levels <- unique(as.character(d_imp$pid))
    
    r_idx <- match(pid_levels, r_pid)
    mpl_draws_pid <- r_draws[, r_idx, drop = FALSE]
    
    mpl_stats <- t(vapply(seq_along(pid_levels), function(j) {
      summ(mpl_draws_pid[, j])
    }, numeric(4L)))
    
    part_rows <- vector("list", length(pid_levels))
    
    for (i in seq_along(pid_levels)) {
      
      pid_i <- pid_levels[i]
      dd    <- d_imp[pid == pid_i]
      
      s_obj  <- summ(dd$r_obj05)
      s_rq4  <- summ(dd$r_rq4)
      s_self <- summ(dd$r_self)
      
      part_rows[[i]] <- data.table(
        dataset = ds,
        treatment = tr,
        pid = pid_i,
        n_trials = nrow(dd),
        
        mpl_r_median = mpl_stats[i, "median"],
        mpl_r_mean   = mpl_stats[i, "mean"],
        mpl_r_q025   = mpl_stats[i, "q025"],
        mpl_r_q975   = mpl_stats[i, "q975"],
        
        r_obj05_median = s_obj["median"],
        r_obj05_mean   = s_obj["mean"],
        r_obj05_q025   = s_obj["q025"],
        r_obj05_q975   = s_obj["q975"],
        
        r_rq4_median = s_rq4["median"],
        r_rq4_mean   = s_rq4["mean"],
        r_rq4_q025   = s_rq4["q025"],
        r_rq4_q975   = s_rq4["q975"],
        
        r_self_median = s_self["median"],
        r_self_mean   = s_self["mean"],
        r_self_q025   = s_self["q025"],
        r_self_q975   = s_self["q975"]
      )
    }
    
    part_tbl <- rbindlist(part_rows, fill = TRUE)
    
    # --------------------------------------------------------
    # treatment-level summary across participants
    # --------------------------------------------------------
    
    summary_tbl <- rbindlist(list(
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "mpl_r_median",
        t(summ(part_tbl$mpl_r_median))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "mpl_r_mean",
        t(summ(part_tbl$mpl_r_mean))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "r_obj05_median",
        t(summ(part_tbl$r_obj05_median))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "r_obj05_mean",
        t(summ(part_tbl$r_obj05_mean))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "r_rq4_median",
        t(summ(part_tbl$r_rq4_median))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "r_rq4_mean",
        t(summ(part_tbl$r_rq4_mean))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "r_self_median",
        t(summ(part_tbl$r_self_median))
      ),
      data.table(
        dataset = ds,
        treatment = tr,
        measure = "r_self_mean",
        t(summ(part_tbl$r_self_mean))
      )
    ), fill = TRUE)
    
    fwrite(part_tbl, f_part)
    fwrite(summary_tbl, f_sum)
    
    msg("Saved: ", f_part)
    msg("Saved: ", f_sum)
    
    outputs[[tr]] <- list(
      participants = part_tbl,
      summary = summary_tbl
    )
  }
  
  invisible(outputs)
}

# Example:
# ex3_implied_r(cfg)