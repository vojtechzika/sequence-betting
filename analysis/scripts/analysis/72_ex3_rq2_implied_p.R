# ============================================================
# scripts/analysis/72_ex3_rq2_implied_prob.R
#
# EX3 / RQ2:
# Implied subjective probabilities from positive stakes
# Runs for ALL treatments in cfg$run$treatment
#
# Outputs per treatment:
#   data/clean/<ds>/output/ex3_rq2_trials_<tr>.csv
#   data/clean/<ds>/output/ex3_rq2_participants_incl_<tr>.csv
#   data/clean/<ds>/output/ex3_rq2_participants_excl_<tr>.csv
#   data/clean/<ds>/models/ex3_rq2_<tr>.rds
#
# Additional output:
#   data/clean/<ds>/output/ex3_rq2_treatment_comparison.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

ex3_rq2_implied_p <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$treatment))
  
  ds <- as.character(cfg$run$dataset)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  
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
  
  endowment <- as.numeric(cfg$design$seq$endowment)
  stopifnot(length(endowment) == 1L, is.finite(endowment), endowment > 0)
  
  tr_cfg  <- cfg$design$seq$treatments
  tr_mult <- vapply(tr_cfg, as.numeric, numeric(1))
  tr_run  <- unique(as.character(cfg$run$treatment))
  
  stopifnot(all(tr_run %in% names(tr_mult)))
  
  summ_prob <- function(x) {
    
    if (!length(x) || all(!is.finite(x))) {
      return(c(
        median  = NA_real_,
        mean    = NA_real_,
        q025    = NA_real_,
        q975    = NA_real_,
        p_gt_05 = NA_real_
      ))
    }
    
    q <- quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    
    c(
      median  = median(x, na.rm = TRUE),
      mean    = mean(x, na.rm = TRUE),
      q025    = q[1],
      q975    = q[2],
      p_gt_05 = mean(x > 0.5, na.rm = TRUE)
    )
  }
  
  outputs <- list()
  treat_rows <- list()
  
  for (tr in tr_run) {
    
    msg("Running implied probability analysis for treatment: ", tr)
    
    m <- as.numeric(tr_mult[[tr]])
    stopifnot(is.finite(m), m > 1)
    
    f_trials <- file.path(out_dir, paste0("ex3_rq2_trials_", tr, ".csv"))
    f_pincl  <- file.path(out_dir, paste0("ex3_rq2_participants_incl_", tr, ".csv"))
    f_pexcl  <- file.path(out_dir, paste0("ex3_rq2_participants_excl_", tr, ".csv"))
    f_rds    <- file.path(mod_dir, paste0("ex3_rq2_", tr, ".rds"))
    
    skip_all <- should_skip(
      paths = c(f_trials, f_pincl, f_pexcl, f_rds),
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 RQ2 implied probabilities (", ds, "/", tr, ")")
    )
    
    d <- master[treat == tr & is.finite(stake) & stake > 0]
    
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
    
    if (nrow(d) == 0L) {
      
      trials_tbl <- data.table(
        dataset = character(),
        treatment = character(),
        pid = character(),
        sequence = character(),
        stake = numeric(),
        n_draws = integer(),
        n_flagged = integer(),
        share_flagged = numeric(),
        median = numeric(),
        mean = numeric(),
        q025 = numeric(),
        q975 = numeric(),
        p_gt_05 = numeric()
      )
      
      pincl_tbl <- data.table(
        dataset = character(),
        treatment = character(),
        pid = character(),
        n_trials = integer(),
        n_draws_total = integer(),
        n_flagged = integer(),
        share_flagged = numeric(),
        median = numeric(),
        mean = numeric(),
        q025 = numeric(),
        q975 = numeric(),
        p_gt_05 = numeric()
      )
      
      pexcl_tbl <- copy(pincl_tbl)
      
      if (!skip_all) {
        fwrite(trials_tbl, f_trials)
        fwrite(pincl_tbl,  f_pincl)
        fwrite(pexcl_tbl,  f_pexcl)
        saveRDS(list(
          dataset = ds,
          treatment = tr,
          multiplier = m,
          endowment = endowment,
          trials = trials_tbl,
          participants_incl = pincl_tbl,
          participants_excl = pexcl_tbl,
          p_hat_draws = NULL,
          flagged = NULL
        ), f_rds)
        
        msg("Saved: ", f_trials)
        msg("Saved: ", f_pincl)
        msg("Saved: ", f_pexcl)
        msg("Saved: ", f_rds)
      }
      
      outputs[[tr]] <- list(
        trials = trials_tbl,
        participants_incl = pincl_tbl,
        participants_excl = pexcl_tbl
      )
      
      next
    }
    
    d <- d[pid %in% r_pid]
    if (nrow(d) == 0L) {
      stop("No positive-stake trials with matching MPL r draws for treatment ", tr, ".")
    }
    
    idx_pid <- match(d$pid, r_pid)
    
    R <- nrow(r_draws)
    J <- nrow(d)
    
    r_mat <- r_draws[, idx_pid, drop = FALSE]
    
    gain <- endowment + d$stake * (m - 1)
    loss <- endowment - d$stake
    
    gain_mat  <- matrix(gain, nrow = R, ncol = J, byrow = TRUE)
    loss_mat  <- matrix(loss, nrow = R, ncol = J, byrow = TRUE)
    ratio_mat <- loss_mat / gain_mat
    
    # implied p from FOC under CRRA
    raw <- 1 / (1 + (m - 1) * (ratio_mat ^ r_mat))
    
    invalid <- !is.finite(r_mat) |
      !is.finite(gain_mat) | !is.finite(loss_mat) |
      gain_mat <= 0 | loss_mat <= 0 |
      !is.finite(raw)
    
    p_hat <- pmin(pmax(raw, 0), 1)
    
    # fallback for invalid values
    p_hat[loss_mat <= 0 & !is.finite(raw)] <- 1
    p_hat[gain_mat <= 0 & !is.finite(raw)] <- 0
    p_hat[is.na(p_hat) & !is.na(raw)] <- pmin(pmax(raw[is.na(p_hat) & !is.na(raw)], 0), 1)
    
    flagged <- invalid | raw < 0 | raw > 1
    flagged[is.na(p_hat)] <- TRUE
    
    # --------------------------------------------------------
    # trial-level summaries
    # --------------------------------------------------------
    
    trial_sum <- t(vapply(seq_len(J), function(j) {
      summ_prob(p_hat[, j])
    }, numeric(5L)))
    
    trials_tbl <- data.table(
      dataset = ds,
      treatment = tr,
      pid = as.character(d$pid),
      sequence = as.character(d$seq),
      stake = as.numeric(d$stake),
      n_draws = R,
      n_flagged = colSums(flagged, na.rm = TRUE),
      share_flagged = colMeans(flagged, na.rm = TRUE),
      median = trial_sum[, "median"],
      mean   = trial_sum[, "mean"],
      q025   = trial_sum[, "q025"],
      q975   = trial_sum[, "q975"],
      p_gt_05 = trial_sum[, "p_gt_05"]
    )
    
    # --------------------------------------------------------
    # participant-level summaries
    # --------------------------------------------------------
    
    pid_levels <- unique(as.character(d$pid))
    
    make_pid_tbl <- function(exclude_flagged = FALSE) {
      
      rows <- vector("list", length(pid_levels))
      
      for (i in seq_along(pid_levels)) {
        
        pid_i <- pid_levels[i]
        j_idx <- which(d$pid == pid_i)
        
        vals <- p_hat[, j_idx, drop = FALSE]
        flg  <- flagged[, j_idx, drop = FALSE]
        
        vals_vec <- as.vector(vals)
        flg_vec  <- as.vector(flg)
        
        if (exclude_flagged) {
          vals_use <- vals_vec[!flg_vec]
        } else {
          vals_use <- vals_vec
        }
        
        s <- summ_prob(vals_use)
        
        rows[[i]] <- data.table(
          dataset = ds,
          treatment = tr,
          pid = pid_i,
          n_trials = length(j_idx),
          n_draws_total = length(vals_vec),
          n_flagged = sum(flg_vec, na.rm = TRUE),
          share_flagged = mean(flg_vec, na.rm = TRUE),
          median = s["median"],
          mean   = s["mean"],
          q025   = s["q025"],
          q975   = s["q975"],
          p_gt_05 = s["p_gt_05"]
        )
      }
      
      rbindlist(rows, fill = TRUE)
    }
    
    pincl_tbl <- make_pid_tbl(FALSE)
    pexcl_tbl <- make_pid_tbl(TRUE)
    
    # --------------------------------------------------------
    # treatment comparison row
    # --------------------------------------------------------
    
    p_vec <- as.vector(p_hat)
    s_tr  <- summ_prob(p_vec)
    
    treat_rows[[tr]] <- data.table(
      dataset = ds,
      treatment = tr,
      median = s_tr["median"],
      mean   = s_tr["mean"],
      q025   = s_tr["q025"],
      q975   = s_tr["q975"],
      p_gt_05 = s_tr["p_gt_05"]
    )
    
    # --------------------------------------------------------
    # save
    # --------------------------------------------------------
    
    if (!skip_all) {
      fwrite(trials_tbl, f_trials)
      fwrite(pincl_tbl,  f_pincl)
      fwrite(pexcl_tbl,  f_pexcl)
      
      saveRDS(
        list(
          dataset = ds,
          treatment = tr,
          multiplier = m,
          endowment = endowment,
          trial_data = d,
          trials = trials_tbl,
          participants_incl = pincl_tbl,
          participants_excl = pexcl_tbl,
          p_hat_draws = p_hat,
          flagged = flagged
        ),
        f_rds
      )
      
      msg("Saved: ", f_trials)
      msg("Saved: ", f_pincl)
      msg("Saved: ", f_pexcl)
      msg("Saved: ", f_rds)
    }
    
    outputs[[tr]] <- list(
      trials = trials_tbl,
      participants_incl = pincl_tbl,
      participants_excl = pexcl_tbl
    )
  }
  
  # --------------------------------------------------------
  # treatment comparison table
  # --------------------------------------------------------
  
  treat_tbl <- rbindlist(treat_rows, fill = TRUE)
  
  f_treat <- file.path(out_dir, "ex3_rq2_treatment_comparison.csv")
  fwrite(treat_tbl, f_treat)
  msg("Saved: ", f_treat)
  
  invisible(outputs)
}