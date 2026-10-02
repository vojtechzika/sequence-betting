# ============================================================
# 63_ex2_tables.R
#
# PURPOSE
#   Produces two output tables per treatment from the fitted
#   EX2 Stan model:
#
#   (1) ex2_summary_<tr>.csv  -- pooled and outcome-level
#       posterior summaries of beta_opt, beta_rt, tau, sigma.
#
#   (2) ex2_participants_<tr>.csv  -- participant-level
#       posterior summaries of RQ1--RQ4 outcomes (mu_b, mu_a,
#       mu_c, mu_h), with mu_a suppressed for participants
#       below the minimum-bets threshold.
#
# INPUT
#   path_out/lotr_scored.csv
#   path_src/master_sequences.csv
#   path_mod/ex2_fit_<tr>.rds
#   path_mod/rq{1..4}_fit_sequences_<tr>_full.rds
#   path_mod/rq{1..4}_pid_levels_<tr>_full.rds
#
# OUTPUT
#   path_out/ex2_summary_<tr>.csv
#   path_out/ex2_participants_<tr>.csv
#
# NOTES
#   - Participant outcomes in ex2_participants are summarized
#     from the upstream RQ models, not conditional on EX2
#     covariate availability.  The RQ2 eligibility mask
#     (n_bets < rq2_min_bets) is applied only to mu_a columns;
#     participants remain in the table for other outcomes.
#   - Directories are created by 00_setup.R; no dir.create here.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex2_tables <- function(cfg) {
  
  # ----------------------------------------------------------
  # 0. Configuration
  # ----------------------------------------------------------
  ds           <- as.character(cfg$run$data_folder)
  treatments   <- unique(as.character(cfg$run$treatment))
  rq2_min_bets <- as.integer(cfg$design$rq2$min_bets)
  
  # ----------------------------------------------------------
  # 1. Load master data (stake filter for bet count)
  # ----------------------------------------------------------
  f_master <- file.path(path_src, "master_sequences.csv")
  if (!file.exists(f_master))
    stop("master_sequences.csv not found at: ", f_master)
  master <- fread(f_master, encoding = "UTF-8")
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  
  # ----------------------------------------------------------
  # 2. Helpers
  # ----------------------------------------------------------
  
  # Posterior summary: median, mean, 95% CrI, P(>0)
  summ <- function(x, prob = TRUE) {
    q <- quantile(x, c(0.025, 0.975), names = FALSE)
    c(
      median = median(x),
      mean   = mean(x),
      q025   = q[1L],
      q975   = q[2L],
      p_gt0  = if (isTRUE(prob)) mean(x > 0) else NA_real_
    )
  }
  
  # Participant-level posteriors from a draws x participants matrix.
  # Returns a data.table with columns pid and <prefix>_{median,mean,q025,q975}.
  mat_tbl <- function(pid_vec, mat, prefix) {
    tmp <- t(apply(mat, 2L, summ, prob = FALSE))
    out <- data.table(pid = pid_vec)
    out[, (paste0(prefix, "_median")) := tmp[, "median"]]
    out[, (paste0(prefix, "_mean"))   := tmp[, "mean"]]
    out[, (paste0(prefix, "_q025"))   := tmp[, "q025"]]
    out[, (paste0(prefix, "_q975"))   := tmp[, "q975"]]
    out
  }
  
  # ----------------------------------------------------------
  # 3. Loop over treatments
  # ----------------------------------------------------------
  outputs <- list()
  
  for (tr in treatments) {
    
    msg("EX2 tables: building for treatment '", tr, "'")
    
    f_ex2_fit <- file.path(path_mod, paste0("ex2_fit_", tr, ".rds"))
    
    rq_files <- lapply(1:4, function(k) list(
      fit = file.path(path_mod, paste0("rq", k, "_fit_sequences_", tr, "_full.rds")),
      pid = file.path(path_mod, paste0("rq", k, "_pid_levels_",    tr, "_full.rds"))
    ))
    
    f_summary <- file.path(path_out, paste0("ex2_summary_",      tr, ".csv"))
    f_parts   <- file.path(path_out, paste0("ex2_participants_", tr, ".csv"))
    
    skip_summary <- should_skip(f_summary, cfg, "output",
                                paste0("EX2 summary (", ds, "/", tr, ")"))
    skip_parts   <- should_skip(f_parts,   cfg, "output",
                                paste0("EX2 participants (", ds, "/", tr, ")"))
    
    if (skip_summary && skip_parts) next
    
    # ============================================================
    # 3a. EX2 MODEL SUMMARY
    # ============================================================
    if (!skip_summary) {
      
      if (!file.exists(f_ex2_fit))
        stop("EX2 fit not found for tr='", tr, "': ", f_ex2_fit)
      
      fit  <- readRDS(f_ex2_fit)
      post <- rstan::extract(fit)
      
      req_params <- c("beta_opt_bar", "beta_rt_bar",
                      "tau_opt", "tau_rt",
                      "beta_opt_k", "beta_rt_k", "sigma_k")
      miss <- setdiff(req_params, names(post))
      if (length(miss) > 0L)
        stop("Stan posterior missing parameters: ", paste(miss, collapse = ", "))
      
      rows <- list(
        data.table(block = "pooled", outcome = NA_character_, term = "beta_opt_bar",
                   t(summ(post$beta_opt_bar))),
        data.table(block = "pooled", outcome = NA_character_, term = "beta_rt_bar",
                   t(summ(post$beta_rt_bar))),
        data.table(block = "pooled", outcome = NA_character_, term = "tau_opt",
                   t(summ(post$tau_opt, prob = FALSE))),
        data.table(block = "pooled", outcome = NA_character_, term = "tau_rt",
                   t(summ(post$tau_rt,  prob = FALSE)))
      )
      
      outcome_labels <- c("b", "a", "c", "h")
      for (k in seq_along(outcome_labels)) {
        lbl <- outcome_labels[k]
        rows[[length(rows) + 1L]] <- data.table(
          block = "outcome", outcome = lbl, term = "beta_opt",
          t(summ(post$beta_opt_k[, k]))
        )
        rows[[length(rows) + 1L]] <- data.table(
          block = "outcome", outcome = lbl, term = "beta_rt",
          t(summ(post$beta_rt_k[, k]))
        )
        rows[[length(rows) + 1L]] <- data.table(
          block = "outcome", outcome = lbl, term = "sigma",
          t(summ(post$sigma_k[, k], prob = FALSE))
        )
      }
      
      tbl_summary <- rbindlist(rows, fill = TRUE)
      tbl_summary[, `:=`(dataset = ds, treatment = tr)]
      setcolorder(tbl_summary, c(
        "dataset", "treatment", "block", "outcome", "term",
        "median", "mean", "q025", "q975", "p_gt0"
      ))
      
      fwrite(tbl_summary, f_summary)
      msg("Saved: ", f_summary)
    }
    
    # ============================================================
    # 3b. PARTICIPANT TABLE
    # ============================================================
    if (!skip_parts) {
      
      # Betting trial counts for the RQ2 eligibility mask
      bet_n <- master[treat == tr & stake > 0, .(n_bets = .N), by = pid]
      
      # Parameter names for each RQ
      outcome_params <- c("mu_b_i", "mu_a_i", "mu_c_i", "mu_h_i")
      outcome_prefix <- c("mu_b",   "mu_a",   "mu_c",   "mu_h")
      
      tbl_list <- vector("list", 4L)
      
      for (k in 1:4) {
        f_fit_k <- rq_files[[k]]$fit
        f_pid_k <- rq_files[[k]]$pid
        if (!file.exists(f_fit_k)) stop("Missing upstream fit: ",        f_fit_k)
        if (!file.exists(f_pid_k)) stop("Missing upstream pid levels: ", f_pid_k)
        
        pid_k  <- as.character(readRDS(f_pid_k))
        post_k <- rstan::extract(readRDS(f_fit_k))
        
        param <- outcome_params[k]
        if (is.null(post_k[[param]]))
          stop("Posterior parameter '", param, "' not found in RQ", k, " fit.")
        
        tbl_list[[k]] <- mat_tbl(pid_k, post_k[[param]], outcome_prefix[k])
      }
      
      # Apply RQ2 minimum-bets mask on mu_a columns only
      tbl_list[[2L]] <- merge(tbl_list[[2L]], bet_n, by = "pid", all.x = TRUE)
      tbl_list[[2L]][is.na(n_bets), n_bets := 0L]
      
      mask_cols <- c("mu_a_median", "mu_a_mean", "mu_a_q025", "mu_a_q975")
      tbl_list[[2L]][
        n_bets < rq2_min_bets,
        (mask_cols) := lapply(mask_cols, function(.) NA_real_)
      ]
      
      # Union of participant ids across all four RQ models
      pid_all   <- sort(unique(unlist(lapply(tbl_list, `[[`, "pid"))))
      tbl_parts <- data.table(pid = pid_all)
      
      # Merge RQ1
      tbl_parts <- merge(tbl_parts, tbl_list[[1L]], by = "pid", all.x = TRUE)
      
      # Merge RQ2 (carry n_bets and masked mu_a columns)
      tbl_parts <- merge(
        tbl_parts,
        tbl_list[[2L]][, c("pid", "n_bets",
                           "mu_a_median", "mu_a_mean",
                           "mu_a_q025",   "mu_a_q975")],
        by = "pid", all.x = TRUE
      )
      
      # Merge RQ3 and RQ4
      tbl_parts <- merge(tbl_parts, tbl_list[[3L]], by = "pid", all.x = TRUE)
      tbl_parts <- merge(tbl_parts, tbl_list[[4L]], by = "pid", all.x = TRUE)
      
      tbl_parts[, `:=`(dataset = ds, treatment = tr, rq2_min_bets = rq2_min_bets)]
      
      setcolorder(tbl_parts, c(
        "dataset", "treatment", "pid", "rq2_min_bets", "n_bets",
        "mu_b_median", "mu_b_mean", "mu_b_q025", "mu_b_q975",
        "mu_a_median", "mu_a_mean", "mu_a_q025", "mu_a_q975",
        "mu_c_median", "mu_c_mean", "mu_c_q025", "mu_c_q975",
        "mu_h_median", "mu_h_mean", "mu_h_q025", "mu_h_q975"
      ))
      
      setorder(tbl_parts, pid)
      fwrite(tbl_parts, f_parts)
      msg("Saved: ", f_parts)
      
      outputs[[tr]] <- list(summary = f_summary, participants = f_parts)
    }
  }
  
  invisible(outputs)
}