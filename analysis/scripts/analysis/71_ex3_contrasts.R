# ============================================================
# scripts/analysis/71_ex3_rq_contrasts.R
#
# PURPOSE
#   Computes FN vs FP posterior contrasts for EX3 (RQ1-RQ3).
#   Sequence-level and group-level participant contrasts are
#   computed draw-wise from the fitted Stan objects.
#
# SAMPLE:
#   Full sample within each treatment (no HL-consistency or
#   normative-better filtering). EX3 is exploratory throughout.
#
# INPUT
#   path_mod/rq1_fit_sequences_<tr>_full.rds
#   path_mod/rq1_pid_levels_<tr>_full.rds
#   path_mod/rq1_seq_levels_<tr>_full.rds
#   (analogously for rq2, rq3)
#
# OUTPUT
#   path_out/ex3_rq1_sequence_contrasts.csv
#   path_out/ex3_rq1_participant_contrasts.csv
#   path_out/ex3_rq1_spearman.csv
#   (analogously for rq2, rq3)
#
# NOTES
#   - RQ3 Spearman correlation between FN and FP vectors {mu_s^c}
#     is preregistered. Added for RQ1 and RQ2 for completeness.
#   - RQ3 reports sensitivity thresholds at both 0.03 and 0.05
#     per prereg (rho_Delta in {0.03, 0.05}).
#   - RQ4 contrasts are commented out; to be added later.
#   - Participant-level contrasts are group-level (between-subject).
#
# CALL ORDER IN PIPELINE:
#   rq1_stan(cfg)         -- must be run first for all treatments
#   rq2_stan(cfg)         -- must be run first for all treatments
#   rq3_stan(cfg)         -- must be run first for all treatments
#   ex3_contrasts(cfg)    -- this script
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex3_contrasts <- function(cfg) {
  
  fn <- "m25"
  fp <- "m19"
  
  # ------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------
  
  summ_delta <- function(x, thresholds = 0.05) {
    q <- quantile(x, c(.025, .975), names = FALSE)
    base <- c(
      median = median(x),
      mean   = mean(x),
      q025   = q[1],
      q975   = q[2],
      p_gt0  = mean(x > 0)
    )
    thresh_stats <- sapply(thresholds, function(thr) mean(abs(x) > thr))
    names(thresh_stats) <- paste0(
      "p_abs_gt_", gsub("\\.", "", sprintf("%05.3f", thresholds))
    )
    c(base, thresh_stats)
  }
  
  summ_spearman <- function(rho_draws) {
    q <- quantile(rho_draws, c(.025, .975), names = FALSE)
    c(
      median   = median(rho_draws),
      mean     = mean(rho_draws),
      q025     = q[1],
      q975     = q[2],
      p_gt0    = mean(rho_draws > 0)
    )
  }
  
  # ------------------------------------------------------------
  # RQ map
  # ------------------------------------------------------------
  
  rq_map <- list(
    list(rq = "rq1", mu_seq = "mu_b", mu_pid = "mu_b_i", thresholds = 0.05),
    list(rq = "rq2", mu_seq = "mu_a", mu_pid = "mu_a_i", thresholds = 0.05),
    list(rq = "rq3", mu_seq = "mu_c", mu_pid = "mu_c_i", thresholds = c(0.03, 0.05))
    # list(rq = "rq4", mu_seq = "mu_h", mu_pid = "mu_h_i", thresholds = 0.05)
    # RQ4 to be added later.
  )
  
  for (rq in rq_map) {
    
    rq_name     <- rq$rq
    mu_seq_name <- rq$mu_seq
    mu_pid_name <- rq$mu_pid
    thresholds  <- rq$thresholds
    
    msg("Running EX3 contrasts for ", toupper(rq_name), " (FN vs FP)")
    
    # ----------------------------------------------------------
    # paths
    # ----------------------------------------------------------
    
    f_fn_fit <- file.path(path_mod, paste0(rq_name, "_fit_sequences_", fn, "_full.rds"))
    f_fp_fit <- file.path(path_mod, paste0(rq_name, "_fit_sequences_", fp, "_full.rds"))
    
    f_fn_pid <- file.path(path_mod, paste0(rq_name, "_pid_levels_", fn, "_full.rds"))
    f_fp_pid <- file.path(path_mod, paste0(rq_name, "_pid_levels_", fp, "_full.rds"))
    
    f_fn_seq <- file.path(path_mod, paste0(rq_name, "_seq_levels_", fn, "_full.rds"))
    f_fp_seq <- file.path(path_mod, paste0(rq_name, "_seq_levels_", fp, "_full.rds"))
    
    f_seq_out  <- file.path(path_out, paste0("ex3_", rq_name, "_sequence_contrasts.csv"))
    f_pid_out  <- file.path(path_out, paste0("ex3_", rq_name, "_participant_contrasts.csv"))
    f_spm_out  <- file.path(path_out, paste0("ex3_", rq_name, "_spearman.csv"))
    
    stopifnot(
      file.exists(f_fn_fit), file.exists(f_fp_fit),
      file.exists(f_fn_pid), file.exists(f_fp_pid),
      file.exists(f_fn_seq), file.exists(f_fp_seq)
    )
    
    skip_seq <- should_skip(
      paths = f_seq_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 ", toupper(rq_name), " sequence contrasts")
    )
    
    skip_pid <- should_skip(
      paths = f_pid_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 ", toupper(rq_name), " participant contrasts")
    )
    
    skip_spm <- should_skip(
      paths = f_spm_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 ", toupper(rq_name), " Spearman correlation")
    )
    
    if (skip_seq && skip_pid && skip_spm) next
    
    # ----------------------------------------------------------
    # load fits
    # ----------------------------------------------------------
    
    post_fn <- rstan::extract(readRDS(f_fn_fit))
    post_fp <- rstan::extract(readRDS(f_fp_fit))
    
    mu_seq_fn <- post_fn[[mu_seq_name]]
    mu_seq_fp <- post_fp[[mu_seq_name]]
    
    mu_pid_fn <- post_fn[[mu_pid_name]]
    mu_pid_fp <- post_fp[[mu_pid_name]]
    
    pid_fn <- as.character(readRDS(f_fn_pid))
    pid_fp <- as.character(readRDS(f_fp_pid))
    
    seq_fn <- as.character(readRDS(f_fn_seq))
    seq_fp <- as.character(readRDS(f_fp_seq))
    
    stopifnot(identical(seq_fn, seq_fp))
    
    # ==========================================================
    # sequence contrasts
    # ==========================================================
    
    if (!skip_seq) {
      
      msg("Computing ", toupper(rq_name), " sequence contrasts")
      
      S    <- ncol(mu_seq_fn)
      rows <- vector("list", S)
      
      for (s in seq_len(S)) {
        delta_draw <- mu_seq_fn[, s] - mu_seq_fp[, s]
        rows[[s]]  <- data.table(
          sequence = seq_fn[s],
          seq_id   = s,
          t(summ_delta(delta_draw, thresholds))
        )
      }
      
      tbl_seq <- rbindlist(rows)
      tbl_seq[, `:=`(rq = rq_name, treatment_fn = fn, treatment_fp = fp)]
      
      setcolorder(tbl_seq, c(
        "rq", "sequence", "seq_id", "treatment_fn", "treatment_fp",
        "median", "mean", "q025", "q975", "p_gt0",
        names(tbl_seq)[grepl("^p_abs_gt_", names(tbl_seq))]
      ))
      
      fwrite(tbl_seq, f_seq_out)
      msg("Saved: ", f_seq_out)
    }
    
    # ==========================================================
    # participant contrasts (group-level, draw-wise)
    # ==========================================================
    
    if (!skip_pid) {
      
      msg("Computing ", toupper(rq_name), " participant-level group contrasts")
      
      stopifnot(is.matrix(mu_pid_fn), is.matrix(mu_pid_fp))
      
      fn_draw <- rowMeans(mu_pid_fn, na.rm = TRUE)
      fp_draw <- rowMeans(mu_pid_fp, na.rm = TRUE)
      
      Tn         <- min(length(fn_draw), length(fp_draw))
      delta_draw <- fn_draw[seq_len(Tn)] - fp_draw[seq_len(Tn)]
      
      tbl_pid <- data.table(
        rq           = rq_name,
        treatment_fn = fn,
        treatment_fp = fp,
        n_pid_fn     = length(pid_fn),
        n_pid_fp     = length(pid_fp),
        t(summ_delta(delta_draw, thresholds))
      )
      
      fwrite(tbl_pid, f_pid_out)
      msg("Saved: ", f_pid_out)
    }
    
    # ==========================================================
    # Spearman correlation between FN and FP sequence-level vectors
    # ==========================================================
    
    if (!skip_spm) {
      
      msg("Computing ", toupper(rq_name), " Spearman correlation (FN vs FP)")
      
      Tn <- min(nrow(mu_seq_fn), nrow(mu_seq_fp))
      
      rho_draws <- vapply(seq_len(Tn), function(t) {
        cor(mu_seq_fn[t, ], mu_seq_fp[t, ], method = "spearman")
      }, numeric(1L))
      
      tbl_spm <- data.table(
        rq           = rq_name,
        treatment_fn = fn,
        treatment_fp = fp,
        t(summ_spearman(rho_draws))
      )
      
      fwrite(tbl_spm, f_spm_out)
      msg("Saved: ", f_spm_out)
    }
  }
  
  invisible(TRUE)
}