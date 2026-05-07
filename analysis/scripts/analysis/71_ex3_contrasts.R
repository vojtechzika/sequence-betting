# ============================================================
# scripts/analysis/71_ex3_rq_contrasts.R
#
# PURPOSE
#   Computes FN vs FP posterior contrasts for EX3 (RQ1–RQ3).
#   Sequence-level and group-level participant contrasts are
#   computed draw-wise from the fitted Stan objects.
#
# SAMPLE
#   Full sample within each treatment (no HL-consistency or
#   normative-better filtering). EX3 is exploratory throughout.
#
# INPUT
#   path_mod/rq{1,2,3}_fit_sequences_<fn>_full.rds
#   path_mod/rq{1,2,3}_pid_levels_<fn>_full.rds
#   path_mod/rq{1,2,3}_seq_levels_<fn>_full.rds
#   (analogously for fp treatment)
#
# OUTPUT
#   path_out/ex3_rq{1,2,3}_sequence_contrasts.csv
#   path_out/ex3_rq{1,2,3}_participant_contrasts.csv
#   path_out/ex3_rq{1,2,3}_spearman.csv
#
# NOTES
#   - Treatments, thresholds, and sample tag are read from cfg;
#     no values are hardcoded.
#   - FN = cfg$run$treatment[1] (confirmatory); FP = [2].
#   - RQ3 reports sensitivity thresholds at rho_Delta in {0.03, 0.05}
#     per prereg; these are drawn from cfg$design$rq3$rho[1:2].
#   - RQ4 contrasts are omitted (selection confound across frames;
#     see prereg EX3 §RQ4).
#   - Participant-level contrasts are group-level (between-subject):
#     draw-wise treatment means are differenced.
#
# CALL ORDER IN PIPELINE
#   rq1_stan(cfg)      -- must be run first for all treatments
#   rq2_stan(cfg)      -- must be run first for all treatments
#   rq3_stan(cfg)      -- must be run first for all treatments
#   ex3_contrasts(cfg) -- this script
# ============================================================

ex3_contrasts <- function(cfg) {
  
  # ----------------------------------------------------------
  # Resolve treatments and sample tag from cfg
  # ----------------------------------------------------------
  
  stopifnot(length(cfg$run$treatment) >= 2L)
  fn  <- cfg$run$treatment[1L]   # confirmatory (FN)
  fp  <- cfg$run$treatment[2L]   # exploratory  (FP)
  tag <- "full"                  # EX3 always uses the full (unfiltered) sample
  
  # ----------------------------------------------------------
  # Resolve per-RQ thresholds from cfg$design
  # ----------------------------------------------------------
  # cfg$design$rq{k}$rho: first element = primary threshold,
  # remaining = sensitivity. For RQ3 the prereg specifies
  # rho_Delta in {0.03, 0.05}, which correspond to positions [2,1].
  
  rq_map <- list(
    list(
      rq         = "rq1",
      mu_seq     = "mu_b",
      mu_pid     = "mu_b_i",
      thresholds = cfg$design$rq1$rho[1L]
    ),
    list(
      rq         = "rq2",
      mu_seq     = "mu_a",
      mu_pid     = "mu_a_i",
      thresholds = cfg$design$rq2$rho[1L]
    ),
    list(
      rq         = "rq3",
      mu_seq     = "mu_c",
      mu_pid     = "mu_c_i",
      thresholds = sort(unique(cfg$design$rq3$rho[1:2]))  # {0.03, 0.05}
    )
    # RQ4: FN–FP contrasts not computed; selection into conditioning
    # set (b_is = 1) differs sharply across frames (prereg §EX3 RQ4).
  )
  
  # ----------------------------------------------------------
  # Helpers
  # ----------------------------------------------------------
  
  summ_delta <- function(x, thresholds) {
    q    <- quantile(x, c(0.025, 0.975), names = FALSE)
    base <- c(
      median = median(x),
      mean   = mean(x),
      q025   = q[1L],
      q975   = q[2L],
      p_gt0  = mean(x > 0)
    )
    thresh_stats        <- vapply(thresholds, function(thr) mean(abs(x) > thr), numeric(1L))
    names(thresh_stats) <- paste0("p_abs_gt_", gsub("\\.", "", sprintf("%05.3f", thresholds)))
    c(base, thresh_stats)
  }
  
  summ_spearman <- function(rho_draws) {
    q <- quantile(rho_draws, c(0.025, 0.975), names = FALSE)
    c(
      median = median(rho_draws),
      mean   = mean(rho_draws),
      q025   = q[1L],
      q975   = q[2L],
      p_gt0  = mean(rho_draws > 0)
    )
  }
  
  # ----------------------------------------------------------
  # Main loop
  # ----------------------------------------------------------
  
  for (rq in rq_map) {
    
    rq_name    <- rq$rq
    mu_seq_nm  <- rq$mu_seq
    mu_pid_nm  <- rq$mu_pid
    thresholds <- rq$thresholds
    
    msg("EX3 contrasts: ", toupper(rq_name), " (", fn, " vs ", fp, ")")
    
    # ---- file paths -------------------------------------------------
    
    make_mod_path <- function(stem, tr)
      file.path(path_mod, paste0(rq_name, "_", stem, "_", tr, "_", tag, ".rds"))
    
    f_fn_fit <- make_mod_path("fit_sequences", fn)
    f_fp_fit <- make_mod_path("fit_sequences", fp)
    f_fn_pid <- make_mod_path("pid_levels",    fn)
    f_fp_pid <- make_mod_path("pid_levels",    fp)
    f_fn_seq <- make_mod_path("seq_levels",    fn)
    f_fp_seq <- make_mod_path("seq_levels",    fp)
    
    stopifnot(
      file.exists(f_fn_fit), file.exists(f_fp_fit),
      file.exists(f_fn_pid), file.exists(f_fp_pid),
      file.exists(f_fn_seq), file.exists(f_fp_seq)
    )
    
    f_seq_out <- file.path(path_out, paste0("ex3_", rq_name, "_sequence_contrasts.csv"))
    f_pid_out <- file.path(path_out, paste0("ex3_", rq_name, "_participant_contrasts.csv"))
    f_spm_out <- file.path(path_out, paste0("ex3_", rq_name, "_spearman.csv"))
    
    skip_seq <- should_skip(f_seq_out, cfg, "output",
                            paste0("EX3 ", toupper(rq_name), " sequence contrasts"))
    skip_pid <- should_skip(f_pid_out, cfg, "output",
                            paste0("EX3 ", toupper(rq_name), " participant contrasts"))
    skip_spm <- should_skip(f_spm_out, cfg, "output",
                            paste0("EX3 ", toupper(rq_name), " Spearman correlation"))
    
    if (skip_seq && skip_pid && skip_spm) {
      msg("  All outputs exist and overwrite_outputs = FALSE; skipping ", toupper(rq_name))
      next
    }
    
    # ---- load posteriors (shared across all three outputs) ----------
    
    post_fn   <- rstan::extract(readRDS(f_fn_fit))
    post_fp   <- rstan::extract(readRDS(f_fp_fit))
    
    mu_seq_fn <- post_fn[[mu_seq_nm]]   # draws x sequences
    mu_seq_fp <- post_fp[[mu_seq_nm]]
    
    seq_fn <- as.character(readRDS(f_fn_seq))
    seq_fp <- as.character(readRDS(f_fp_seq))
    stopifnot(identical(seq_fn, seq_fp))
    
    # ================================================================
    # (1) Sequence contrasts
    # ================================================================
    
    if (!skip_seq) {
      
      msg("  Computing sequence contrasts")
      
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
        grep("^p_abs_gt_", names(tbl_seq), value = TRUE)
      ))
      
      fwrite(tbl_seq, f_seq_out)
      msg("  Saved: ", f_seq_out)
    }
    
    # ================================================================
    # (2) Participant-level group contrasts (between-subject)
    # ================================================================
    
    if (!skip_pid) {
      
      msg("  Computing participant-level group contrasts")
      
      mu_pid_fn <- post_fn[[mu_pid_nm]]   # draws x participants (FN)
      mu_pid_fp <- post_fp[[mu_pid_nm]]   # draws x participants (FP)
      
      stopifnot(is.matrix(mu_pid_fn), is.matrix(mu_pid_fp))
      
      pid_fn <- as.character(readRDS(f_fn_pid))
      pid_fp <- as.character(readRDS(f_fp_pid))
      
      fn_draw    <- rowMeans(mu_pid_fn, na.rm = TRUE)
      fp_draw    <- rowMeans(mu_pid_fp, na.rm = TRUE)
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
      msg("  Saved: ", f_pid_out)
    }
    
    # ================================================================
    # (3) Spearman correlation between FN and FP sequence-level vectors
    # ================================================================
    
    if (!skip_spm) {
      
      msg("  Computing Spearman correlation (FN vs FP sequence vectors)")
      
      Tn        <- min(nrow(mu_seq_fn), nrow(mu_seq_fp))
      rho_draws <- vapply(seq_len(Tn), function(t)
        cor(mu_seq_fn[t, ], mu_seq_fp[t, ], method = "spearman"),
        numeric(1L)
      )
      
      tbl_spm <- data.table(
        rq           = rq_name,
        treatment_fn = fn,
        treatment_fp = fp,
        t(summ_spearman(rho_draws))
      )
      
      fwrite(tbl_spm, f_spm_out)
      msg("  Saved: ", f_spm_out)
    }
    
  } # end rq loop
  
  invisible(TRUE)
}