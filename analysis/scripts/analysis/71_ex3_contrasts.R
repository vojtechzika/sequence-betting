# ============================================================
# scripts/analysis/71_ex3_rq_contrasts.R
#
# EX3 – FN vs FP contrasts for RQ1–RQ4
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex3_contrasts <- function(cfg) {
  
  ds <- as.character(cfg$run$dataset)
  
  fn <- "m25"
  fp <- "m19"
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------------------------------------------------
  # helper
  # ------------------------------------------------------------
  
  summ_delta <- function(x) {
    q <- quantile(x, c(.025, .975), names = FALSE)
    c(
      median = median(x),
      mean   = mean(x),
      q025   = q[1],
      q975   = q[2],
      p_gt0  = mean(x > 0),
      p_abs_gt_005 = mean(abs(x) > 0.05)
    )
  }
  
  # ------------------------------------------------------------
  # RQ map
  # ------------------------------------------------------------
  
  rq_map <- list(
    list(rq = "rq1", mu_seq = "mu_b", mu_pid = "mu_b_i"),
    list(rq = "rq2", mu_seq = "mu_a", mu_pid = "mu_a_i"),
    list(rq = "rq3", mu_seq = "mu_c", mu_pid = "mu_c_i"),
    list(rq = "rq4", mu_seq = "mu_h", mu_pid = "mu_h_i")
  )
  
  outputs <- list()
  
  for (rq in rq_map) {
    
    rq_name <- rq$rq
    mu_seq_name <- rq$mu_seq
    mu_pid_name <- rq$mu_pid
    
    msg("Running EX3 contrasts for ", toupper(rq_name), " (FN vs FP)")
    
    # ----------------------------------------------------------
    # inputs
    # ----------------------------------------------------------
    
    f_fn_fit <- file.path(mod_dir, paste0(rq_name, "_fit_sequences_", fn, "_full.rds"))
    f_fp_fit <- file.path(mod_dir, paste0(rq_name, "_fit_sequences_", fp, "_full.rds"))
    
    f_fn_pid <- file.path(mod_dir, paste0(rq_name, "_pid_levels_", fn, "_full.rds"))
    f_fp_pid <- file.path(mod_dir, paste0(rq_name, "_pid_levels_", fp, "_full.rds"))
    
    f_fn_seq <- file.path(mod_dir, paste0(rq_name, "_seq_levels_", fn, "_full.rds"))
    f_fp_seq <- file.path(mod_dir, paste0(rq_name, "_seq_levels_", fp, "_full.rds"))
    
    # ----------------------------------------------------------
    # outputs
    # ----------------------------------------------------------
    
    f_seq_out <- file.path(out_dir, paste0("ex3_", rq_name, "_sequence_contrasts.csv"))
    f_pid_out <- file.path(out_dir, paste0("ex3_", rq_name, "_participant_contrasts.csv"))
    
    skip_seq <- should_skip(
      paths = f_seq_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 ", toupper(rq_name), " sequence contrasts (", ds, ")")
    )
    
    skip_pid <- should_skip(
      paths = f_pid_out,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX3 ", toupper(rq_name), " participant contrasts (", ds, ")")
    )
    
    if (skip_seq && skip_pid) next
    
    # ----------------------------------------------------------
    # load fits
    # ----------------------------------------------------------
    
    fit_fn <- readRDS(f_fn_fit)
    fit_fp <- readRDS(f_fp_fit)
    
    post_fn <- rstan::extract(fit_fn)
    post_fp <- rstan::extract(fit_fp)
    
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
      
      S <- ncol(mu_seq_fn)
      rows <- vector("list", S)
      
      for (s in seq_len(S)) {
        
        delta_draw <- mu_seq_fn[, s] - mu_seq_fp[, s]
        
        rows[[s]] <- data.table(
          sequence = seq_fn[s],
          seq_id   = s,
          t(summ_delta(delta_draw))
        )
      }
      
      tbl_seq <- rbindlist(rows)
      
      tbl_seq[, `:=`(
        dataset = ds,
        rq = rq_name,
        treatment_fn = fn,
        treatment_fp = fp
      )]
      
      setcolorder(
        tbl_seq,
        c(
          "dataset",
          "rq",
          "sequence",
          "seq_id",
          "treatment_fn",
          "treatment_fp",
          "median",
          "mean",
          "q025",
          "q975",
          "p_gt0",
          "p_abs_gt_005"
        )
      )
      
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
      
      Tn <- min(length(fn_draw), length(fp_draw))
      delta_draw <- fn_draw[seq_len(Tn)] - fp_draw[seq_len(Tn)]
      
      q <- quantile(delta_draw, c(.025, .975), names = FALSE)
      
      tbl_pid <- data.table(
        dataset = ds,
        rq = rq_name,
        treatment_fn = fn,
        treatment_fp = fp,
        n_pid_fn = length(pid_fn),
        n_pid_fp = length(pid_fp),
        median = median(delta_draw),
        mean   = mean(delta_draw),
        q025   = q[1],
        q975   = q[2],
        p_gt0  = mean(delta_draw > 0),
        p_abs_gt_005 = mean(abs(delta_draw) > 0.05)
      )
      
      fwrite(tbl_pid, f_pid_out)
      msg("Saved: ", f_pid_out)
    }
    
    outputs[[rq_name]] <- list(
      sequence_file = if (!skip_seq) f_seq_out else NULL,
      participant_file = if (!skip_pid) f_pid_out else NULL
    )
  }
  
  invisible(outputs)
}