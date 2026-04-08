# ============================================================
# scripts/analysis/63_ex2_tables.R
#
# EX2 tables
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

ex2_tables <- function(cfg) {
  
  ds <- as.character(cfg$run$dataset)
  treatments <- as.character(cfg$run$treatment)
  rq2_min_bets <- as.integer(cfg$design$rq2$min_bets)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  master <- fread(f_master)
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  
  # ------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------
  
  summ <- function(x, prob = TRUE) {
    q <- quantile(x, c(.025, .975), names = FALSE)
    c(
      median = median(x),
      mean   = mean(x),
      q025   = q[1],
      q975   = q[2],
      p_gt0  = if (prob) mean(x > 0) else NA_real_
    )
  }
  
  mat_tbl <- function(pid, mat, prefix) {
    
    tmp <- t(apply(mat, 2, summ, prob = FALSE))
    
    out <- data.table(pid = pid)
    
    out[, paste0(prefix, "_median") := tmp[, "median"]]
    out[, paste0(prefix, "_mean")   := tmp[, "mean"]]
    out[, paste0(prefix, "_q025")   := tmp[, "q025"]]
    out[, paste0(prefix, "_q975")   := tmp[, "q975"]]
    
    out
  }
  
  outputs <- list()
  
  for (tr in treatments) {
    
    msg("Building EX2 tables for treatment: ", tr)
    
    f_ex2_fit <- file.path(mod_dir, paste0("ex2_fit_", tr, ".rds"))
    
    f_rq1_fit <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, "_full.rds"))
    f_rq1_pid <- file.path(mod_dir, paste0("rq1_pid_levels_", tr, "_full.rds"))
    
    f_rq2_fit <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, "_full.rds"))
    f_rq2_pid <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, "_full.rds"))
    
    f_rq3_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, "_full.rds"))
    f_rq3_pid <- file.path(mod_dir, paste0("rq3_pid_levels_", tr, "_full.rds"))
    
    f_rq4_fit <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_rq4_pid <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, "_full.rds"))
    
    f_summary <- file.path(out_dir, paste0("ex2_summary_", tr, ".csv"))
    f_parts   <- file.path(out_dir, paste0("ex2_participants_", tr, ".csv"))
    
    skip_summary <- should_skip(
      paths = f_summary,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX2 summary (", ds, "/", tr, ")")
    )
    
    skip_parts <- should_skip(
      paths = f_parts,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX2 participants (", ds, "/", tr, ")")
    )
    
    if (skip_summary && skip_parts) next
    
    # ============================================================
    # 1) EX2 SUMMARY
    # ============================================================
    
    if (!skip_summary) {
      
      fit  <- readRDS(f_ex2_fit)
      post <- rstan::extract(fit)
      
      rows <- list(
        
        data.table(
          block="pooled", outcome=NA, term="beta_opt_bar",
          t(summ(post$beta_opt_bar))
        ),
        
        data.table(
          block="pooled", outcome=NA, term="beta_rt_bar",
          t(summ(post$beta_rt_bar))
        ),
        
        data.table(
          block="pooled", outcome=NA, term="tau_opt",
          t(summ(post$tau_opt, FALSE))
        ),
        
        data.table(
          block="pooled", outcome=NA, term="tau_rt",
          t(summ(post$tau_rt, FALSE))
        )
      )
      
      outcomes <- c("b","a","c","h")
      
      for (k in seq_along(outcomes)) {
        
        rows[[length(rows)+1]] <- data.table(
          block="outcome",
          outcome=outcomes[k],
          term="beta_opt",
          t(summ(post$beta_opt_k[,k]))
        )
        
        rows[[length(rows)+1]] <- data.table(
          block="outcome",
          outcome=outcomes[k],
          term="beta_rt",
          t(summ(post$beta_rt_k[,k]))
        )
        
        rows[[length(rows)+1]] <- data.table(
          block="outcome",
          outcome=outcomes[k],
          term="sigma",
          t(summ(post$sigma_k[,k], FALSE))
        )
      }
      
      tbl_summary <- rbindlist(rows, fill = TRUE)
      
      tbl_summary[,`:=`(
        dataset=ds,
        treatment=tr
      )]
      
      setcolorder(tbl_summary,c(
        "dataset","treatment","block","outcome","term",
        "median","mean","q025","q975","p_gt0"
      ))
      
      fwrite(tbl_summary,f_summary)
      
      msg("Saved EX2 summary table: ", f_summary)
    }
    
    # ============================================================
    # 2) PARTICIPANT TABLE
    # ============================================================
    
    if (!skip_parts) {
      
      betN <- master[
        treat==tr & stake>0,
        .(n_bets=.N),
        by=pid
      ]
      
      post1 <- rstan::extract(readRDS(f_rq1_fit))
      post2 <- rstan::extract(readRDS(f_rq2_fit))
      post3 <- rstan::extract(readRDS(f_rq3_fit))
      post4 <- rstan::extract(readRDS(f_rq4_fit))
      
      pid1 <- as.character(readRDS(f_rq1_pid))
      pid2 <- as.character(readRDS(f_rq2_pid))
      pid3 <- as.character(readRDS(f_rq3_pid))
      pid4 <- as.character(readRDS(f_rq4_pid))
      
      rq1_tbl <- mat_tbl(pid1, post1$mu_b_i, "mu_b")
      rq2_tbl <- mat_tbl(pid2, post2$mu_a_i, "mu_a")
      rq3_tbl <- mat_tbl(pid3, post3$mu_c_i, "mu_c")
      rq4_tbl <- mat_tbl(pid4, post4$mu_h_i, "mu_h")
      
      rq2_tbl <- merge(rq2_tbl, betN, by="pid", all.x=TRUE)
      rq2_tbl[is.na(n_bets), n_bets := 0]
      
      rq2_tbl[n_bets < rq2_min_bets,
              c("mu_a_median","mu_a_mean","mu_a_q025","mu_a_q975") :=
                .(NA_real_,NA_real_,NA_real_,NA_real_)]
      
      pid_all <- sort(unique(c(pid1,pid2,pid3,pid4)))
      
      tbl_parts <- data.table(pid=pid_all)
      
      tbl_parts <- merge(tbl_parts,rq1_tbl,by="pid",all.x=TRUE)
      
      tbl_parts <- merge(
        tbl_parts,
        rq2_tbl[,c("pid","n_bets",
                   "mu_a_median","mu_a_mean",
                   "mu_a_q025","mu_a_q975")],
        by="pid",
        all.x=TRUE
      )
      
      tbl_parts <- merge(tbl_parts,rq3_tbl,by="pid",all.x=TRUE)
      tbl_parts <- merge(tbl_parts,rq4_tbl,by="pid",all.x=TRUE)
      
      tbl_parts[,`:=`(
        dataset=ds,
        treatment=tr,
        rq2_min_bets=rq2_min_bets
      )]
      
      setcolorder(tbl_parts,c(
        "dataset","treatment","pid","rq2_min_bets","n_bets",
        "mu_b_median","mu_b_mean","mu_b_q025","mu_b_q975",
        "mu_a_median","mu_a_mean","mu_a_q025","mu_a_q975",
        "mu_c_median","mu_c_mean","mu_c_q025","mu_c_q975",
        "mu_h_median","mu_h_mean","mu_h_q025","mu_h_q975"
      ))
      
      setorder(tbl_parts,pid)
      
      fwrite(tbl_parts,f_parts)
      
      msg("Saved EX2 participant table: ", f_parts)
    }
    
  }
  
  invisible(outputs)
}