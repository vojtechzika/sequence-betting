# ============================================================
# scripts/analysis/32_rq3_tables.R
# RQ3 tables (exploratory + confirmatory) from *_full / *_conf fits
# Naming: rq3_<tr>_<type>_<unit>.csv
#   type = exploratory (tag full) | confirmatory (tag conf)
#   unit = sequences | participants
# Uses generated quantities:
#   mu_c   (iters x S): per-sequence population mean welfare loss share
#   mu_c_i (iters x N): per-participant mean welfare loss share
# ============================================================

library(data.table)
library(rstan)

rq3_tables <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  stopifnot(!is.null(design$rq3), !is.null(design$rq3$rho))
  
  rho_vec <- as.numeric(design$rq3$rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  
  rho_main <- rho_vec[1]
  rho_nm_main <- gsub("\\.", "", sprintf("%.2f", rho_main))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  tags <- c("full", "conf")
  outputs <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      # If you only produce *_conf for betting_normative==TRUE, skip early:
      if (tag == "conf" && !isTRUE(design$a_flags$betting_normative[[tr]])) next
      
      f_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, "_", tag, ".rds"))
      f_pid <- file.path(mod_dir, paste0("rq3_pid_levels_",  tr, "_", tag, ".rds"))
      f_seq <- file.path(mod_dir, paste0("rq3_seq_levels_",  tr, "_", tag, ".rds"))
      
      if (!file.exists(f_fit) || !file.exists(f_pid) || !file.exists(f_seq)) next
      
      pid_levels <- as.character(readRDS(f_pid))
      seq_levels <- as.character(readRDS(f_seq))
      stopifnot(length(pid_levels) >= 1L, length(seq_levels) >= 1L)
      
      # Tag -> output type label
      type_nm <- if (tag == "full") "exploratory" else "confirmatory"
      
      # Counts from master, restricted to pid_levels/seq_levels actually used in fit
      infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
      stopifnot(file.exists(infile))
      dt <- fread(infile, encoding = "UTF-8")
      stopifnot(all(c("pid","treat","seq","stake") %in% names(dt)))
      
      dt[, pid := as.character(pid)]
      dt[, treat := as.character(treat)]
      dt[, seq := as.character(seq)]
      
      d_counts <- dt[treat == tr & pid %in% pid_levels & seq %in% seq_levels]
      n_trials_by_seq <- d_counts[, .(n_trials = .N), by = seq]
      setkey(n_trials_by_seq, seq)
      
      # Extract posterior
      fit  <- readRDS(f_fit)
      post <- rstan::extract(fit)
      
      if (is.null(post$mu_c))   stop("RQ3 fit missing generated quantity mu_c: ", f_fit)
      if (is.null(post$mu_c_i)) stop("RQ3 fit missing generated quantity mu_c_i: ", f_fit)
      
      mu_s_draws <- post$mu_c     # iters x S
      mu_i_draws <- post$mu_c_i   # iters x N
      
      stopifnot(is.matrix(mu_s_draws), ncol(mu_s_draws) == length(seq_levels))
      stopifnot(is.matrix(mu_i_draws), ncol(mu_i_draws) == length(pid_levels))
      
      # ============================================================
      # SEQUENCES
      # L_s(rho) = P(mu_s^c > rho)
      # labels: strong/moderate/weak/neutral using rho_main
      # ============================================================
      seq_tbl <- data.table(
        sequence    = seq_levels,
        n_trials    = as.integer(n_trials_by_seq[.(seq_levels), n_trials]),
        mu_c_median = apply(mu_s_draws, 2, median),
        mu_c_mean   = apply(mu_s_draws, 2, mean),
        mu_c_q025   = apply(mu_s_draws, 2, quantile, probs = 0.025),
        mu_c_q975   = apply(mu_s_draws, 2, quantile, probs = 0.975)
      )
      seq_tbl[is.na(n_trials), n_trials := 0L]
      
      for (rho in rho_vec) {
        nm <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        seq_tbl[, (nm) := apply(mu_s_draws, 2, function(x) mean(x > rho))]
      }
      
      colL_main <- paste0("L_rho_", rho_nm_main)
      seq_tbl[, loss_label :=
                fifelse(get(colL_main) >= 0.95, "strong",
                        fifelse(get(colL_main) >= 0.80, "moderate",
                                fifelse(get(colL_main) >= 0.50, "weak", "neutral")))]
      
      setorder(seq_tbl, sequence)
      
      f_seq_csv <- file.path(out_dir, paste0("rq3_", tr, "_", type_nm, "_sequences.csv"))
      if (!should_skip(
        paths = f_seq_csv,
        cfg   = cfg,
        type  = "output",
        label = paste0("RQ3 sequences ", type_nm, " (", ds, "/", tr, ")")
      )) {
        fwrite(seq_tbl, f_seq_csv)
        msg("Saved: ", f_seq_csv)
      }
      
      # ============================================================
      # PARTICIPANTS
      # L_i(rho) = P(mu_i^c > rho)
      # labels: solid/likely/leaning/neutral using rho_main
      # ============================================================
      part_tbl <- data.table(
        pid         = pid_levels,
        mu_c_median = apply(mu_i_draws, 2, median),
        mu_c_mean   = apply(mu_i_draws, 2, mean),
        mu_c_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
        mu_c_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975)
      )
      
      for (rho in rho_vec) {
        nm <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
        part_tbl[, (nm) := apply(mu_i_draws, 2, function(x) mean(x > rho))]
      }
      
      colLi_main <- paste0("L_rho_", rho_nm_main)
      part_tbl[, loss_label :=
                 fifelse(get(colLi_main) >= 0.95, "solid",
                         fifelse(get(colLi_main) >= 0.90, "likely",
                                 fifelse(get(colLi_main) >= 0.75, "leaning", "neutral")))]
      
      setorder(part_tbl, pid)
      
      f_pid_csv <- file.path(out_dir, paste0("rq3_", tr, "_", type_nm, "_participants.csv"))
      if (!should_skip(
        paths = f_pid_csv,
        cfg   = cfg,
        type  = "output",
        label = paste0("RQ3 participants ", type_nm, " (", ds, "/", tr, ")")
      )) {
        fwrite(part_tbl, f_pid_csv)
        msg("Saved: ", f_pid_csv)
      }
      
      outputs[[paste(tr, tag, sep = "_")]] <- list(
        fit = f_fit, pid = f_pid, seq = f_seq,
        sequences_csv = f_seq_csv,
        participants_csv = f_pid_csv
      )
    }
  }
  
  invisible(outputs)
}