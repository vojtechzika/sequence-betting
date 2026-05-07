# ============================================================
# 04_pids_with_positive_a_star.R
#
# PURPOSE
#   Flags participants by normative stake optimality per treatment
#   using cached a_star_draws. Computes posterior probabilities of
#   zero and positive optimal stake and classifies participants as
#   normative non-betters or normative betters. When consistent_only
#   = TRUE, additionally processes full-sample draws for descriptive
#   purposes.
#
# INPUT
#   path_mod/a_star_draws_<tr>.rds       -- primary (consistent-only or full)
#   path_mod/a_star_draws_<tr>_full.rds  -- full-sample descriptive,
#                                           only if consistent_only = TRUE
#
# OUTPUT
#   path_mod/a_star_pid_flags_<tr>.rds       -- primary results + pid sets
#   path_out/a_star_pid_flags_<tr>.csv       -- primary classification table
#   path_mod/a_star_pid_flags_<tr>_full.rds  -- full-sample results,
#   path_out/a_star_pid_flags_<tr>_full.csv  -- only if consistent_only = TRUE
#
# NOTES
#   - betting_normative[[tr]] == TRUE  -> exclude normative non-betters
#   - betting_normative[[tr]] == FALSE -> flag normative betters
#   - Classification thresholds from cfg$design$a_flags$tau
#   - First tau value is primary; rest are sensitivity checks
#   - RDS and CSV are skipped independently via overwrite flags
# ============================================================
pids_with_positive_a_star <- function(cfg) {
  
  stopifnot(!is.null(cfg$design$seq), !is.null(cfg$design$seq$treatments))
  tr_names <- names(cfg$design$seq$treatments)
  stopifnot(length(tr_names) > 0L, all(nzchar(tr_names)))
  
  stopifnot(!is.null(cfg$design$a_flags))
  af <- cfg$design$a_flags
  bn <- af$betting_normative
  stopifnot(!is.null(bn), all(tr_names %in% names(bn)))
  
  tau_vec <- as.numeric(af$tau)
  stopifnot(length(tau_vec) >= 1L, all(is.finite(tau_vec)),
            all(tau_vec > 0), all(tau_vec < 1))
  
  consistent_only <- isTRUE(cfg$run$consistent_only)
  outputs         <- list()
  
  for (tr in tr_names) {
    
    suffix_map <- if (consistent_only) c("", "_full") else c("")
    
    for (sfx in suffix_map) {
      
      f_in  <- file.path(path_mod, paste0("a_star_draws_",    tr, sfx, ".rds"))
      f_rds <- file.path(path_mod, paste0("a_star_pid_flags_", tr, sfx, ".rds"))
      f_csv <- file.path(path_out, paste0("a_star_pid_flags_", tr, sfx, ".csv"))
      
      if (!file.exists(f_in)) {
        msg("Skipping (input not found): ", f_in)
        next
      }
      
      skip_rds <- should_skip(paths = f_rds, cfg = cfg, type = "model",
                              label = paste0("a* pid flags RDS (", tr, sfx, ")"))
      skip_csv <- should_skip(paths = f_csv, cfg = cfg, type = "output",
                              label = paste0("a* pid flags CSV (", tr, sfx, ")"))
      
      if (skip_rds && skip_csv) next
      
      obj <- readRDS(f_in)
      stopifnot(is.list(obj), all(c("pid", "a_star_draws") %in% names(obj)))
      
      pid          <- as.character(obj$pid)
      a_star_draws <- obj$a_star_draws
      stopifnot(is.matrix(a_star_draws), ncol(a_star_draws) == length(pid))
      
      P0    <- colMeans(a_star_draws == 0L, na.rm = TRUE)
      Pplus <- 1 - P0
      
      tbl <- data.table(
        treatment         = tr,
        pid               = pid,
        P0                = P0,
        Pplus             = Pplus,
        betting_normative = isTRUE(bn[[tr]])
      )
      
      pid_sets <- list()
      
      for (tau in tau_vec) {
        tau_nm <- gsub("\\.", "", sprintf("%.2f", tau))
        
        if (isTRUE(bn[[tr]])) {
          col_flag <- paste0("exclude_nonbetters_tau", tau_nm)
          tbl[, (col_flag) := (P0 >= tau)]
          pid_sets[[paste0("pid_keep_tau",    tau_nm)]] <- tbl[get(col_flag) == FALSE, pid]
          pid_sets[[paste0("pid_exclude_tau", tau_nm)]] <- tbl[get(col_flag) == TRUE,  pid]
        } else {
          col_flag <- paste0("flag_betters_tau", tau_nm)
          tbl[, (col_flag) := (Pplus >= tau)]
          pid_sets[[paste0("pid_flag_tau",   tau_nm)]] <- tbl[get(col_flag) == TRUE,  pid]
          pid_sets[[paste0("pid_unflag_tau", tau_nm)]] <- tbl[get(col_flag) == FALSE, pid]
        }
      }
      
      setorder(tbl, pid)
      
      if (!skip_csv) {
        fwrite(tbl, f_csv)
        msg("Saved: ", f_csv)
      }
      
      if (!skip_rds) {
        saveRDS(list(
          treatment           = tr,
          sample              = if (sfx == "") "primary" else "full",
          tau                 = tau_vec,
          betting_normative   = isTRUE(bn[[tr]]),
          table               = tbl,
          pid_sets            = pid_sets,
          source_a_star_draws = f_in
        ), f_rds)
        msg("Saved: ", f_rds)
      }
      
      outputs[[paste0(tr, sfx)]] <- list(csv = f_csv, rds = f_rds)
    }
  }
  
  invisible(outputs)
}