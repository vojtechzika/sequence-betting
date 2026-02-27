# ============================================================
# scripts/indices/05_pids_with_positive_a_star.R
#
# Flags participants by normative stake optimality a*(r_i; m) per treatment
# using cached a_star_draws_<tr>.rds (computed earlier from mpl_r_draws.rds).
#
# Per treatment tr:
#   - computes P0(tr)    = P(a* == 0)
#   - computes Pplus(tr) = P(a*  > 0)
#   - applies prereg-style classification using cfg$design$a_flags:
#       betting_normative[[tr]] == TRUE  -> exclude "normative non-betters" if P0   >= tau
#       betting_normative[[tr]] == FALSE -> flag   "normative betters"     if Pplus>= tau
#
# Outputs (per treatment):
#   models/a_star_pid_flags_<tr>.rds    (table + pid vectors per tau)
#   output/a_star_pid_flags_<tr>.csv    (table only)
# ============================================================

library(data.table)

pids_with_positive_a_star <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  
  ds <- as.character(cfg$run$dataset)
  
  # ---- design: treatments & flags ----
  stopifnot(!is.null(cfg$design$seq), !is.null(cfg$design$seq$treatments))
  tr_names <- names(cfg$design$seq$treatments)
  stopifnot(length(tr_names) > 0L, all(nzchar(tr_names)))
  
  stopifnot(!is.null(cfg$design$a_flags))
  af <- cfg$design$a_flags
  
  stopifnot(!is.null(af$betting_normative))
  bn <- af$betting_normative
  stopifnot(all(tr_names %in% names(bn)))
  
  # thresholds: first is main, rest sensitivity (your convention)
  stopifnot(!is.null(af$tau))
  tau_vec <- as.numeric(af$tau)
  stopifnot(length(tau_vec) >= 1L, all(is.finite(tau_vec)), all(tau_vec > 0), all(tau_vec < 1))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_names) {
    
    f_in <- file.path(mod_dir, paste0("a_star_draws_", tr, ".rds"))
    stopifnot(file.exists(f_in))
    
    f_rds <- file.path(mod_dir, paste0("a_star_pid_flags_", tr, ".rds"))
    f_csv <- file.path(out_dir, paste0("a_star_pid_flags_", tr, ".csv"))
    
    # Skip each output independently (so you can delete just one and regenerate it)
    if (should_skip(paths = f_rds, cfg = cfg, type = "output",
                    label = paste0("a* pid flags RDS (", ds, "/", tr, ")"))) {
      f_rds <- NULL
    }
    if (should_skip(paths = f_csv, cfg = cfg, type = "output",
                    label = paste0("a* pid flags CSV (", ds, "/", tr, ")"))) {
      f_csv <- NULL
    }
    
    # If both outputs are being skipped, move on
    if (is.null(f_rds) && is.null(f_csv)) next
    
    obj <- readRDS(f_in)
    stopifnot(is.list(obj), all(c("pid", "a_star_draws") %in% names(obj)))
    
    pid <- as.character(obj$pid)
    a_star_draws <- obj$a_star_draws
    stopifnot(is.matrix(a_star_draws), ncol(a_star_draws) == length(pid))
    
    # Compute P0 / Pplus from draws
    P0    <- colMeans(a_star_draws == 0L, na.rm = TRUE)
    Pplus <- 1 - P0
    
    # Build base table
    tbl <- data.table(
      dataset = ds,
      treatment = tr,
      pid = pid,
      P0 = P0,
      Pplus = Pplus,
      betting_normative = isTRUE(bn[[tr]])
    )
    
    # Classification per tau
    for (tau in tau_vec) {
      tau_nm <- gsub("\\.", "", sprintf("%.2f", tau))
      
      if (isTRUE(bn[[tr]])) {
        # betting benchmark applies -> exclude normative non-betters
        col_flag <- paste0("exclude_nonbetters_tau", tau_nm)
        tbl[, (col_flag) := (P0 >= tau)]
      } else {
        # no-bet benchmark -> flag normative betters (do not exclude by default)
        col_flag <- paste0("flag_betters_tau", tau_nm)
        tbl[, (col_flag) := (Pplus >= tau)]
      }
    }
    
    setorder(tbl, pid)
    
    # Also store reusable pid vectors per tau (source of truth for downstream)
    pid_sets <- list()
    for (tau in tau_vec) {
      tau_nm <- gsub("\\.", "", sprintf("%.2f", tau))
      
      if (isTRUE(bn[[tr]])) {
        col_flag <- paste0("exclude_nonbetters_tau", tau_nm)
        pid_sets[[paste0("pid_keep_tau", tau_nm)]]    <- tbl[get(col_flag) == FALSE, pid]
        pid_sets[[paste0("pid_exclude_tau", tau_nm)]] <- tbl[get(col_flag) == TRUE,  pid]
      } else {
        col_flag <- paste0("flag_betters_tau", tau_nm)
        pid_sets[[paste0("pid_flag_tau", tau_nm)]]    <- tbl[get(col_flag) == TRUE,  pid]
        pid_sets[[paste0("pid_unflag_tau", tau_nm)]]  <- tbl[get(col_flag) == FALSE, pid]
      }
    }
    
    if (!is.null(f_csv)) {
      fwrite(tbl, f_csv)
      msg("Saved: ", f_csv)
    }
    
    if (!is.null(f_rds)) {
      saveRDS(list(
        dataset = ds,
        treatment = tr,
        tau = tau_vec,
        betting_normative = isTRUE(bn[[tr]]),
        table = tbl,
        pid_sets = pid_sets,
        source_a_star_draws = f_in
      ), f_rds)
      msg("Saved: ", f_rds)
    }
    
    outputs[[tr]] <- list(csv = f_csv, rds = f_rds)
  }
  
  invisible(outputs)
}

# Example:
# source(here::here("scripts","indices","05_pids_with_positive_a_star.R"))
# pids_with_positive_a_star(cfg)