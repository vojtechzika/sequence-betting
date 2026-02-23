# scripts/analysis/23_rq2_participants.R
# ------------------------------------------------------------
# RQ2 participant summaries (post-processing only):
# - Requires outputs produced by 21_rq2_stan.R (merged preprocessing + Stan):
#     models/rq2_fit_sequences_<tr>.rds
#     models/rq2_pid_levels_<tr>.rds
#     output/rq2_prepared_<tr>.csv
#
# - Maps posterior draws back to Δa/e units:
#     mu_i^a(k) = (1/e) * [ delta_bar_i + (alpha(k) + u_i(k)) * sd_star_i ]
#   (because mean_s b_s = 0 by sum-to-zero identification)
#
# - Computes:
#     U^a_i = P(mu_i^a < -rho) and O^a_i = P(mu_i^a > rho)
#   for rho in design$rhos$rq2_rho
#
# - Writes:
#     output/rq2_participants_<tr>.csv
# ------------------------------------------------------------

library(data.table)
library(rstan)

rq2_participants <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  run    <- cfg$run
  design <- cfg$design
  
  ds <- as.character(run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # endowment
  stopifnot(!is.null(design$seq), !is.null(design$seq$endowment))
  e <- as.numeric(design$seq$endowment)
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  
  # rho thresholds
  stopifnot(!is.null(design$rhos), !is.null(design$rhos$rq2_rho))
  rho_vec <- as.numeric(design$rhos$rq2_rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  rho_main <- rho_vec[1]
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # Files produced by 21_rq2_stan.R (merged preprocessing + Stan)
    f_fit        <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_prep_csv   <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    
    if (!file.exists(f_fit)) {
      stop("RQ2 fit not found: ", f_fit, "\nRun rq2_stan(cfg) first.")
    }
    if (!file.exists(f_pid_levels)) {
      stop("RQ2 pid_levels not found: ", f_pid_levels, "\nRun rq2_stan(cfg) first.")
    }
    if (!file.exists(f_prep_csv)) {
      stop("RQ2 prepared csv not found: ", f_prep_csv, "\nRun rq2_stan(cfg) first.")
    }
    
    pid_levels <- readRDS(f_pid_levels)
    pid_levels <- as.character(pid_levels)
    stopifnot(length(pid_levels) >= 1L, all(nzchar(pid_levels)))
    
    prep <- fread(f_prep_csv, encoding = "UTF-8")
    prep[, pid := as.character(pid)]
    
    # Expect these columns from the merged preprocessing in 21_
    required_cols <- c("pid", "delta_bar", "sd_star", "n_bets")
    missing_cols <- setdiff(required_cols, names(prep))
    if (length(missing_cols) > 0) {
      stop(
        "Prepared file missing columns: ", paste(missing_cols, collapse = ", "),
        "\nFile: ", f_prep_csv,
        "\nExpected at least: ", paste(required_cols, collapse = ", ")
      )
    }
    
    # participant constants used in mapping back
    pid_map <- prep[, .(
      n_rows = .N,
      
      delta_bar_first = delta_bar[1],
      delta_bar_min   = min(delta_bar, na.rm = TRUE),
      delta_bar_max   = max(delta_bar, na.rm = TRUE),
      
      sd_star_first   = sd_star[1],
      sd_star_min     = min(sd_star, na.rm = TRUE),
      sd_star_max     = max(sd_star, na.rm = TRUE),
      
      n_bets_first    = n_bets[1],
      n_bets_min      = min(n_bets, na.rm = TRUE),
      n_bets_max      = max(n_bets, na.rm = TRUE)
    ), by = pid]
    
    # tolerate tiny numeric jitter; fail if it is substantively inconsistent
    tol <- 1e-10
    
    bad <- pid_map[
      (delta_bar_max - delta_bar_min) > tol |
        (sd_star_max   - sd_star_min)   > tol |
        (n_bets_max    - n_bets_min)    > 0
    ]
    
    if (nrow(bad) > 0) {
      stop(
        "Prepared file has inconsistent participant constants for ",
        nrow(bad), " pid(s). Example(s):\n",
        paste0(
          head(bad$pid, 10),
          " | delta_bar range=",
          sprintf("[%.12g, %.12g]", bad$delta_bar_min, bad$delta_bar_max),
          " | sd_star range=",
          sprintf("[%.12g, %.12g]", bad$sd_star_min, bad$sd_star_max),
          " | n_bets range=",
          sprintf("[%d, %d]", bad$n_bets_min, bad$n_bets_max)
          ,
          collapse = "\n")
      )
    }
    
    # collapse to the constants we actually need
    pid_map <- pid_map[, .(
      pid       = pid,
      delta_bar = delta_bar_first,
      sd_star   = sd_star_first,
      n_bets    = as.integer(n_bets_first)
    )]
    
    setkey(pid_map, pid)
    pid_map <- pid_map[.(pid_levels)]
    if (anyNA(pid_map$pid)) {
      bad <- pid_levels[is.na(pid_map$pid)]
      stop("Could not align pid_map to pid_levels. Example missing pid: ", paste(head(bad, 10), collapse = ", "))
    }
    
    delta_bar_vec <- as.numeric(pid_map$delta_bar)
    sd_star_vec   <- as.numeric(pid_map$sd_star)
    n_bets_vec    <- as.integer(pid_map$n_bets)
    
    # Load fit + extract draws
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    if (is.null(post$alpha)) stop("Fit does not contain 'alpha'. File: ", f_fit)
    if (is.null(post$u))     stop("Fit does not contain 'u'. File: ", f_fit)
    
    alpha <- post$alpha          # iters
    u     <- post$u              # iters x N
    
    iters <- length(alpha)
    N <- length(pid_levels)
    
    stopifnot(is.matrix(u), nrow(u) == iters, ncol(u) == N)
    stopifnot(length(delta_bar_vec) == N, length(sd_star_vec) == N)
    
    # mu_i^a draws in Δa/e units
    mu_i_a_draws <- matrix(NA_real_, nrow = iters, ncol = N)
    for (i in 1:N) {
      mu_i_a_draws[, i] <- (delta_bar_vec[i] + (alpha + u[, i]) * sd_star_vec[i]) / e
    }
    
    part_tbl <- data.table(
      pid         = pid_levels,
      n_bets      = n_bets_vec,
      mu_a_median = apply(mu_i_a_draws, 2, median),
      mu_a_mean   = apply(mu_i_a_draws, 2, mean),
      mu_a_q025   = apply(mu_i_a_draws, 2, quantile, probs = 0.025),
      mu_a_q975   = apply(mu_i_a_draws, 2, quantile, probs = 0.975)
    )
    
    for (rho in rho_vec) {
      nmU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      nmO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      part_tbl[, (nmU) := apply(mu_i_a_draws, 2, function(x) mean(x < -rho))]
      part_tbl[, (nmO) := apply(mu_i_a_draws, 2, function(x) mean(x >  rho))]
    }
    
    colU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    colO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    
    part_tbl[, under_label :=
               fifelse(get(colU) >= 0.95, "solid",
                       fifelse(get(colU) >= 0.90, "likely",
                               fifelse(get(colU) >= 0.75, "leaning", "neutral")))]
    
    part_tbl[, over_label :=
               fifelse(get(colO) >= 0.95, "solid",
                       fifelse(get(colO) >= 0.90, "likely",
                               fifelse(get(colO) >= 0.75, "leaning", "neutral")))]
    
    setorder(part_tbl, pid)
    
    f_part <- file.path(out_dir, paste0("rq2_participants_", tr, ".csv"))
    fwrite(part_tbl, f_part)
    
    msg("Saved: ", f_part)
    
    outputs[[tr]] <- list(
      participants = part_tbl,
      mu_i_a_draws = mu_i_a_draws,
      fit_file = f_fit,
      prep_file = f_prep_csv,
      pid_levels_file = f_pid_levels
    )
  }
  
  invisible(outputs)
}

# Example:
# rq2_participants(cfg)