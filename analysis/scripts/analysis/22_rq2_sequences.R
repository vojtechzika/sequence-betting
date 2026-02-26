# scripts/analysis/22_rq2_sequences.R
# ------------------------------------------------------------
# RQ2 sequences post-processing:
# - Load prepared data (rq2_prepared_<tr>.csv) to recover participant constants:
#     delta_bar_i, sd_star_i (computed on betting trials after exclusions)
# - Load Stan fit (rq2_fit_sequences_<tr>.rds) and level files:
#     rq2_pid_levels_<tr>.rds, rq2_seq_levels_<tr>.rds
# - Map posterior draws back to absolute scale (Δa/e units):
#     mu_s^a(k) = (1/e) * mean_i[ delta_bar_i + (alpha(k)+u_i(k)+b_s(k)) * sd_star_i ]
# - Compute U^a_s = P(mu_s^a < -rho) and O^a_s = P(mu_s^a > rho)
#   for rho in design$rhos$rq2_rho
# - Save: rq2_sequences_calibration_<tr>.csv
# ------------------------------------------------------------

library(data.table)
library(rstan)

rq2_sequences <- function(cfg) {
  
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
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ----- expected artifacts from 21_rq2_stan.R (NO dataset tag in filenames)
    f_fit        <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    f_pid_levels <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f_seq_levels <- file.path(mod_dir, paste0("rq2_seq_levels_", tr, ".rds"))
    f_prep_csv   <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    
    if (!file.exists(f_fit))        stop("RQ2 fit not found: ", f_fit, "\nRun rq2_stan(cfg) first.")
    if (!file.exists(f_pid_levels)) stop("RQ2 pid_levels not found: ", f_pid_levels, "\nRun rq2_stan(cfg) first.")
    if (!file.exists(f_seq_levels)) stop("RQ2 seq_levels not found: ", f_seq_levels, "\nRun rq2_stan(cfg) first.")
    if (!file.exists(f_prep_csv))   stop("RQ2 prepared csv not found: ", f_prep_csv, "\nRun rq2_stan(cfg) first.")
    
    pid_levels <- as.character(readRDS(f_pid_levels))
    seq_levels <- as.character(readRDS(f_seq_levels))
    stopifnot(length(pid_levels) >= 1L, length(seq_levels) >= 1L)
    
    prep <- fread(f_prep_csv, encoding = "UTF-8")
    stopifnot(all(c("pid","seq","delta_bar","sd_star","n_bets") %in% names(prep)))
    prep[, pid := as.character(pid)]
    prep[, seq := as.character(seq)]
    
    # participant constants used in mapping back (must be unique per pid)
    pid_map <- prep[, .(
      delta_bar = unique(delta_bar),
      sd_star   = unique(sd_star),
      n_bets    = unique(n_bets)
    ), by = pid]
    
    # sanity: each pid should have exactly one delta_bar / sd_star / n_bets
    if (pid_map[, any(length(delta_bar) != 1L | length(sd_star) != 1L | length(n_bets) != 1L), by = pid][, any(V1)]) {
      stop("Prepared file has non-unique delta_bar/sd_star/n_bets within pid (should be constant within pid).")
    }
    
    setkey(pid_map, pid)
    pid_map <- pid_map[.(pid_levels)]
    if (anyNA(pid_map$pid)) stop("pid_map could not be aligned to pid_levels for treatment='", tr, "'.")
    
    delta_bar_vec <- as.numeric(pid_map$delta_bar)
    sd_star_vec   <- as.numeric(pid_map$sd_star)
    
    # ----- load fit + extract
    fit <- readRDS(f_fit)
    post <- rstan::extract(fit)
    
    if (is.null(post$alpha)) stop("Fit does not contain alpha.")
    if (is.null(post$u))     stop("Fit does not contain u (participant effects).")
    if (is.null(post$b))     stop("Fit does not contain b (sequence effects).")
    
    alpha <- post$alpha   # iters
    u     <- post$u       # iters x N
    b     <- post$b       # iters x S
    
    iters <- length(alpha)
    N <- length(pid_levels)
    S <- length(seq_levels)
    
    stopifnot(is.matrix(u), is.matrix(b))
    stopifnot(nrow(u) == iters, ncol(u) == N)
    stopifnot(nrow(b) == iters, ncol(b) == S)
    stopifnot(length(delta_bar_vec) == N, length(sd_star_vec) == N)
    
    # ------------------------------------------------------------
    # mu_s^a draws in Δa/e units
    # mu_s^a(k) = (1/e) * mean_i[ delta_bar_i + (alpha+u_i+b_s)*sd_star_i ]
    # ------------------------------------------------------------
    mu_s_a_draws <- matrix(NA_real_, nrow = iters, ncol = S)
    
    # eta_base[k,i] = alpha[k] + u[k,i]
    eta_base <- sweep(u, 1, alpha, FUN = "+")  # iters x N
    
    for (s in 1:S) {
      # eta_mat[k,i] = alpha[k] + u[k,i] + b[k,s]
      eta_mat <- eta_base + b[, s]  # vector recycled across columns -> iters x N
      
      # mapped[k,i] = delta_bar_i + eta_mat[k,i] * sd_star_i
      mapped <- sweep(eta_mat, 2, sd_star_vec, FUN = "*")
      mapped <- sweep(mapped, 2, delta_bar_vec, FUN = "+")
      
      mu_s_a_draws[, s] <- rowMeans(mapped) / e
    }
    
    seq_tbl <- data.table(
      sequence    = seq_levels,
      mu_a_median = apply(mu_s_a_draws, 2, median),
      mu_a_mean   = apply(mu_s_a_draws, 2, mean),
      mu_a_q025   = apply(mu_s_a_draws, 2, quantile, probs = 0.025),
      mu_a_q975   = apply(mu_s_a_draws, 2, quantile, probs = 0.975)
    )
    
    for (rho in rho_vec) {
      nmU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      nmO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      seq_tbl[, (nmU) := apply(mu_s_a_draws, 2, function(x) mean(x < -rho))]
      seq_tbl[, (nmO) := apply(mu_s_a_draws, 2, function(x) mean(x >  rho))]
    }
    
    rho_main <- rho_vec[1]
    colU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    colO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    
    seq_tbl[, under_label :=
              fifelse(get(colU) >= 0.95, "strong",
                      fifelse(get(colU) >= 0.80, "moderate",
                              fifelse(get(colU) >= 0.50, "weak", "neutral")))]
    
    seq_tbl[, over_label :=
              fifelse(get(colO) >= 0.95, "strong",
                      fifelse(get(colO) >= 0.80, "moderate",
                              fifelse(get(colO) >= 0.50, "weak", "neutral")))]
    
    setorder(seq_tbl, sequence)
    
    f_seq_cal <- file.path(out_dir, paste0("rq2_sequences_", tr, ".csv"))
    fwrite(seq_tbl, f_seq_cal)
    
    msg("Saved: ", f_seq_cal)
    
    outputs[[tr]] <- list(
      sequences    = seq_tbl,
      mu_s_a_draws = mu_s_a_draws,
      fit_file     = f_fit,
      prep_file    = f_prep_csv,
      pid_levels_file = f_pid_levels,
      seq_levels_file = f_seq_levels
    )
  }
  
  invisible(outputs)
}

# Example:
# rq2_sequences(cfg)