# scripts/analysis/33_rq3_participants.R
# ------------------------------------------------------------
# RQ3 participant summaries (works for hurdle-gamma and gaussian fits)
#
# - Loads rq3 fit + pid_levels + seq_levels
# - Detects fitted model by checking parameter names in rstan::extract()
# - Computes participant-level expected welfare loss share:
#     mu_i^c = E_s[ E(y | i,s) ]
#   hurdle-gamma: E(y|i,s) = (1 - pi0_is) * mean_pos_is
#     pi0_is = inv_logit(a0 + u0_i + b0_s)
#     mean_pos_is = exp(ap + up_i + bp_s)
#   gaussian: E(y|i,s) = alpha + u_i + b_s
# - Computes L_i(rho) = P(mu_i^c > rho) for rho in design$rhos$rq3_rho
# - Saves: data/clean/<ds>/output/rq3_participants_<tr>.csv
# ------------------------------------------------------------

library(data.table)
library(rstan)

rq3_participants <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  stopifnot(!is.null(cfg$design$rhos), !is.null(cfg$design$rhos$rq3_rho))
  rho_vec <- as.numeric(cfg$design$rhos$rq3_rho)
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0), all(rho_vec < 1))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, ".rds"))
    f_pid <- file.path(mod_dir, paste0("rq3_pid_levels_", tr, ".rds"))
    f_seq <- file.path(mod_dir, paste0("rq3_seq_levels_", tr, ".rds"))
    
    if (!file.exists(f_fit)) stop("RQ3 fit not found: ", f_fit)
    if (!file.exists(f_pid)) stop("RQ3 pid_levels not found: ", f_pid)
    if (!file.exists(f_seq)) stop("RQ3 seq_levels not found: ", f_seq)
    
    fit <- readRDS(f_fit)
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    
    post <- rstan::extract(fit)
    
    # ------------------------------------------------------------
    # Detect which RQ3 model was fitted
    # ------------------------------------------------------------
    is_hurdle <- all(c("a0","u0","b0","ap","up","bp") %in% names(post))
    is_gauss  <- all(c("alpha","u","b") %in% names(post))
    
    if (!is_hurdle && !is_gauss) {
      stop(
        "RQ3 fit missing expected parameter sets.\n",
        "Found names: ", paste(sort(names(post)), collapse = ", "), "\n",
        "Expected either hurdle-gamma {a0,u0,b0,ap,up,bp} or gaussian {alpha,u,b}."
      )
    }
    
    # ------------------------------------------------------------
    # Compute mu_i draws
    # ------------------------------------------------------------
    if (is_hurdle) {
      # Dimensions: u0 iters x N, b0 iters x S, etc.
      a0 <- post$a0           # iters
      ap <- post$ap           # iters
      u0 <- post$u0           # iters x N
      up <- post$up           # iters x N
      b0 <- post$b0           # iters x S
      bp <- post$bp           # iters x S
      
      stopifnot(is.matrix(u0), is.matrix(up), is.matrix(b0), is.matrix(bp))
      iters <- length(a0)
      N <- length(pid_levels)
      S <- length(seq_levels)
      
      stopifnot(nrow(u0) == iters, ncol(u0) == N)
      stopifnot(nrow(up) == iters, ncol(up) == N)
      stopifnot(nrow(b0) == iters, ncol(b0) == S)
      stopifnot(nrow(bp) == iters, ncol(bp) == S)
      
      mu_i_draws <- matrix(NA_real_, nrow = iters, ncol = N)
      
      # For each draw k and participant i:
      # mu_i(k) = mean_s (1 - inv_logit(a0+u0_i+b0_s)) * exp(ap+up_i+bp_s)
      for (i in 1:N) {
        logit_pi0 <- sweep(b0, 1, a0 + u0[, i], FUN = "+")  # iters x S
        pi0 <- 1 / (1 + exp(-logit_pi0))
        mean_pos <- exp(sweep(bp, 1, ap + up[, i], FUN = "+")) # iters x S
        Ey <- (1 - pi0) * mean_pos
        mu_i_draws[, i] <- rowMeans(Ey)
      }
      
    } else {
      alpha <- post$alpha   # iters
      u     <- post$u       # iters x N
      b     <- post$b       # iters x S
      
      stopifnot(is.matrix(u), is.matrix(b))
      iters <- length(alpha)
      N <- length(pid_levels)
      S <- length(seq_levels)
      
      stopifnot(nrow(u) == iters, ncol(u) == N)
      stopifnot(nrow(b) == iters, ncol(b) == S)
      
      mu_i_draws <- matrix(NA_real_, nrow = iters, ncol = N)
      
      # E_s[alpha + u_i + b_s] = alpha + u_i + mean_s(b_s)
      # With sum-to-zero b, mean_s(b_s) ~ 0, but compute it explicitly.
      b_mean <- rowMeans(b)  # iters
      for (i in 1:N) {
        mu_i_draws[, i] <- alpha + u[, i] + b_mean
      }
    }
    
    # ------------------------------------------------------------
    # Build participant table
    # ------------------------------------------------------------
    part_tbl <- data.table(
      pid        = pid_levels,
      mu_c_median = apply(mu_i_draws, 2, median),
      mu_c_mean   = apply(mu_i_draws, 2, mean),
      mu_c_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
      mu_c_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975)
    )
    
    for (rho in rho_vec) {
      nm <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      part_tbl[, (nm) := apply(mu_i_draws, 2, function(x) mean(x > rho))]
    }
    
    # Labels using main rho (first element)
    rho_main <- rho_vec[1]
    colL <- paste0("L_rho_", gsub("\\.", "", sprintf("%.2f", rho_main)))
    
    part_tbl[, loss_label :=
               fifelse(get(colL) >= 0.95, "solid",
                       fifelse(get(colL) >= 0.90, "likely",
                               fifelse(get(colL) >= 0.75, "leaning", "neutral")))]
    
    setorder(part_tbl, pid)
    
    f_out <- file.path(out_dir, paste0("rq3_participants_", tr, ".csv"))
    fwrite(part_tbl, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- list(
      participants = part_tbl,
      mu_i_draws = mu_i_draws,
      model = if (is_hurdle) "hurdle_gamma" else "gaussian"
    )
  }
  
  invisible(outputs)
}

# Example:
# rq3_participants(cfg)