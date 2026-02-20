# scripts/analysis/23_rq2_sequences_process.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(rstan)

rq2_postprocess_sequences <- function(dataset = "pilot",
                                      treat_fn = "m25",
                                      e = 100,
                                      rho_vec = c(0.05, 0.03, 0.08),
                                      seed = 12345) {
  
  # Load prepared data (for a_bar, s_star, pid levels)
  f_dat <- file.path(path_mod_ds(dataset), paste0("rq2_data_", treat_fn, ".rds"))
  stopifnot(file.exists(f_dat))
  obj <- readRDS(f_dat)
  
  pid_levels <- obj$pid_levels
  seq_levels <- obj$seq_levels
  per_i <- obj$per_i[pid %in% pid_levels]
  
  # Ensure per_i aligned to pid_levels
  setkey(per_i, pid)
  per_i <- per_i[J(pid_levels)]
  stopifnot(all(per_i$pid == pid_levels))
  
  a_bar <- per_i$a_bar
  s_star <- per_i$s_star
  N <- length(pid_levels)
  S <- length(seq_levels)
  
  # Load RQ2 fit
  f_fit <- file.path(path_mod_ds(dataset), paste0("rq2_fit_", treat_fn, ".rds"))
  stopifnot(file.exists(f_fit))
  fit <- readRDS(f_fit)
  post <- rstan::extract(fit)
  
  alpha_draws <- post$alpha        # iters_rq2
  u_draws     <- post$u            # iters_rq2 x N
  b_draws     <- post$b            # iters_rq2 x S
  
  iters_rq2 <- length(alpha_draws)
  stopifnot(ncol(u_draws) == N, ncol(b_draws) == S)
  
  # Load a* draws from HL
  f_astar <- file.path(path_mod_ds(dataset), "a_star_draws_fn.rds")
  stopifnot(file.exists(f_astar))
  astar <- readRDS(f_astar)
  pid_levels_hl <- astar$pid
  a_star_draws <- astar$a_star_draws_fn  # iters_hl x N_hl
  
  # Align HL pids to RQ2 pids
  idx <- match(pid_levels, pid_levels_hl)
  if (anyNA(idx)) {
    bad <- pid_levels[is.na(idx)]
    stop("Some RQ2 pids missing in HL a* draws. Example: ", paste(head(bad, 10), collapse = ", "))
  }
  a_star_draws <- a_star_draws[, idx, drop = FALSE]  # iters_hl x N
  
  iters_hl <- nrow(a_star_draws)
  
  # Match draw counts by resampling HL draws to RQ2 iterations (simple + reproducible)
  set.seed(seed)
  hl_idx <- sample.int(iters_hl, size = iters_rq2, replace = TRUE)
  a_star_k <- a_star_draws[hl_idx, , drop = FALSE]   # iters_rq2 x N
  
  # Compute mu_s^a draws:
  # mu_s^a(k) = (1/e) * mean_i [ (a_bar_i - a*_i(k)) + (alpha(k)+u_i(k)+b_s(k))*s_star_i ]
  mu_sa_draws <- matrix(NA_real_, nrow = iters_rq2, ncol = S)
  
  # Precompute term_i(k) = a_bar_i - a*_i(k)
  base_term <- a_star_k
  for (i in 1:N) base_term[, i] <- a_bar[i] - base_term[, i]  # overwrite as base term
  
  for (s in 1:S) {
    # eta_i(k) = alpha(k) + u_i(k) + b_s(k)
    # For each k, compute mean over i of [ base_term(k,i) + eta_i(k)*s_star_i ]
    # Vectorized over i using sweep:
    eta <- sweep(u_draws, 1, alpha_draws + b_draws[, s], FUN = "+")  # iters x N
    scaled <- sweep(eta, 2, s_star, FUN = "*")                       # iters x N
    mu_sa_draws[, s] <- rowMeans(base_term + scaled) / e
  }
  
  # Summaries
  seq_tbl <- data.table(
    seq = seq_levels,
    mu_a_median = apply(mu_sa_draws, 2, median),
    mu_a_mean   = apply(mu_sa_draws, 2, mean),
    mu_a_q025   = apply(mu_sa_draws, 2, quantile, probs = 0.025),
    mu_a_q975   = apply(mu_sa_draws, 2, quantile, probs = 0.975)
  )
  
  # Under/over probs for each rho
  for (rho in rho_vec) {
    nm_u <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
    nm_o <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
    seq_tbl[, (nm_u) := apply(mu_sa_draws, 2, function(x) mean(x < -rho))]
    seq_tbl[, (nm_o) := apply(mu_sa_draws, 2, function(x) mean(x >  rho))]
  }
  
  setorder(seq_tbl, seq)
  
  out_csv <- file.path(path_out_ds(dataset), paste0("rq2_sequences_", dataset, "_", treat_fn, ".csv"))
  fwrite(seq_tbl, out_csv)
  msg("Saved RQ2 sequence table: ", out_csv)
  
  invisible(list(sequences = seq_tbl, mu_sa_draws = mu_sa_draws))
}

# Example:
# rq2_postprocess_sequences("pilot")
# rq2_postprocess_sequences("main")