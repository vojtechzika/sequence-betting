# ============================================================
# 61_ex2-tables.R
#
# PURPOSE
#   Posterior summary of the EX2 model (61_ex2-stan-associations.R) for
#   each treatment: the pooled slopes of optimism and response time and
#   their between-outcome SDs, and the outcome-specific slopes and residual
#   SDs for betting (b), stake deviation (a), welfare loss (c) and Heads
#   choice (h). Each row gives the posterior median, mean, 95% credible
#   interval and, for slopes, P(slope > 0).
#
# INPUT
#   model "ex2": fit-<tr>.rds, data-<tr>.rds
#
# OUTPUT (data/output/csv/)
#   ex2-<tr>-full-summary.csv   block (pooled, outcome), outcome, term,
#                               median, mean, q025, q975, p_gt0,
#                               n (participants of the outcome block)
#
# NOTES
#   - Participant-level outcomes are not tabulated here; they are in
#     rq<k>-<tr>-full-participants.csv of the RQ tables
#   - Slopes are in standard-deviation units of the standardised outcome
#     per standard deviation of the predictor
#   - Output files are overwritten at every run
# ============================================================

ex2_tables <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  labels <- c("b", "a", "c", "h")

  summ <- function(x, slope = TRUE) {
    q <- quantile(x, c(0.025, 0.975), names = FALSE)
    list(median = median(x), mean = mean(x), q025 = q[1], q975 = q[2],
         p_gt0 = if (slope) mean(x > 0) else NA_real_)
  }

  for (tr in tr_vec) {

    fit  <- readRDS(model_path("ex2", paste0("fit-", tr, ".rds")))
    dat  <- readRDS(model_path("ex2", paste0("data-", tr, ".rds")))
    post <- rstan::extract(fit)
    n_k  <- vapply(dat$blocks, function(b) length(b$pid), integer(1L))

    rows <- list(
      cbind(data.table(block = "pooled", outcome = NA_character_, term = "beta_opt_bar", n = NA_integer_),
            as.data.table(summ(post$beta_opt_bar))),
      cbind(data.table(block = "pooled", outcome = NA_character_, term = "beta_rt_bar", n = NA_integer_),
            as.data.table(summ(post$beta_rt_bar))),
      cbind(data.table(block = "pooled", outcome = NA_character_, term = "tau_opt", n = NA_integer_),
            as.data.table(summ(post$tau_opt, FALSE))),
      cbind(data.table(block = "pooled", outcome = NA_character_, term = "tau_rt", n = NA_integer_),
            as.data.table(summ(post$tau_rt, FALSE)))
    )
    for (k in seq_along(labels)) {
      nk <- unname(n_k[labels[k]])
      for (tm in c("beta_opt", "beta_rt", "sigma")) {
        x <- switch(tm, beta_opt = post$beta_opt_k[, k], beta_rt = post$beta_rt_k[, k],
                    sigma = post$sigma_k[, k])
        rows[[length(rows) + 1L]] <- cbind(
          data.table(block = "outcome", outcome = labels[k], term = tm, n = nk),
          as.data.table(summ(x, tm != "sigma")))
      }
    }

    tbl <- rbindlist(rows)
    tbl[, `:=`(treatment = tr, sample = "full")]
    setcolorder(tbl, c("treatment", "sample", "block", "outcome", "term", "n",
                       "median", "mean", "q025", "q975", "p_gt0"))
    f <- file.path(path_csv, paste0("ex2-", tr, "-full-summary.csv"))
    fwrite(tbl, f);  msg("Saved: ", f)
  }

  invisible(TRUE)
}
