# ============================================================
# 61_ex2-tables.R
#
# PURPOSE
#   Posterior summary of the EX2 model (61_ex2-stan-associations.R) for
#   each treatment, and the correlation of the residuals between outcomes: the pooled slopes of optimism and response time and
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
#   ex2-<tr>-full-residual-corr.csv
#                               posterior of the correlation of the model
#                               residuals for each pair of outcomes (term
#                               "corr", n = participants in both blocks) and
#                               of the factor by which the standard error of
#                               an equally weighted mean of the four slopes
#                               exceeds its value under independent residuals
#                               (term "se_inflation")
#
# NOTES
#   - Participant-level outcomes are not tabulated here; they are in
#     rq<k>-<tr>-full-participants.csv of the RQ tables
#   - Slopes are in standard-deviation units of the standardised outcome
#     per standard deviation of the predictor
#   - Residual correlations are computed after estimation and do not enter
#     the model, which assumes independent residuals. For each of Sres
#     posterior draws, the residual of participant i on outcome k is the
#     standardised outcome draw minus the fitted value; draw t of the
#     outcome (cycled over the Trep stored draws) is paired across outcomes.
#     The outcomes come from separate fits, so draws are independent given
#     the participant, and the correlations describe the residuals as the
#     model sees them, including the upstream estimation noise
#   - se_inflation = sqrt(sum of all entries of the 4 x 4 residual
#     correlation matrix / 4); it holds for slopes of equal sign and is an upper
#     bound otherwise
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

    # ---- Residual correlations between outcomes (post-processing) ----
    blk <- dat$blocks
    if (is.null(blk[[1]]$y)) stop("data-", tr, ".rds has no outcome draws; rerun 61_ex2-stan-associations.R.", call. = FALSE)
    S    <- dim(post$alpha_k)[1]
    Sres <- min(S, 1000L)
    sidx <- round(seq(1, S, length.out = Sres))
    pairs <- t(combn(seq_along(labels), 2L))
    cor_draw <- matrix(NA_real_, Sres, nrow(pairs))
    infl     <- numeric(Sres)
    common   <- vapply(seq_len(nrow(pairs)), function(j)
      length(intersect(blk[[labels[pairs[j, 1]]]]$pid, blk[[labels[pairs[j, 2]]]]$pid)), integer(1L))

    for (s in seq_len(Sres)) {
      ss <- sidx[s]
      res <- lapply(seq_along(labels), function(k) {
        b <- blk[[labels[k]]]
        tt <- ((s - 1L) %% nrow(b$y)) + 1L
        fit_k <- post$alpha_k[ss, k] + post$beta_opt_k[ss, k] * b$Zopt + post$beta_rt_k[ss, k] * b$Zrt
        setNames(b$y[tt, ] - fit_k, b$pid)
      })
      R <- diag(length(labels))
      for (j in seq_len(nrow(pairs))) {
        a <- res[[pairs[j, 1]]];  b <- res[[pairs[j, 2]]]
        ids <- intersect(names(a), names(b))
        r <- cor(a[ids], b[ids])
        cor_draw[s, j] <- r
        R[pairs[j, 1], pairs[j, 2]] <- r;  R[pairs[j, 2], pairs[j, 1]] <- r
      }
      infl[s] <- sqrt(sum(R) / length(labels))
    }

    rc <- rbindlist(c(
      lapply(seq_len(nrow(pairs)), function(j)
        cbind(data.table(outcome = paste(labels[pairs[j, ]], collapse = "-"), term = "corr", n = common[j]),
              as.data.table(summ(cor_draw[, j], FALSE)))),
      list(cbind(data.table(outcome = "pooled", term = "se_inflation", n = NA_integer_),
                 as.data.table(summ(infl, FALSE))))))
    rc[, `:=`(treatment = tr, sample = "full")]
    setcolorder(rc, c("treatment", "sample", "outcome", "term", "n", "median", "mean", "q025", "q975"))
    rc[, p_gt0 := NULL]
    f <- file.path(path_csv, paste0("ex2-", tr, "-full-residual-corr.csv"))
    fwrite(rc, f);  msg("Saved: ", f)
  }

  invisible(TRUE)
}
