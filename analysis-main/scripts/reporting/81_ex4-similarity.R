# ============================================================
# 81_ex4-similarity.R
#
# PURPOSE
#   EX4, preregistered comparison: similarity between the canonical
#   rule-based classifications (rules and scores as in
#   82_ex4-canonical-rules.R) and the anchor-based directional indices of
#   EX1, chi_s and chi_i, with posterior uncertainty propagated through the
#   draws of the indices. For every rule and posterior draw the script
#   computes
#     - the Pearson correlation between the rule scores and chi (across the
#       64 sequences and across participants),
#     - the sign agreement (share of units with a nonzero rule score and a
#       nonzero index on which both have the same sign),
#     - the R-squared of the regression of chi_s on Q_s, for each rule and
#       for all rules jointly,
#   and summarises them by the posterior median, mean, 95% credible interval
#   and, for correlations, P(> 0). For sequences it also reports the
#   posterior probability that the signs of Q_s and chi_s agree.
#   Non-posterior Pearson correlations of the participant rule scores with
#   the posterior median of chi_i, LOT-R and mean log response time are
#   reported alongside as exploratory descriptives.
#
# SAMPLE
#   Full sample of RQ4 (betting trials), the sample on which the draws of
#   chi were computed.
#
# INPUT
#   model "ex1-chi": chi-draws-<tr>-eps<main>.rds (chi_s, chi_i, seq_levels,
#                    pid_levels); main eps = first value of cfg$design$ex1$eps
#   model "mpl": sample/analysis-trials.csv, through analysis_trials()
#   path_der/lotr-scored.csv
#
# OUTPUT (data/output/csv/)
#   ex4-<tr>-sequence-similarity.csv     correlation and agreement by rule
#   ex4-<tr>-sequence-diagnostics.csv    Q_s and P(sign match) by rule, sequence
#   ex4-<tr>-participant-similarity.csv  correlation and agreement by rule
#   ex4-<tr>-sequence-r2.csv             single-rule and joint R-squared
#   ex4-<tr>-rule-predictors.csv         exploratory correlations of Q_i
#
# NOTES
#   - A draw in which a statistic is undefined (constant input) is dropped
#     from its summary; n_draws gives the number of defined draws
#   - The EX1 anchor framework turned out to be uninformative (pure
#     sequences behaviourally neutral), so these results are reported only
#     as the preregistered comparison
#   - Output files are overwritten at every run
# ============================================================

ex4_similarity <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  heads  <- as.character(design$seq$side_labels$heads)
  tails  <- as.character(design$seq$side_labels$tails)
  eps_nm <- gsub("\\.", "", sprintf("%.2f", as.numeric(design$ex1$eps)[1L]))
  rules  <- c("last", "run2", "run3", "imb", "alt")

  summ_post <- function(x, prob = TRUE) {
    x <- x[is.finite(x)]
    if (!length(x)) return(list(median = NA_real_, mean = NA_real_, q025 = NA_real_,
                                q975 = NA_real_, p_gt0 = NA_real_, n_draws = 0L))
    q <- quantile(x, c(0.025, 0.975), names = FALSE)
    list(median = median(x), mean = mean(x), q025 = q[1], q975 = q[2],
         p_gt0 = if (prob) mean(x > 0) else NA_real_, n_draws = length(x))
  }
  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || sd(x[ok]) <= 0 || sd(y[ok]) <= 0) return(NA_real_)
    cor(x[ok], y[ok])
  }
  safe_agree <- function(x, y) {
    ok <- is.finite(x) & is.finite(y) & x != 0 & y != 0
    if (!any(ok)) return(NA_real_)
    mean(sign(x[ok]) == sign(y[ok]))
  }
  safe_r2 <- function(X, y) {
    ok <- is.finite(y) & rowSums(!is.finite(X)) == 0L
    X <- X[ok, , drop = FALSE]; y <- y[ok]
    if (nrow(X) < ncol(X) + 2L || sd(y) <= 0 || any(apply(X, 2L, sd) <= 0)) return(NA_real_)
    fit <- tryCatch(lm(y ~ X), error = function(e) NULL)
    if (is.null(fit)) NA_real_ else summary(fit)$r.squared
  }
  row_of <- function(level, rule, metric, x, prob = TRUE)
    data.table(level = level, rule = rule, metric = metric, as.data.table(summ_post(x, prob)))

  seq_features <- function(s) {
    ch <- strsplit(s, "", fixed = TRUE)[[1L]]
    stopifnot(length(ch) == 6L, all(ch %in% c(heads, tails)))
    x <- ifelse(ch == heads, 1L, -1L)
    list(x6 = x[6L], run = tail(rle(x)$lengths, 1L), B = sum(x), A = sum(x[-1L] != x[-6L]))
  }

  lotr <- fread(file.path(path_der, "lotr-scored.csv"), encoding = "UTF-8")
  lotr <- unique(lotr[, .(pid = as.character(pid), lotr_score = as.numeric(lotr_score))])

  for (tr in tr_vec) {

    chi <- readRDS(model_path("ex1-chi", paste0("chi-draws-", tr, "-eps", eps_nm, ".rds")))
    seq_levels <- as.character(chi$seq_levels); pid_levels <- as.character(chi$pid_levels)
    stopifnot(ncol(chi$chi_s) == length(seq_levels), ncol(chi$chi_i) == length(pid_levels))

    d <- analysis_trials(tr, "rq4-full")
    d <- d[is.finite(stake) & stake > 0 & side %in% c(heads, tails)]
    d[, `:=`(pid = as.character(pid), seq = as.character(seq),
             c_is = fifelse(side == heads, 1L, -1L))]
    stopifnot(setequal(unique(d$pid), pid_levels))

    ft <- lapply(seq_levels, seq_features)
    sf <- data.table(seq = seq_levels, x6 = vapply(ft, `[[`, integer(1L), "x6"),
                     run = vapply(ft, `[[`, integer(1L), "run"),
                     B = vapply(ft, `[[`, integer(1L), "B"), A = vapply(ft, `[[`, integer(1L), "A"))
    d <- merge(d, sf, by = "seq")

    # Trial-level classification: +1 hot-handish, -1 gamblerish, NA undefined
    d[, q_last := as.numeric(c_is * x6)]
    d[, q_run2 := fifelse(run >= 2L, as.numeric(c_is * x6), NA_real_)]
    d[, q_run3 := fifelse(run >= 3L, as.numeric(c_is * x6), NA_real_)]
    d[, q_imb  := fifelse(B != 0L, as.numeric(-c_is * sign(B)), NA_real_)]
    d[, q_alt  := fcase(A <= 2L, as.numeric(c_is * x6), A >= 4L, -as.numeric(c_is * x6),
                        default = NA_real_)]

    # Mean rule score by unit, in the order of the index draws
    score_by <- function(unit, levels_) {
      sapply(rules, function(r) {
        qc <- paste0("q_", r)
        v  <- d[!is.na(get(qc)), .(Q = mean(get(qc))), by = c(unit)]
        v$Q[match(levels_, v[[unit]])]
      })
    }
    Qs <- score_by("seq", seq_levels)   # sequences x rules
    Qi <- score_by("pid", pid_levels)   # participants x rules

    # ---- sequence level ----
    seq_rows <- list(); diag_rows <- list(); r2_rows <- list()
    for (r in rules) {
      qs <- Qs[, r]
      seq_rows[[paste0(r, "c")]] <- row_of("sequence", r, "correlation",
                                           apply(chi$chi_s, 1L, function(z) safe_cor(qs, z)), TRUE)
      seq_rows[[paste0(r, "a")]] <- row_of("sequence", r, "agreement",
                                           apply(chi$chi_s, 1L, function(z) safe_agree(qs, z)), FALSE)
      p_match <- vapply(seq_along(seq_levels), function(j) {
        if (!is.finite(qs[j]) || qs[j] == 0) return(NA_real_)
        mean(sign(qs[j]) == sign(chi$chi_s[, j]) & chi$chi_s[, j] != 0)
      }, numeric(1L))
      diag_rows[[r]] <- data.table(rule = r, sequence = seq_levels, seq_id = seq_along(seq_levels),
                                   Q_seq = qs, p_match = p_match)
      r2_rows[[r]] <- cbind(data.table(level = "sequence", model = "single_rule", rule = r,
                                       metric = "r_squared"),
                            as.data.table(summ_post(apply(chi$chi_s, 1L,
                              function(z) safe_r2(Qs[, r, drop = FALSE], z)), FALSE)))
    }
    r2_rows[["joint"]] <- cbind(data.table(level = "sequence", model = "joint_rules",
                                           rule = "all", metric = "r_squared"),
                                as.data.table(summ_post(apply(chi$chi_s, 1L,
                                  function(z) safe_r2(Qs, z)), FALSE)))

    # ---- participant level ----
    pid_rows <- list()
    for (r in rules) {
      qi <- Qi[, r]
      pid_rows[[paste0(r, "c")]] <- row_of("participant", r, "correlation",
                                           apply(chi$chi_i, 1L, function(z) safe_cor(qi, z)), TRUE)
      pid_rows[[paste0(r, "a")]] <- row_of("participant", r, "agreement",
                                           apply(chi$chi_i, 1L, function(z) safe_agree(qi, z)), FALSE)
    }

    # ---- exploratory correlates of the participant rule scores ----
    rt <- d[is.finite(screen_ms) & screen_ms > 0, .(rt_log_mean = mean(log(screen_ms))), by = pid]
    pt <- merge(data.table(pid = pid_levels), merge(rt, lotr, by = "pid", all = TRUE),
                by = "pid", all.x = TRUE)
    chi_med <- apply(chi$chi_i, 2L, median)
    pred <- rbindlist(lapply(rules, function(r) {
      qi <- Qi[, r]
      rbindlist(list(
        data.table(rule = r, predictor = "chi_i_median", estimate = safe_cor(qi, chi_med)),
        data.table(rule = r, predictor = "lotr_score",   estimate = safe_cor(qi, pt$lotr_score)),
        data.table(rule = r, predictor = "rt_log_mean",  estimate = safe_cor(qi, pt$rt_log_mean))))
    }))

    tag <- function(x) { x[, `:=`(treatment = tr, sample = "full")]; x }
    out <- list(
      "sequence-similarity"    = tag(rbindlist(seq_rows)),
      "sequence-diagnostics"  = tag(rbindlist(diag_rows)),
      "participant-similarity" = tag(rbindlist(pid_rows)),
      "sequence-r2"           = tag(rbindlist(r2_rows)),
      "rule-predictors"       = tag(pred)
    )
    for (nm in names(out)) {
      f <- file.path(path_csv, paste0("ex4-", tr, "-", nm, ".csv"))
      fwrite(out[[nm]], f);  msg("Saved: ", f)
    }
  }

  invisible(TRUE)
}
