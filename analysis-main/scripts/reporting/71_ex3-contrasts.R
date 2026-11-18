# ============================================================
# 71_ex3-contrasts.R
#
# PURPOSE
#   EX3: posterior contrasts between samples for betting (RQ1) and welfare
#   loss (RQ3), and the rank correlation of the sequence means for side
#   choice (RQ4), the contrasts specified in the preregistration (Section
#   2.2.3), extended by a robustness comparison within the FN treatment:
#     fnfull-fpfull   the full sample of the FN treatment minus the full
#                     sample of the FP treatment (the preregistered contrast)
#     fnconf-fnfull   the confirmatory sample of the FN treatment (consistent
#                     Holt-Laury choices) minus the full sample of the FN
#                     treatment, which shows the effect of the restriction to
#                     consistent participants (the FP treatment has no
#                     confirmatory sample, so the confirmatory sample is not
#                     contrasted with FP)
#   The stake deviation (RQ2) is contrasted only within the FN treatment
#   (fnconf-fnfull), because its benchmark is determined by the treatment
#   (a* = 0 for nearly all participants in FP). RQ4 is not contrasted: it
#   describes direction conditional on betting, which differs in selection
#   between treatments, so that contrasts of its levels would compare
#   differently selected bettors. Only the rank correlation of the sequence
#   means (3) is computed for RQ4 (fnfull-fpfull), which asks whether the
#   sequences are ordered alike in the two treatments and does not compare
#   levels. For each research
#   question and comparison the script computes, draw by draw,
#     (1) the contrast of the sequence means, Delta_s = mu_s(first) - mu_s(second),
#         for each of the 64 sequences,
#     (2) the contrast of the averages of the participant means, a contrast
#         defined at the group level,
#     (3) the Spearman correlation between the sequence means of the two
#         samples across the 64 sequences.
#   Each contrast is summarised by its posterior median, mean, 95% credible
#   interval, P(Delta > 0) and P(|Delta| > rho) for the practical-equivalence
#   threshold rho of the research question (first value of its design entry).
#   The thresholds preregistered for EX3 are reported in addition when they
#   differ from rho: 0.05 for betting (rho = 0.10) and 0.03 for welfare loss
#   (rho = 0.05, which is also preregistered); P(|Delta| > rho_alt) is given
#   by rho_alt and p_abs_gt_rho_alt (NA for the stake deviation, for which the
#   preregistration states no EX3 threshold).
#
# SAMPLE
#   Full and confirmatory samples as named above; EX3 is exploratory.
#
# MODELS
#   RQ1 "rq1": the model selected in selected-models.csv.
#   RQ2 "rq2": the model selected in selected-models.csv (mu_a, mu_a_i).
#   RQ4 "rq4": the model selected in selected-models.csv (mu_h, mu_h_i), full
#   sample, rank correlation only.
#   RQ3: the hurdle-Gamma model in all samples (folder rq3-hurdle-gamma),
#   whatever is selected for "rq3", so that the estimand E[y] is the same in
#   all arms. The agreement of the hurdle-Gamma model with the Gamma-only and
#   Gaussian models is reported in 35_rq3-likelihood-comparison.R.
#
# INPUT
#   fit-, pid-levels-, seq-levels-<tr>-<sample>.rds of the models above
#
# OUTPUT (data/output/csv/), <tag> = fnfull-fpfull, fnconf-fnfull
#   ex3-rq<k>-sequence-contrasts-<tag>.csv      one row per sequence
#   ex3-rq<k>-participant-contrasts-<tag>.csv   group-level contrast
#   ex3-rq<k>-spearman-<tag>.csv                rank correlation of sequence means
#                                               k = 1, 2, 3 (k = 2 only fnconf-fnfull)
#   ex3-rq4-spearman-fnfull-fpfull.csv          rank correlation only (no contrasts)
#
# NOTES
#   - The fits of the samples are independent, so draws are paired by row
#     index over the smaller number of draws
#   - The confirmatory sample is a subset of the full sample of the FN
#     treatment, but the two fits are independent, so the contrast between
#     them describes the difference of the estimates, not of the data
#   - Spearman correlations are computed within each draw and summarised
#     over draws
#   - Output files are overwritten at every run
# ============================================================

ex3_contrasts <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) >= 2L)
  fn <- tr_vec[1L]; fp <- tr_vec[2L]
  design <- cfg$design

  specs <- list(
    list(rq = "rq1", model = "rq1",              mu_s = "mu_b", mu_i = "mu_b_i",
         rho = design$rq1$rho[1L], rho_alt = 0.05),
    list(rq = "rq2", model = "rq2",              mu_s = "mu_a", mu_i = "mu_a_i",
         rho = design$rq2$rho[1L], rho_alt = NA_real_, only = "fnconf-fnfull"),
    list(rq = "rq3", model = "rq3-hurdle-gamma", mu_s = "mu_c", mu_i = "mu_c_i",
         rho = design$rq3$rho[1L], rho_alt = 0.03),
    # RQ4: rank correlation of the sequence means only, FN against FP
    list(rq = "rq4", model = "rq4",              mu_s = "mu_h", mu_i = "mu_h_i",
         rho = NA_real_, rho_alt = NA_real_, only = "fnfull-fpfull", rank_only = TRUE)
  )

  # comparisons: first arm minus second arm
  comparisons <- list(
    list(tag = "fnfull-fpfull", a = c(fn, "full"),         b = c(fp, "full")),
    list(tag = "fnconf-fnfull", a = c(fn, "confirmatory"), b = c(fn, "full"))
  )

  # posterior summary of a contrast; rho_alt is the additional preregistered threshold (NA if none)
  summ_delta <- function(x, rho, rho_alt = NA_real_) {
    q <- quantile(x, c(0.025, 0.975), names = FALSE)
    data.table(median = median(x), mean = mean(x), q025 = q[1], q975 = q[2],
               p_gt0 = mean(x > 0), rho = rho, p_abs_gt_rho = mean(abs(x) > rho),
               rho_alt = rho_alt,
               p_abs_gt_rho_alt = if (is.na(rho_alt)) NA_real_ else mean(abs(x) > rho_alt))
  }

  for (sp in specs) {

    # posterior draws of one arm (treatment and sample)
    load_arm <- function(tr, smp) {
      f   <- function(stem) model_path(sp$model, paste0(stem, "-", tr, "-", smp, ".rds"))
      fit <- readRDS(f("fit"))
      post <- rstan::extract(fit, pars = c(sp$mu_s, sp$mu_i))
      if (sp$rq == "rq3") stopifnot(!is.null(rstan::extract(fit, pars = "a0")$a0))
      list(mu_s = post[[sp$mu_s]], mu_i = post[[sp$mu_i]],
           seq = as.character(readRDS(f("seq-levels"))),
           pid = as.character(readRDS(f("pid-levels"))))
    }
    arms <- list()
    arm  <- function(x) {
      key <- paste(x, collapse = "-")
      if (is.null(arms[[key]])) arms[[key]] <<- load_arm(x[1L], x[2L])
      arms[[key]]
    }

    for (cp in comparisons) {

      if (!is.null(sp$only) && cp$tag != sp$only) next

      a <- arm(cp$a); b <- arm(cp$b)
      if (!identical(a$seq, b$seq))
        stop("EX3 ", sp$rq, " ", cp$tag, ": sequence order differs between the arms.", call. = FALSE)
      Tn <- min(nrow(a$mu_s), nrow(b$mu_s))
      meta <- data.table(rq = sp$rq, contrast = cp$tag, treatment_a = cp$a[1L], sample_a = cp$a[2L],
                         treatment_b = cp$b[1L], sample_b = cp$b[2L], model = sp$model)

      rank_only <- isTRUE(sp$rank_only)

      # (1) sequence contrasts
      if (!rank_only) {
      seq_tbl <- rbindlist(lapply(seq_along(a$seq), function(s)
        cbind(data.table(sequence = a$seq[s], seq_id = s),
              summ_delta(a$mu_s[seq_len(Tn), s] - b$mu_s[seq_len(Tn), s], sp$rho, sp$rho_alt))))
      seq_tbl[, n_draws := Tn]
      seq_tbl <- cbind(meta[rep(1L, nrow(seq_tbl))], seq_tbl)
      f1 <- file.path(path_csv, paste0("ex3-", sp$rq, "-sequence-contrasts-", cp$tag, ".csv"))
      fwrite(seq_tbl, f1);  msg("Saved: ", f1)

      # (2) group-level participant contrast
      d_pid <- rowMeans(a$mu_i)[seq_len(Tn)] - rowMeans(b$mu_i)[seq_len(Tn)]
      pid_tbl <- cbind(meta, n_pid_a = length(a$pid), n_pid_b = length(b$pid), n_draws = Tn,
                       summ_delta(d_pid, sp$rho, sp$rho_alt))
      f2 <- file.path(path_csv, paste0("ex3-", sp$rq, "-participant-contrasts-", cp$tag, ".csv"))
      fwrite(pid_tbl, f2);  msg("Saved: ", f2)
      }

      # (3) rank correlation of the sequence means
      rho_d <- vapply(seq_len(Tn), function(t) cor(a$mu_s[t, ], b$mu_s[t, ], method = "spearman"), numeric(1))
      q <- quantile(rho_d, c(0.025, 0.975), names = FALSE)
      spm <- cbind(meta, data.table(n_draws = Tn, median = median(rho_d), mean = mean(rho_d),
                                    q025 = q[1], q975 = q[2], p_gt0 = mean(rho_d > 0)))
      f3 <- file.path(path_csv, paste0("ex3-", sp$rq, "-spearman-", cp$tag, ".csv"))
      fwrite(spm, f3);  msg("Saved: ", f3)
    }
  }

  invisible(TRUE)
}
