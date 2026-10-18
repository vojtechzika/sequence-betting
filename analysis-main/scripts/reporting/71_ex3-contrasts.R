# ============================================================
# 71_ex3-contrasts.R
#
# PURPOSE
#   EX3: posterior contrasts between the treatments (first treatment of
#   cfg$run$treatment minus the second; FN minus FP in the design) for
#   betting (RQ1) and welfare loss (RQ3), the contrasts specified in the
#   preregistration (Section 2.2.3). Contrasts for RQ2 and RQ4 (Appendix C.8
#   of the preregistration) are not computed: the RQ2 benchmark is
#   determined by the treatment (a* = 0 for nearly all participants in FP)
#   and RQ4 describes direction conditional on betting, which differs in
#   selection between treatments. For each research question the script
#   computes, draw by draw,
#     (1) the contrast of the sequence means, Delta_s = mu_s(FN) - mu_s(FP),
#         for each of the 64 sequences,
#     (2) the contrast of the treatment averages of the participant means,
#         a between-subject contrast defined at the group level,
#     (3) the Spearman correlation between the FN and FP sequence means
#         across the 64 sequences.
#   Each contrast is summarised by its posterior median, mean, 95% credible
#   interval, P(Delta > 0) and P(|Delta| > rho) for the practical-equivalence
#   threshold rho of the research question (first value of its design entry).
#
# SAMPLE
#   The full sample of both treatments (no Holt-Laury or normative-better
#   restriction); EX3 is exploratory.
#
# MODELS
#   RQ1 "rq1": the model selected in selected-models.csv.
#   RQ3: the hurdle-Gamma model in both treatments
#   (folder rq3-hurdle-gamma), whatever is selected for "rq3", so that the
#   estimand E[y] is the same in both arms. The agreement of the hurdle-Gamma
#   model with the Gamma-only and Gaussian models is reported in
#   35_rq3-likelihood-comparison.R.
#
# INPUT
#   fit-, pid-levels-, seq-levels-<tr>-full.rds of the models above
#
# OUTPUT (data/output/csv/)
#   ex3-rq<k>-sequence-contrasts.csv      one row per sequence
#   ex3-rq<k>-participant-contrasts.csv   group-level contrast
#   ex3-rq<k>-spearman.csv                FN-FP rank correlation of sequence means
#                                         k = 1, 3
#
# NOTES
#   - The FN and FP fits are independent, so draws are paired by row index
#     over the smaller number of draws
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
         rho = design$rq1$rho[1L]),
    list(rq = "rq3", model = "rq3-hurdle-gamma", mu_s = "mu_c", mu_i = "mu_c_i",
         rho = design$rq3$rho[1L])
  )

  summ_delta <- function(x, rho) {
    q <- quantile(x, c(0.025, 0.975), names = FALSE)
    data.table(median = median(x), mean = mean(x), q025 = q[1], q975 = q[2],
               p_gt0 = mean(x > 0), rho = rho, p_abs_gt_rho = mean(abs(x) > rho))
  }

  for (sp in specs) {

    load_tr <- function(tr) {
      f <- function(stem) model_path(sp$model, paste0(stem, "-", tr, "-full.rds"))
      post <- rstan::extract(readRDS(f("fit")), pars = c(sp$mu_s, sp$mu_i))
      list(mu_s = post[[sp$mu_s]], mu_i = post[[sp$mu_i]],
           seq = as.character(readRDS(f("seq-levels"))),
           pid = as.character(readRDS(f("pid-levels"))))
    }
    a <- load_tr(fn); b <- load_tr(fp)
    if (!identical(a$seq, b$seq))
      stop("EX3 ", sp$rq, ": sequence order differs between treatments.", call. = FALSE)
    if (sp$rq == "rq3") {
      chk <- rstan::extract(readRDS(model_path(sp$model, paste0("fit-", fn, "-full.rds"))), pars = "a0")
      stopifnot(!is.null(chk$a0))
    }
    Tn <- min(nrow(a$mu_s), nrow(b$mu_s))

    # (1) sequence contrasts
    seq_tbl <- rbindlist(lapply(seq_along(a$seq), function(s)
      cbind(data.table(sequence = a$seq[s], seq_id = s),
            summ_delta(a$mu_s[seq_len(Tn), s] - b$mu_s[seq_len(Tn), s], sp$rho))))
    seq_tbl[, `:=`(rq = sp$rq, treatment_fn = fn, treatment_fp = fp, model = sp$model, n_draws = Tn)]
    setcolorder(seq_tbl, c("rq", "sequence", "seq_id", "treatment_fn", "treatment_fp", "model", "n_draws"))
    f1 <- file.path(path_csv, paste0("ex3-", sp$rq, "-sequence-contrasts.csv"))
    fwrite(seq_tbl, f1);  msg("Saved: ", f1)

    # (2) group-level participant contrast
    d_pid <- rowMeans(a$mu_i)[seq_len(Tn)] - rowMeans(b$mu_i)[seq_len(Tn)]
    pid_tbl <- cbind(data.table(rq = sp$rq, treatment_fn = fn, treatment_fp = fp, model = sp$model,
                                n_pid_fn = length(a$pid), n_pid_fp = length(b$pid), n_draws = Tn),
                     summ_delta(d_pid, sp$rho))
    f2 <- file.path(path_csv, paste0("ex3-", sp$rq, "-participant-contrasts.csv"))
    fwrite(pid_tbl, f2);  msg("Saved: ", f2)

    # (3) rank correlation of the sequence means
    rho_d <- vapply(seq_len(Tn), function(t) cor(a$mu_s[t, ], b$mu_s[t, ], method = "spearman"), numeric(1))
    q <- quantile(rho_d, c(0.025, 0.975), names = FALSE)
    spm <- data.table(rq = sp$rq, treatment_fn = fn, treatment_fp = fp, model = sp$model, n_draws = Tn,
                      median = median(rho_d), mean = mean(rho_d), q025 = q[1], q975 = q[2],
                      p_gt0 = mean(rho_d > 0))
    f3 <- file.path(path_csv, paste0("ex3-", sp$rq, "-spearman.csv"))
    fwrite(spm, f3);  msg("Saved: ", f3)
  }

  invisible(TRUE)
}
