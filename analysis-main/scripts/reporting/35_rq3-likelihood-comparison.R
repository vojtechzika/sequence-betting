# ============================================================
# 35_rq3-likelihood-comparison.R
#
# PURPOSE
#   Compares the RQ3 sequence-level welfare loss across the three likelihoods
#   fitted in the modelling stage: the hurdle-Gamma model (the preregistered
#   primary model, used as the reference), the Gamma-only model and the
#   Gaussian model. For every treatment and sample and for each comparator
#   the script reports
#     - the grand mean welfare loss of the reference and of the comparator
#       (posterior median and 95% interval),
#     - the Spearman correlation of the sequence means,
#     - the share of sequence pairs ordered in the same direction,
#     - the share of sequences on the same side of the grand mean,
#     - the largest absolute difference between the sequence means,
#     - the model-implied mean probability of a zero loss in the sample.
#
# INTERPRETATION
#   The sequence mean mu_c is the unconditional mean loss, zeros included, in
#   the hurdle-Gamma and Gaussian models. In the Gamma-only model it is the
#   mean of the positive losses only, because the model is fitted to the
#   trials with a positive loss. Levels of the Gamma-only model are therefore
#   comparable with the hurdle-Gamma model only where zero losses are rare;
#   in other samples (see zero_share_model) only the ordering statistics are
#   informative. The Gaussian model ignores the boundary at zero and the
#   right skew, so it serves as a check on the sign and ordering of the
#   sequence effects.
#
# FUNCTION
#   rq3_likelihood_comparison(cfg)
#
# INPUT (model folders are read by name, independent of selected-models.csv)
#   rq3-hurdle-gamma, rq3-gamma, rq3-gaussian:
#     fit-<tr>-<sample>.rds, seq-levels-<tr>-<sample>.rds  (draws of mu_c)
#   rq3-hurdle-gamma: hurdle-summary.csv
#
# OUTPUT
#   data/output/csv/rq3-likelihood-comparison.csv
#       one row per treatment, sample and comparator
#
# NOTES
#   - Samples without fits in all three folders are skipped
#   - Output files are overwritten at every run
# ============================================================

rq3_likelihood_comparison <- function(cfg) {

  tr_vec  <- unique(as.character(cfg$run$treatment))
  samples <- c("full", "confirmatory")

  ref_folder  <- "rq3-hurdle-gamma"
  comparators <- c(gamma_only = "rq3-gamma", gaussian = "rq3-gaussian")

  hs <- fread(model_path(ref_folder, "hurdle-summary.csv"))

  # Draws of the sequence means mu_c (draws x sequences) and sequence labels
  get_fit <- function(folder, tr, smp) {
    f_fit <- model_path(folder, paste0("fit-", tr, "-", smp, ".rds"), must_exist = FALSE)
    f_seq <- model_path(folder, paste0("seq-levels-", tr, "-", smp, ".rds"), must_exist = FALSE)
    if (!file.exists(f_fit) || !file.exists(f_seq)) return(NULL)
    list(mu  = rstan::extract(readRDS(f_fit), pars = "mu_c")$mu_c,
         seq = readRDS(f_seq))
  }

  qs <- function(x) quantile(x, c(0.025, 0.975), names = FALSE)
  # direction (-1, 0, 1) of every pair of sequences
  pair_dir <- function(m) {
    d <- sign(outer(m, m, "-"))
    d[upper.tri(d)]
  }

  rows <- list()
  for (tr in tr_vec) for (smp in samples) {

    ref <- get_fit(ref_folder, tr, smp)
    if (is.null(ref)) next

    for (cmp in names(comparators)) {
      alt <- get_fit(comparators[[cmp]], tr, smp)
      if (is.null(alt)) next
      stopifnot(identical(as.character(ref$seq), as.character(alt$seq)))

      m_ref <- colMeans(ref$mu)
      m_alt <- colMeans(alt$mu)
      g_ref <- rowMeans(ref$mu)
      g_alt <- rowMeans(alt$mu)

      rows[[length(rows) + 1L]] <- data.table(
        treatment         = tr,
        sample            = smp,
        comparator        = cmp,
        n_sequences       = length(m_ref),
        zero_share_model  = hs[treatment == tr & sample == smp, pi0_mean],
        grand_ref_median  = median(g_ref),
        grand_ref_q025    = qs(g_ref)[1],
        grand_ref_q975    = qs(g_ref)[2],
        grand_comp_median = median(g_alt),
        grand_comp_q025   = qs(g_alt)[1],
        grand_comp_q975   = qs(g_alt)[2],
        spearman_seq      = cor(m_ref, m_alt, method = "spearman"),
        pair_agreement    = mean(pair_dir(m_ref) == pair_dir(m_alt)),
        side_agreement    = mean(sign(m_ref - mean(m_ref)) == sign(m_alt - mean(m_alt))),
        max_abs_diff_seq  = max(abs(m_ref - m_alt))
      )
    }
  }

  stopifnot(length(rows) > 0L)
  out <- rbindlist(rows)

  f <- file.path(path_csv, "rq3-likelihood-comparison.csv")
  fwrite(out, f)
  msg("Saved: ", f)

  invisible(out)
}
