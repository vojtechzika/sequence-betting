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
#   For the sequences that the reference model places above or below the grand
#   mean with at least 95% posterior probability (the rule of 31_rq3-tables.R),
#   the script also reports how many keep the same side and the same class under
#   the comparator, how many pure sequences (one symbol repeated) are classed as
#   above the grand mean by each model, and how many of the five sequences with
#   the highest median loss the two models share. The estimates and classes of
#   every sequence are written to a second file.
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
#   data/output/csv/rq3-likelihood-comparison-sequences.csv
#       one row per treatment, sample, comparator and sequence: median and 95%
#       interval of the sequence mean under the reference and the comparator,
#       probability of exceeding the grand mean, and the class of each model
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

  # Probability that a sequence mean exceeds / falls below the grand mean of the
  # same draw, and the class derived from it (thresholds of 31_rq3-tables.R)
  p_above <- function(mu) {
    g <- rowMeans(mu)
    apply(mu, 2, function(x) mean(x > g))
  }
  p_below <- function(mu) {
    g <- rowMeans(mu)
    apply(mu, 2, function(x) mean(x < g))
  }
  classify <- function(mu) {
    pa <- p_above(mu)
    pb <- p_below(mu)
    fifelse(pa >= 0.95, "above",
    fifelse(pa >= 0.80, "likely_above",
    fifelse(pb >= 0.95, "below",
    fifelse(pb >= 0.80, "likely_below", "neutral"))))
  }
  top5 <- function(m) order(m, decreasing = TRUE)[1:5]

  rows     <- list()
  seq_rows <- list()
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

      # Sequence-level classes under the reference and the comparator
      seq_chr <- as.character(ref$seq)
      is_pure <- !grepl("H", seq_chr) | !grepl("T", seq_chr)
      lab_ref <- classify(ref$mu)
      lab_alt <- classify(alt$mu)
      pa_ref  <- p_above(ref$mu)
      pa_alt  <- p_above(alt$mu)
      med_ref <- apply(ref$mu, 2, median)
      med_alt <- apply(alt$mu, 2, median)
      flagged <- lab_ref %in% c("above", "below")

      seq_rows[[length(seq_rows) + 1L]] <- data.table(
        treatment  = tr,
        sample     = smp,
        comparator = cmp,
        sequence   = seq_chr,
        pure       = is_pure,
        ref_median = med_ref,
        ref_q025   = apply(ref$mu, 2, quantile, probs = 0.025),
        ref_q975   = apply(ref$mu, 2, quantile, probs = 0.975),
        comp_median = med_alt,
        comp_q025  = apply(alt$mu, 2, quantile, probs = 0.025),
        comp_q975  = apply(alt$mu, 2, quantile, probs = 0.975),
        ref_p_above_grand  = pa_ref,
        comp_p_above_grand = pa_alt,
        ref_class  = lab_ref,
        comp_class = lab_alt
      )

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
        max_abs_diff_seq  = max(abs(m_ref - m_alt)),
        n_flagged_ref        = sum(flagged),
        n_flagged_same_side  = sum(flagged & (pa_ref >= 0.5) == (pa_alt >= 0.5)),
        n_flagged_same_class = sum(flagged & lab_alt == lab_ref),
        n_pure_above_ref     = sum(is_pure & lab_ref == "above"),
        n_pure_above_comp    = sum(is_pure & lab_alt == "above"),
        top5_shared          = length(intersect(top5(med_ref), top5(med_alt)))
      )
    }
  }

  stopifnot(length(rows) > 0L)
  out <- rbindlist(rows)

  f <- file.path(path_csv, "rq3-likelihood-comparison.csv")
  fwrite(out, f)
  msg("Saved: ", f)

  f_seq_out <- file.path(path_csv, "rq3-likelihood-comparison-sequences.csv")
  fwrite(rbindlist(seq_rows), f_seq_out)
  msg("Saved: ", f_seq_out)

  invisible(out)
}
