# ============================================================
# 91_synthesis-sequences.R
#
# PURPOSE
#   Synthesis, sequence level. For every sequence, each of the four headline
#   estimands is expressed as a paired posterior difference from its own
#   grand mean, d_s(t) = mu_s(t) - grand(t), read off the same posterior
#   draw t of the same fit (no independence assumption and no normal
#   approximation). Per sequence and estimand the script reports the
#   posterior median of d_s, the equal-tailed credible interval and
#   sig = TRUE if the interval excludes zero (a plain zero-exclusion rule,
#   not the practical-equivalence band used for the labels of the RQ
#   tables). The output feeds the mirror-pair figure, the mechanism checks
#   and the sequence tables of the Synthesis section.
#
#   Estimands (generated quantities of the selected models)
#     h  Heads share, RQ4        mu_h,  grand mean hbar
#     c  welfare loss, RQ3       mu_c,  grand mean = mean of mu_c over sequences
#     b  betting rate, RQ1       mu_b,  grand mean grand_mean_bet
#     a  stake deviation, RQ2    mu_a,  grand mean = mean of mu_a over sequences
#
# SAMPLES (runs)
#   confirmatory   first treatment of cfg$run$treatment, if it has a
#                  confirmatory sample: RQ1-RQ3 from the confirmatory fit,
#                  RQ4 from the full fit (RQ4 has no confirmatory sample)
#   full           every treatment in cfg$run$treatment: all four from the
#                  full-sample fits
#   The models are those selected under the keys rq1, rq2, rq3 and rq4.
#
# INPUT
#   models "rq1".."rq4": fit-, seq-levels-<tr>-<sample>.rds
#   cfg$design$a_flags$betting_normative (which treatments have a
#   confirmatory sample)
#
# OUTPUT (data/output/csv/)
#   synthesis-<tr>-<sample>-sequence-diffs.csv   one row per sequence:
#       <e>_d_median, <e>_d_lo, <e>_d_hi, <e>_sig for e in h, c, b, a
#   synthesis-<tr>-<sample>-grand-means.csv       per estimand: model folder,
#       grand median, interval, mean of the sequence medians
#
# NOTES
#   - The all-sequence run is required by the mirror-pair figure and the
#     mechanism checks, which need all 32 mirror pairs
#   - A warning is issued if the grand median and the mean of the sequence
#     medians differ by more than 0.01 for an estimand
#   - A missing fit for one estimand is skipped with a warning; if no fit
#     can be loaded the script stops
#   - Output files are overwritten at every run
# ============================================================

synthesis_sequences <- function(cfg, ci = 0.95) {

  design <- cfg$design
  tr_all <- unique(as.character(cfg$run$treatment))
  fn     <- tr_all[1L]

  runs <- list()
  if (isTRUE(design$a_flags$betting_normative[[fn]])) runs[[1L]] <- c(fn, "confirmatory")
  for (tr in tr_all) runs[[length(runs) + 1L]] <- c(tr, "full")

  # Paired difference of each sequence mean from the grand mean, by draw
  diff_from_grand <- function(mu, grand, seq_levels) {
    stopifnot(is.matrix(mu), length(grand) == nrow(mu))
    d   <- mu - grand
    lo  <- apply(d, 2L, quantile, probs = (1 - ci) / 2)
    hi  <- apply(d, 2L, quantile, probs = 1 - (1 - ci) / 2)
    data.table(sequence = seq_levels, d_median = apply(d, 2L, median),
               d_lo = lo, d_hi = hi, sig = (lo > 0) | (hi < 0))
  }

  # Draws of an estimand from the model selected under `key`
  load_draws <- function(key, tr, sample, mu_par, grand_par = NULL) {
    f <- function(stem) model_path(key, paste0(stem, "-", tr, "-", sample, ".rds"))
    if (!all(file.exists(f("fit"), f("seq-levels")))) return(NULL)
    post <- rstan::extract(readRDS(f("fit")), pars = c(mu_par, grand_par))
    stopifnot(is.matrix(post[[mu_par]]))
    grand <- if (is.null(grand_par)) rowMeans(post[[mu_par]]) else as.numeric(post[[grand_par]])
    list(mu = post[[mu_par]], grand = grand,
         seq_levels = as.character(readRDS(f("seq-levels"))),
         model = basename(model_path(key)))
  }

  for (r in runs) {
    tr <- r[1L]; tag <- r[2L]
    rq_tag <- tag
    specs <- list(
      h = load_draws("rq4", tr, "full", "mu_h", "hbar"),
      c = load_draws("rq3", tr, rq_tag, "mu_c"),
      b = load_draws("rq1", tr, rq_tag, "mu_b", "grand_mean_bet"),
      a = load_draws("rq2", tr, rq_tag, "mu_a")
    )
    for (nm in names(specs))
      if (is.null(specs[[nm]])) warning("Synthesis sequences: no fit for estimand '", nm,
                                        "' (", tr, "/", tag, "); skipped.", call. = FALSE)
    loaded <- specs[!vapply(specs, is.null, logical(1L))]
    if (!length(loaded)) stop("Synthesis sequences: no fits for ", tr, "/", tag, call. = FALSE)

    sets      <- lapply(loaded, `[[`, "seq_levels")
    sequences <- sets[[1L]][sets[[1L]] %in% Reduce(intersect, sets)]
    if (!length(sequences)) stop("No sequence common to all estimands (", tr, "/", tag, ").", call. = FALSE)

    out <- NULL; grand_rows <- list()
    for (nm in names(loaded)) {
      sp  <- loaded[[nm]]
      idx <- match(sequences, sp$seq_levels)
      tbl <- diff_from_grand(sp$mu[, idx, drop = FALSE], sp$grand, sequences)
      setnames(tbl, c("d_median", "d_lo", "d_hi", "sig"), paste0(nm, c("_d_median", "_d_lo", "_d_hi", "_sig")))
      out <- if (is.null(out)) tbl else merge(out, tbl, by = "sequence", sort = FALSE)
      grand_rows[[nm]] <- data.table(
        estimand = nm, model = sp$model,
        grand_median = median(sp$grand),
        grand_lo = unname(quantile(sp$grand, (1 - ci) / 2)),
        grand_hi = unname(quantile(sp$grand, 1 - (1 - ci) / 2)),
        grand_mean_of_sequence_medians = mean(apply(sp$mu, 2L, median)))
    }
    out <- out[match(sequences, sequence)]
    grand_tbl <- rbindlist(grand_rows)

    gap <- abs(grand_tbl$grand_median - grand_tbl$grand_mean_of_sequence_medians)
    if (any(gap > 0.01))
      warning("Synthesis sequences (", tr, "/", tag, "): grand median and mean of the sequence medians differ by > 0.01 for: ",
              paste(grand_tbl$estimand[gap > 0.01], collapse = ", "), call. = FALSE)

    f1 <- file.path(path_csv, paste0("synthesis-", tr, "-", tag, "-sequence-diffs.csv"))
    f2 <- file.path(path_csv, paste0("synthesis-", tr, "-", tag, "-grand-means.csv"))
    fwrite(out, f1);        msg("Saved: ", f1)
    fwrite(grand_tbl, f2);  msg("Saved: ", f2)
  }

  invisible(TRUE)
}
