# ============================================================
# 93_synthesis_mechanism_checks.R
#
# PURPOSE
#   Checks two things about the sequence-level betting-rate / side-choice
#   pattern surfaced by 92_synthesis_mirror_plot.R:
#
#   (1) MECHANICAL-EXPLANATION (SELECTION) CHECK. Side choice is only
#       observed among participants who bet, so sequences with lower
#       betting rates have their heads-share estimated on a smaller,
#       self-selected sample. If selective bettors are also more
#       directional, low-betting sequences would show inflated side bias.
#       (a) does heads-share CI width relate to betting rate? (expected:
#           yes, mechanically -- fewer bettors = more posterior
#           uncertainty); (b) does heads-share MAGNITUDE (|h_d_median|)
#           relate to betting rate?; (c) within mirror pairs: does the
#           lower-betting member show a wider CI / larger |h|?
#       (d) NEW -- CONDITIONAL test. Imbalance raises betting AND lowers
#           side bias, and the terminal run raises side bias without
#           moving betting, so the unconditional correlations in (b)
#           cannot rule selection in or out: sequence structure confounds
#           them. Selection predicts a NEGATIVE association between |h|
#           and the betting delta once imbalance class and the run are
#           held fixed. Tested by OLS (HC1) of |h| on b_d + imbalance
#           class + run, and by Spearman rho within each
#           imbalance-class x run stratum with >= min_stratum sequences.
#           (The participant-level analogue -- the joint model's
#           correlation between betting propensity and side-choice slopes
#           -- is exported by 95_synthesis_stan_outputs.R.)
#
#   (2) BIMODALITY CHECK. Does the betting-rate delta (b_d_median) look
#       like it's drawn from two distinct sequence "types"? Checked via a
#       1- vs 2-component Gaussian mixture (BIC, mclust::Mclust) and a
#       bimodality coefficient, plus |h| vs b_d among the credibly
#       significant (b_sig = TRUE) sequences.
#
#   Reads only the already-computed CSV from 91_synthesis_sequences.R;
#   does not touch Stan fits.
#
# INPUT
#   path_out/synthesis_sequence_diffs_<tr>_<tag>.csv (ALL-SEQUENCES run)
#
# OUTPUT
#   path_out/synthesis_mechanism_checks_summary_<tr>_<tag>.csv -- one row
#     per scalar test (now including the conditional tests in (1d))
#   path_out/synthesis_mechanism_checks_pairs_<tr>_<tag>.csv -- 32 pairs
#   path_out/synthesis_mechanism_checks_strata_<tr>_<tag>.csv -- (1d)
#     within-stratum correlations
#
# NOTES
#   - Requires mclust (mixture check) and sandwich + lmtest (HC1 SEs).
#   - mclust's BIC is HIGHER-is-better.
#   - Skewness/kurtosis use the biased (population) convention to match
#     an earlier scipy cross-check; BC > 0.555 is a heuristic threshold.
#   - The pure sequences (0-6) form a 2-sequence stratum and never enter
#     the within-stratum correlations; they do enter the OLS.
# ============================================================
synthesis_mechanism_checks <- function(cfg, tr = "m25", tag = "confirmatory",
                                       quiet = TRUE, min_stratum = 6L) {
  
  for (pkg in c("mclust", "sandwich", "lmtest")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("synthesis_mechanism_checks: package '", pkg, "' is required. ",
           "Install it with install.packages(\"", pkg, "\") and re-run.")
    }
  }
  # Mclust() calls mclustBIC() unqualified, so mclust must be attached.
  suppressPackageStartupMessages(library(mclust))
  
  mirror_seq <- function(s) chartr("HT", "TH", s)
  
  biased_skewness <- function(x) {
    m <- mean(x); m2 <- mean((x - m)^2); m3 <- mean((x - m)^3); m3 / m2^1.5
  }
  biased_excess_kurtosis <- function(x) {
    m <- mean(x); m2 <- mean((x - m)^2); m4 <- mean((x - m)^4); m4 / m2^2 - 3
  }
  bimodality_coefficient <- function(x) {
    n <- length(x)
    (biased_skewness(x)^2 + 1) /
      (biased_excess_kurtosis(x) + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
  }
  
  # ---- Main body ----
  f_in <- file.path(path_out, paste0("synthesis_sequence_diffs_", tr, "_", tag, ".csv"))
  if (!file.exists(f_in)) {
    stop("synthesis_mechanism_checks: expected input not found: ", f_in,
         " -- run synthesis_sequences(cfg) first (91_synthesis_sequences.R).")
  }
  dt <- fread(f_in)
  
  required_cols <- c("sequence", "h_d_median", "h_d_lo", "h_d_hi", "h_sig",
                     "b_d_median", "b_sig")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("synthesis_mechanism_checks: input CSV is missing column(s): ",
         paste(missing_cols, collapse = ", "), ".")
  }
  
  dt <- copy(dt)
  dt[, h_ci_width := h_d_hi - h_d_lo]
  dt[, h_abs      := abs(h_d_median)]
  
  # sequence structure (outcome-blind, from the string itself)
  n_H <- nchar(gsub("T", "", dt$sequence))
  dt[, imb_class := factor(abs(n_H - 3), levels = 0:3, labels = c("3-3", "2-4", "1-5", "0-6"))]
  dt[, run       := as.integer(substr(sequence, 5, 5) == substr(sequence, 6, 6))]
  
  # ---- (1a)/(1b): aggregate (unconditional) checks ----
  sp <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  cor_ciw_b    <- sp(dt$h_ci_width, dt$b_d_median)
  cor_habs_b   <- sp(dt$h_abs,      dt$b_d_median)
  cor_ciw_habs <- sp(dt$h_ci_width, dt$h_abs)
  
  low  <- dt[b_d_median < 0]
  high <- dt[b_d_median > 0]
  tt_ciw  <- t.test(low$h_ci_width, high$h_ci_width)
  tt_habs <- t.test(low$h_abs,      high$h_abs)
  
  aggregate <- list(
    cor_ciw_vs_b = cor_ciw_b, cor_habs_vs_b = cor_habs_b, cor_ciw_vs_habs = cor_ciw_habs,
    ttest_ciw_low_vs_high_bet = tt_ciw, ttest_habs_low_vs_high_bet = tt_habs,
    mean_ciw_low_bet  = mean(low$h_ci_width), mean_ciw_high_bet  = mean(high$h_ci_width),
    mean_habs_low_bet = mean(low$h_abs),      mean_habs_high_bet = mean(high$h_abs),
    n_low_bet = nrow(low), n_high_bet = nrow(high)
  )
  
  # ---- (1c): within-mirror-pair asymmetry ----
  d2 <- copy(dt)
  d2[, mirror := mirror_seq(sequence)]
  missing_mirror <- setdiff(d2$mirror, d2$sequence)
  if (length(missing_mirror) > 0L && !quiet) {
    warning("synthesis_mechanism_checks: mirror sequence(s) missing, pair(s) excluded: ",
            paste(unique(missing_mirror), collapse = ", "))
  }
  d2 <- d2[mirror %in% sequence]
  d2[, pair_key := paste(pmin(sequence, mirror), pmax(sequence, mirror))]
  
  pairs <- rbindlist(lapply(unique(d2$pair_key), function(pk) {
    two <- d2[pair_key == pk]
    if (nrow(two) != 2L) return(NULL)
    ord <- order(two$b_d_median)
    lo <- two[ord[1]]; hi <- two[ord[2]]
    data.table(
      sequence_lo_bet = lo$sequence, sequence_hi_bet = hi$sequence,
      b_lo_bet = lo$b_d_median, b_hi_bet = hi$b_d_median,
      b_gap = hi$b_d_median - lo$b_d_median,
      h_abs_lo_bet = lo$h_abs, h_abs_hi_bet = hi$h_abs,
      ciw_lo_bet = lo$h_ci_width, ciw_hi_bet = hi$h_ci_width
    )
  }))
  pairs[, ciw_gap  := ciw_lo_bet - ciw_hi_bet]
  pairs[, habs_gap := h_abs_lo_bet - h_abs_hi_bet]
  
  wilcox_ciw      <- suppressWarnings(wilcox.test(pairs$ciw_lo_bet, pairs$ciw_hi_bet, paired = TRUE))
  wilcox_habs     <- suppressWarnings(wilcox.test(pairs$h_abs_lo_bet, pairs$h_abs_hi_bet, paired = TRUE))
  cor_bgap_ciwgap <- sp(pairs$b_gap, pairs$ciw_gap)
  
  pair_tests <- list(
    wilcox_ciw_lo_vs_hi_bet = wilcox_ciw, wilcox_habs_lo_vs_hi_bet = wilcox_habs,
    cor_bgap_vs_ciwgap = cor_bgap_ciwgap,
    frac_lo_bet_wider_ci = mean(pairs$ciw_gap > 0),
    frac_lo_bet_larger_habs = mean(pairs$habs_gap > 0),
    n_pairs = nrow(pairs)
  )
  
  # ---- (1d): CONDITIONAL selection test (see header) ----
  fit_cond <- lm(h_abs ~ b_d_median + imb_class + run, data = dt)
  ct_cond  <- lmtest::coeftest(fit_cond, vcov. = sandwich::vcovHC(fit_cond, type = "HC1"))
  ci_cond  <- lmtest::coefci(fit_cond, vcov. = sandwich::vcovHC(fit_cond, type = "HC1"))
  
  strata <- dt[, {
    if (.N >= min_stratum) {
      r <- sp(h_abs, b_d_median)
      .(rho = unname(r$estimate), p_value = r$p.value, n = .N)
    } else .(rho = NA_real_, p_value = NA_real_, n = .N)
  }, by = .(imb_class, run)]
  setorder(strata, imb_class, run)
  
  # ---- (2): bimodality of b_d_median ----
  mix1 <- mclust::Mclust(dt$b_d_median, G = 1, verbose = FALSE)
  mix2 <- mclust::Mclust(dt$b_d_median, G = 2, verbose = FALSE)
  mix  <- if (mix2$bic > mix1$bic) mix2 else mix1
  bic_gap    <- mix2$bic - mix1$bic
  bimod_coef <- bimodality_coefficient(dt$b_d_median)
  
  sig <- dt[b_sig == TRUE]
  sig_split <- list(n_total = nrow(dt), n_sig = nrow(sig),
                    n_sig_positive = sum(sig$b_d_median > 0),
                    n_sig_negative = sum(sig$b_d_median < 0))
  cor_sig_only <- if (nrow(sig) >= 3L) sp(sig$h_abs, sig$b_d_median) else NULL
  
  row <- function(check, statistic, value, p_value, n)
    data.table(check = check, statistic = statistic, value = value, p_value = p_value, n = n)
  
  summary_tbl <- rbindlist(list(
    row("ciw_vs_b (all 64)", "spearman_rho", unname(cor_ciw_b$estimate), cor_ciw_b$p.value, nrow(dt)),
    row("habs_vs_b (all 64, unconditional)", "spearman_rho", unname(cor_habs_b$estimate), cor_habs_b$p.value, nrow(dt)),
    row("ciw_vs_habs (all 64)", "spearman_rho", unname(cor_ciw_habs$estimate), cor_ciw_habs$p.value, nrow(dt)),
    row("ciw low-bet vs high-bet (t-test diff)", "mean_diff",
        mean(low$h_ci_width) - mean(high$h_ci_width), tt_ciw$p.value, nrow(dt)),
    row("habs low-bet vs high-bet (t-test diff)", "mean_diff",
        mean(low$h_abs) - mean(high$h_abs), tt_habs$p.value, nrow(dt)),
    row("pair: ciw_lo vs ciw_hi (wilcoxon)", "V", unname(wilcox_ciw$statistic), wilcox_ciw$p.value, nrow(pairs)),
    row("pair: habs_lo vs habs_hi (wilcoxon)", "V", unname(wilcox_habs$statistic), wilcox_habs$p.value, nrow(pairs)),
    row("pair: b_gap vs ciw_gap", "spearman_rho", unname(cor_bgap_ciwgap$estimate), cor_bgap_ciwgap$p.value, nrow(pairs)),
    row("CONDITIONAL: habs ~ b_d + imb_class + run, coef on b_d (HC1)", "beta",
        ct_cond["b_d_median", 1], ct_cond["b_d_median", 4], nrow(dt)),
    row("CONDITIONAL: coef on b_d, 95% CI lower", "ci_low", ci_cond["b_d_median", 1], NA_real_, nrow(dt)),
    row("CONDITIONAL: coef on b_d, 95% CI upper", "ci_high", ci_cond["b_d_median", 2], NA_real_, nrow(dt)),
    strata[!is.na(rho), .(check = paste0("CONDITIONAL: habs_vs_b within ", imb_class, ", run=", run),
                          statistic = "spearman_rho", value = rho, p_value = p_value, n = n)],
    row("mixture: chosen G (mclust)", "G", mix$G, NA_real_, nrow(dt)),
    row("mixture: BIC gap (2comp - 1comp, higher favors 2comp)", "BIC_diff", bic_gap, NA_real_, nrow(dt)),
    row("bimodality coefficient", "BC", bimod_coef, NA_real_, nrow(dt)),
    row("habs_vs_b, b_sig=TRUE only", "spearman_rho",
        if (!is.null(cor_sig_only)) unname(cor_sig_only$estimate) else NA_real_,
        if (!is.null(cor_sig_only)) cor_sig_only$p.value else NA_real_, nrow(sig))
  ))
  
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  save <- function(x, stem, label) {
    f <- file.path(path_out, paste0(stem, tr, "_", tag, ".csv"))
    if (!should_skip(f, cfg, "output", paste0(label, " (", tr, "/", tag, ")"))) {
      fwrite(x, f); msg("Saved: ", f)
    }
  }
  save(summary_tbl, "synthesis_mechanism_checks_summary_", "Synthesis mechanism checks summary")
  save(pairs,       "synthesis_mechanism_checks_pairs_",   "Synthesis mechanism checks pairs")
  save(strata,      "synthesis_mechanism_checks_strata_",  "Synthesis mechanism checks strata")
  
  invisible(list(aggregate = aggregate, pairs = pairs, pair_tests = pair_tests,
                 conditional = list(fit = ct_cond, ci = ci_cond, strata = strata),
                 mixture = mix, bimodality_coef = bimod_coef, sig_split = sig_split,
                 cor_sig_only = cor_sig_only, summary = summary_tbl))
}