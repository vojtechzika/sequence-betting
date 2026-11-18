# ============================================================
# 206_mechanism-checks-table.R
#
# PURPOSE
#   Table of the mechanism checks of the Synthesis for the online appendix:
#   (A) whether selection into betting can explain the relation between the
#   betting-rate deviation and the side bias of the sequences, (B) the
#   conditional test of selection given the imbalance class and the terminal
#   run, (C) the bimodality of the betting-rate deviation.
#
# SAMPLE
#   FN, first treatment of cfg$run$treatment, confirmatory sample (the
#   Heads share comes from the full RQ4 fit, as set in 91_synthesis-sequences.R).
#
# INPUT (data/output/csv/)
#   synthesis-mechanism-checks-summary-<treatment>-confirmatory.csv   (93_synthesis-mechanism-checks.R)
#   synthesis-mechanism-checks-strata-<treatment>-confirmatory.csv
#
# OUTPUT (data/output/tables/)
#   synthesis-mechanism-checks.tex   table with panels A to C
#
# NOTES
#   - Strata with fewer sequences than the minimum of 93_synthesis-mechanism-checks.R
#     have no correlation and are not listed
#   - Output files are overwritten at every run
# ============================================================

mechanism_checks_table <- function(cfg) {

  tr <- unique(as.character(cfg$run$treatment))[1L]
  s  <- fread(file.path(path_csv, paste0("synthesis-mechanism-checks-summary-", tr, "-confirmatory.csv")), encoding = "UTF-8")
  st <- fread(file.path(path_csv, paste0("synthesis-mechanism-checks-strata-", tr, "-confirmatory.csv")), encoding = "UTF-8")

  num <- function(x) paste0("$", gsub("-", "-", sprintf("%.3f", x)), "$")
  pf  <- function(p) if (is.na(p)) "--" else if (p < 0.001) "$<.001$" else paste0("$", sub("^0", "", sprintf("%.3f", p)), "$")
  get <- function(chk) { r <- s[check == chk]; stopifnot(nrow(r) == 1L); r }
  line <- function(label, chk, stat) {
    r <- get(chk)
    paste0("\t\t", label, " & ", stat, " & ", num(r$value), " & ", pf(r$p_value), " & ", r$n, " \\\\")
  }

  A <- c(
    line("CI width vs.\\ betting-rate deviation (all sequences)",            "ciw_vs_b (all 64)",                          "Spearman $\\rho$"),
    line("Bias magnitude vs.\\ betting-rate deviation (all sequences)",      "habs_vs_b (all 64, unconditional)",          "Spearman $\\rho$"),
    line("CI width vs.\\ bias magnitude (all sequences)",                    "ciw_vs_habs (all 64)",                       "Spearman $\\rho$"),
    line("CI width, below- vs.\\ above-mean betting rate",                   "ciw low-bet vs high-bet (t-test diff)",      "Mean diff."),
    line("Bias magnitude, below- vs.\\ above-mean betting rate",             "habs low-bet vs high-bet (t-test diff)",     "Mean diff."),
    line("CI width, lower- vs.\\ higher-betting pair member (paired)",       "pair: ciw_lo vs ciw_hi (wilcoxon)",          "Wilcoxon $V$"),
    line("Bias magnitude, lower- vs.\\ higher-betting pair member (paired)", "pair: habs_lo vs habs_hi (wilcoxon)",        "Wilcoxon $V$"),
    line("Within-pair betting-rate gap vs.\\ CI-width gap",                  "pair: b_gap vs ciw_gap",                     "Spearman $\\rho$"),
    line("Bias magnitude vs.\\ betting-rate deviation, credible shifts only","habs_vs_b, b_sig=TRUE only",                 "Spearman $\\rho$")
  )
  # Wilcoxon V is an integer
  A <- sub("\\$(\\d+)\\.000\\$", "$\\1$", A)

  B <- c(
    paste0("\t\tMixture components selected (1 vs.\\ 2) & $G$ & $", get("mixture: chosen G (mclust)")$value, "$ & -- & 64 \\\\"),
    paste0("\t\tBIC gap (2-component $-$ 1-component) & $\\Delta$BIC & $", sprintf("%.2f", get("mixture: BIC gap (2comp - 1comp, higher favors 2comp)")$value), "$ & -- & 64 \\\\"),
    paste0("\t\tBimodality coefficient & $BC$ & $", sprintf("%.3f", get("bimodality coefficient")$value), "$ & -- & 64 \\\\")
  )

  b  <- get("CONDITIONAL: habs ~ b_d + imb_class + run, coef on b_d (HC1)")
  lo <- get("CONDITIONAL: coef on b_d, 95% CI lower")$value
  hi <- get("CONDITIONAL: coef on b_d, 95% CI upper")$value
  f2 <- function(x) gsub("-", "$-$", sprintf("%.2f", x))
  C  <- paste0("\t\tBias magnitude on betting-rate deviation, given imbalance class and terminal run & OLS $\\beta$ [95\\% CI] & $",
               sprintf("%.2f", b$value), "$ [", f2(lo), ", ", f2(hi), "] & ", pf(b$p_value), " & 64 \\\\")
  st <- st[!is.na(rho)]
  C <- c(C, vapply(seq_len(nrow(st)), function(i) {
    paste0("\t\t\\quad Within imbalance class ", st$imb_class[i], ", ",
           if (st$run[i] == 1L) "with" else "without", " terminal run & Spearman $\\rho$ & ", num(st$rho[i]), " & ", pf(st$p_value[i]), " & ", st$n[i], " \\\\")
  }, ""))

  out <- c(
    "\\begin{table}[htbp]", "\t\\centering", "\t\\footnotesize",
    "\t\\begin{tabular}{>{\\raggedright\\arraybackslash}p{0.44\\textwidth}cccc}", "\t\t\\toprule",
    "\t\tCheck & Statistic & Value & $p$-value & $n$ \\\\", "\t\t\\midrule",
    "\t\t\\multicolumn{5}{l}{\\textit{Panel A: Selection into betting}} \\\\", "\t\t\\addlinespace", A,
    "\t\t\\addlinespace", "\t\t\\multicolumn{5}{l}{\\textit{Panel B: Conditional test of selection}} \\\\", "\t\t\\addlinespace", C,
    "\t\t\\addlinespace", "\t\t\\multicolumn{5}{l}{\\textit{Panel C: Bimodality of the betting-rate deviation}} \\\\", "\t\t\\addlinespace", B,
    "\t\t\\bottomrule", "\t\\end{tabular}",
    "\t\\caption{",
    "\t\t\\footnotesize",
    "\t\t\\textbf{Mechanism checks: side bias, betting-rate deviation, and bimodality.}",
    "\t\tCI width is the width of the 95\\% credible interval of the Heads-share deviation of a sequence, and bias magnitude is the absolute value of that deviation. The betting-rate deviation is the difference of the betting rate of a sequence from its grand mean. Mean-difference rows compare sequences with a below- and an above-mean betting rate ($t$-test). Paired rows compare the lower- and the higher-betting member of each of the 32 mirror pairs (Wilcoxon signed-rank test). The mixture model is fit with \\texttt{mclust}, which reports the BIC on a higher-is-better scale, so that a positive gap favors two components; the bimodality coefficient is compared with the heuristic threshold of $0.555$ and is not a formal test. In Panel B, the terminal run marks sequences whose last two tosses are identical, the regression uses heteroskedasticity-robust standard errors (HC1), and strata with too few sequences are omitted. The confirmatory FN sample is used for the betting rate, and the full sample for the Heads share.",
    "\t}",
    "\t\\label{Tab:synthesis-mechanism-checks}", "\\end{table}")

  f <- file.path(path_tab, "synthesis-mechanism-checks.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
