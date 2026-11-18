# ============================================================
# 43_rq4-top-sequences-table.R
#
# PURPOSE
#   Table of the five sequences with the highest and the five with the lowest
#   posterior mean Heads probability (RQ4), side by side, with the 95%
#   credible interval of each.
#
# SAMPLE
#   Full sample of the first treatment in cfg$run$treatment (FN).
#
# INPUT (data/output/csv/)
#   rq4-<treatment>-full-sequences.csv        mu_h_mean, mu_h_q025, mu_h_q975
#
# OUTPUT (data/output/tables/)
#   rq4-top-sequences.tex
#
# NOTES
#   - Rows are paired by rank: the highest Heads probability next to the lowest
#   - The number of sequences per side is top_n
#   - Output file is overwritten at every run
# ============================================================

rq4_top_sequences_table <- function(cfg, top_n = 5L) {

  tr <- unique(as.character(cfg$run$treatment))[1L]

  pct <- function(m, lo, hi) {
    p <- function(x) sprintf("%.1f\\%%", 100 * x)
    paste0(p(m), " [", p(lo), ", ", p(hi), "]")
  }

  d <- fread(file.path(path_csv, paste0("rq4-", tr, "-full-sequences.csv")), encoding = "UTF-8")
  d <- d[order(-mu_h_mean)]
  stopifnot(nrow(d) == 64L)
  hd <- d[seq_len(top_n)]                       # highest Heads probability
  tl <- d[order(mu_h_mean)][seq_len(top_n)]     # lowest Heads probability

  cell <- function(x) pct(x$mu_h_mean, x$mu_h_q025, x$mu_h_q975)
  rows <- paste0("\t\t\\texttt{", hd$sequence, "} & ", cell(hd), " & ", cell(tl),
                 " & \\texttt{", tl$sequence, "} \\\\")

  out <- c(
    "\\begin{table}[htbp]",
    "\t\\footnotesize",
    "\t\\centering",
    "\t\\renewcommand{\\arraystretch}{1.4}",
    "\t\\setlength{\\tabcolsep}{6pt}",
    "\t\\begin{tabular}{lcrr}",
    "\t\t\\hline",
    "\t\t\\textit{Head}-biased & $\\hat{\\mu}_s^h$ & $\\hat{\\mu}_s^h$ & \\textit{Tail}-biased \\\\",
    "\t\t\\hline",
    rows,
    "\t\t\\hline",
    "\t\\end{tabular}",
    paste0("\t\\caption{\\textbf{Sequences with the highest and the lowest Heads probability.} ",
           "The table shows the five sequences with the highest (left) and the five with the lowest (right) posterior mean Heads probability $\\hat{\\mu}^h_s$ in the FN treatment, with 95\\% credible intervals. ",
           "The sequences of a row are paired by rank. ",
           "All ten are classified as strongly biased.}"),
    "\t\\label{tab:rq4-top-sequences}",
    "\\end{table}")

  f <- file.path(path_tab, "rq4-top-sequences.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
