# ============================================================
# 209_synthesis-top-sequences-table.R
#
# PURPOSE
#   Table of the main text of the Synthesis: the deviations of the sequence-level
#   estimates from their grand means (in percentage points, one decimal) for three groups of five sequences: the largest welfare-loss
#   deviations, the most Head-biased, and the most Tail-biased sequences. An
#   asterisk marks a 95% credible interval that excludes zero.
#
# SAMPLE
#   FN, first treatment of cfg$run$treatment: confirmatory sample for the
#   betting rate, the stake deviation, and the welfare loss, full sample for
#   the Heads probability (as set in 91_synthesis-sequences.R).
#
# INPUT (data/output/csv/)
#   synthesis-<treatment>-confirmatory-sequence-diffs.csv   (91_synthesis-sequences.R)
#   synthesis-<treatment>-confirmatory-grand-means.csv      grand medians of the four estimands
#
# OUTPUT (data/output/tables/)
#   synthesis-sequences-top.tex   table with the grand means below the column labels
#
# NOTES
#   - Groups are chosen by the posterior median of the deviation (largest welfare
#     loss, largest and smallest Heads share), not by the credible interval
#   - A deviation that rounds to zero is printed as +0.0
#   - The grand means below the column labels are in percent
#   - Output files are overwritten at every run
# ============================================================

synthesis_top_sequences_table <- function(cfg) {

  tr <- unique(as.character(cfg$run$treatment))[1L]
  d  <- fread(file.path(path_csv, paste0("synthesis-", tr, "-confirmatory-sequence-diffs.csv")), encoding = "UTF-8")
  g  <- fread(file.path(path_csv, paste0("synthesis-", tr, "-confirmatory-grand-means.csv")), encoding = "UTF-8")
  stopifnot(nrow(d) == 64L)

  cell <- function(x, sig) {
    r <- round(100 * x, 1)
    paste0("$", ifelse(r < 0, "-", "+"), sprintf("%.1f", abs(r)), ifelse(sig, "^*", ""), "$")
  }
  row <- function(i) {
    paste0("\\texttt{", d$sequence[i], "} & ",
           paste(cell(unlist(d[i, .(h_d_median, b_d_median, a_d_median, c_d_median)]),
                      as.logical(unlist(d[i, .(h_sig, b_sig, a_sig, c_sig)]))), collapse = " & "),
           " \\\\")
  }
  grp <- list(
    "Top-5 welfare-losing sequences" = order(-d$c_d_median)[1:5],
    "Top-5 head-biased sequences"    = order(-d$h_d_median)[1:5],
    "Top-5 tail-biased sequences"    = order( d$h_d_median)[1:5]
  )
  gm <- function(e) paste0(sprintf("%.1f", 100 * g[estimand == e, grand_median]), "\\%")

  body <- character(0)
  for (k in seq_along(grp)) {
    if (k > 1L) body <- c(body, "\t\t\\midrule")
    body <- c(body, paste0("\t\t\\multicolumn{5}{l}{\\textit{", names(grp)[k], "}} \\\\"),
              paste0("\t\t", vapply(grp[[k]], row, "")))
  }

  out <- c(
    "\\begin{table}[H]", "\t\\centering", "\t\\footnotesize", "\t\\begin{tabular}{lcccc}", "\t\t\\toprule",
    "\t\t& \\textit{Heads-share $\\Delta$} & \\textit{Betting-rate $\\Delta$} & \\textit{Stake-deviation $\\Delta$} & \\textit{Welfare-loss $\\Delta$} \\\\",
    paste0("\t\t& $\\bar\\mu^h = ", gm("h"), "$ & $\\bar\\mu^b = ", gm("b"), "$ & $\\bar\\mu^a = ", gm("a"), "$ & $\\bar\\mu^c = ", gm("c"), "$ \\\\"),
    "\t\t\\midrule", body, "\t\t\\bottomrule", "\t\\end{tabular}",
    "\t\\caption{",
    "\t\tEach cell reports the difference, in percentage points, between the posterior median for the sequence shown in the row (computed over all participants) and the corresponding grand mean for the quantity shown in the column. Grand means are computed across all participants and sequences and are reported in percent below the column labels.",
    "\t\t\\textit{Heads share} and \\textit{Betting rate} are probabilities, whereas \\textit{Stake deviation} and \\textit{Welfare loss} are shares of the endowment.",
    "\t\tAll quantities are estimated using data from the FN treatment. \\textit{Heads share} uses the full sample; the other three quantities use the confirmatory sample.",
    "\t\tAn asterisk denotes a 95\\% credible interval excluding zero.",
    "\t}", "\t\\label{tab:synthesis_sequences}", "\\end{table}")

  f <- file.path(path_tab, "synthesis-sequences-top.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
