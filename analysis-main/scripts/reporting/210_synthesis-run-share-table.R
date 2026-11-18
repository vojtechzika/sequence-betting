# ============================================================
# 210_synthesis-run-share-table.R
#
# PURPOSE
#   Table of the online appendix on the share of bets against the terminal run
#   by imbalance class and degree of imbalance. Participants are classified from
#   their side choices on imbalanced sequences without a terminal run (see
#   95c_synthesis-run-share-by-imbalance.R), so that classification and outcome
#   rest on different trials. Cells show the share of bets against the run on
#   run-ending sequences, with the number of bets in parentheses.
#
# SAMPLE
#   All bets of each treatment (FN = m25, FP = m19).
#
# INPUT (data/output/csv/)
#   synthesis-run-share-by-imbalance.csv   written by 95c_synthesis-run-share-by-imbalance.R
#
# OUTPUT (data/output/tables/)
#   run-share-by-imbalance.tex   table with one panel per treatment
#
# NOTES
#   - Imbalance levels |H - T| = 0, 2, 4, 6 correspond to the 3--3, 4--2, 5--1 and 6--0 splits
#   - Output files are overwritten at every run
# ============================================================

synthesis_run_share_table <- function(cfg) {

  d <- fread(file.path(path_csv, "synthesis-run-share-by-imbalance.csv"), encoding = "UTF-8")
  cls <- c("With Imbalance" = "follower", "Against Imbalance" = "restorer", "Uncertain" = "uncertain")
  lev <- c(0, 2, 4, 6)

  # One row per class: class size and, for each imbalance level, share (number of bets)
  rows <- function(tr) {
    vapply(names(cls), function(l) {
      x <- d[treatment == tr & cls == cls[[l]]]
      cells <- vapply(lev, function(a) {
        r <- x[abs_imb == a]
        if (nrow(r) == 0L) "--" else paste0(sprintf("%.1f", 100 * r$share_against_run), " (", r$n_trials, ")")
      }, "")
      paste0("\t\t\t", l, " & ", x$class_size[1], " & ", paste(cells, collapse = " & "), " \\\\")
    }, "")
  }

  out <- c(
    "% Source: synthesis-run-share-by-imbalance.csv",
    "% Produced by: 210_synthesis-run-share-table.R",
    "\\begin{table}[h]", "\t\\centering", "\t\\footnotesize", "\t\\begin{tabular}{lrrrrr}", "\t\t\\toprule",
    "\t\t\tImbalance class & $n$ & 3--3 & 4--2 & 5--1 & 6--0 \\\\",
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{6}{l}{\\textit{FN Treatment}} \\\\", rows("m25"),
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{6}{l}{\\textit{FP Treatment}} \\\\", rows("m19"),
    "\t\t\t\\bottomrule", "\t\t\\end{tabular}",
    "\t\\caption{\\textbf{Share of bets against the terminal run on run-ending sequences, by imbalance class and degree of imbalance.}",
    "\t\tParticipants are classified from their side choices on imbalanced sequences without a terminal run (binomial test of majority-symbol choice against $0.5$ at the 5\\% level), so that classification and outcome rest on different trials. \\textit{With Imbalance} participants choose the majority symbol more often than chance, \\textit{Against Imbalance} participants the minority symbol; participants without a significant tendency are \\textit{Uncertain}. Cells show the share in percent with the number of bets in parentheses, and $n$ is the class size. These classes are based on raw choices and differ from the model-based classes of Panel~B of Table~\\ref{M-tab:synthesis-heterogeneity}.}",
    "\t\\label{Tab:run-share-by-imbalance}", "\\end{table}")

  f <- file.path(path_tab, "run-share-by-imbalance.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
