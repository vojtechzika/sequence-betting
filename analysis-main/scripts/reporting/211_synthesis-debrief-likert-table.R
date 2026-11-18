# ============================================================
# 211_synthesis-debrief-likert-table.R
#
# PURPOSE
#   Table of the online appendix with the closed debriefing items by the
#   participants' slope classes. Cells show the mean and, in parentheses, the
#   median of each item, for the imbalance classes (With Imbalance, Against
#   Imbalance, Uncertain) and the run classes (Against Run, Uncertain).
#
# SAMPLE
#   All 202 participants who completed the debriefing questionnaire (FN and FP
#   pooled). Classes are the model-based classes of the joint model.
#
# INPUT (data/output/csv/)
#   synthesis-debrief-likert-by-class.csv   written by 96_synthesis-debrief-validation.R
#
# OUTPUT (data/output/tables/)
#   debrief-likert-by-class.tex
#
# NOTES
#   - Items are answered on a 0 to 10 scale, except Q1 (0 to 100)
#   - The run class With Run (n = 2) is omitted
#   - Output files are overwritten at every run
# ============================================================

synthesis_debrief_likert_table <- function(cfg) {

  d <- fread(file.path(path_csv, "synthesis-debrief-likert-by-class.csv"), encoding = "UTF-8")

  # Items in the order of the post-task questionnaire (Q6, boredom or tiredness, is not analyzed by class)
  items <- c(perceived_win_probability  = "Q1 Expected wins in 100 rounds",
             perceived_risk             = "Q2 Riskiness of the betting",
             belief_independence        = "Q3 Belief in independence of tosses",
             reliance_on_sequence       = "Q4 Use of the six-toss history",
             action_seeking             = "Q5 Betting because acting feels better",
             self_risk_tolerance        = "Q7 Willingness to take financial risks")

  specs <- list(c("imb_class", "with_imb"), c("imb_class", "against_imb"), c("imb_class", "uncertain"),
                c("run_class", "against_run"), c("run_class", "uncertain"))

  cell <- function(it, dim, cl) {
    r <- d[item == it & class_dim == dim & class == cl]
    paste0(sprintf("%.1f", r$mean), " (", format(r$median, nsmall = 0, trim = TRUE), ")")
  }
  n_of <- function(dim, cl) d[class_dim == dim & class == cl, n][1]

  body <- vapply(names(items), function(it) {
    paste0("\t\t\t", items[[it]], " & ",
           paste(vapply(specs, function(s) cell(it, s[1], s[2]), ""), collapse = " & "), " \\\\")
  }, "")

  out <- c(
    "% Source: synthesis-debrief-likert-by-class.csv",
    "% Produced by: 211_synthesis-debrief-likert-table.R",
    "\\begin{table}[h]", "\t\\centering", "\t\\footnotesize", "\t\\begin{tabular}{lccccc}", "\t\t\\toprule",
    "\t\t\t& \\multicolumn{3}{c}{Imbalance class} & \\multicolumn{2}{c}{Run class} \\\\",
    "\t\t\t\\cmidrule(lr){2-4} \\cmidrule(lr){5-6}",
    "\t\t\tItem & With & Against & Uncertain & Against & Uncertain \\\\",
    paste0("\t\t\t$n$ & ", paste(vapply(specs, function(s) as.character(n_of(s[1], s[2])), ""), collapse = " & "), " \\\\"),
    "\t\t\t\\midrule", body,
    "\t\t\t\\bottomrule", "\t\t\\end{tabular}",
    "\t\\caption{\\textbf{Closed debriefing items by slope class.} Cells show the mean and, in parentheses, the median. Item numbers and wording are those of the post-task questionnaire (Section~\\ref{App:PostTaskQuestionsEn}); items are answered on a 0 to 10 scale, except Q1, which is a number between 0 and 100. Classes are the model-based classes of Panel~B of Table~\\ref{M-tab:synthesis-heterogeneity} (FN and FP pooled, $N=202$); the two participants with a credible run slope in the direction of the run are omitted.}",
    "\t\\label{Tab:debrief-likert-by-class}", "\\end{table}")

  f <- file.path(path_tab, "debrief-likert-by-class.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
