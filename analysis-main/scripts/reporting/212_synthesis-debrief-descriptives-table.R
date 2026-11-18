# ============================================================
# 212_synthesis-debrief-descriptives-table.R
#
# PURPOSE
#   Table of the online appendix with descriptive statistics of the closed
#   debriefing questions (Q1 to Q7) and the answers to Q8, by treatment and pooled.
#
# INPUT (data/output/csv/)
#   synthesis-debrief-descriptives.csv            written by 96b_synthesis-debrief-descriptives.R
#   synthesis-debrief-used-strategy-shares.csv    written by 96b_synthesis-debrief-descriptives.R
#
# OUTPUT (data/output/tables/)
#   debrief-descriptives.tex
#
# NOTES
#   - Cells show mean (standard deviation) for each treatment and pooled, and the
#     median and range for the pooled sample
#   - Output files are overwritten at every run
# ============================================================

synthesis_debrief_descriptives_table <- function(cfg) {

  d <- fread(file.path(path_csv, "synthesis-debrief-descriptives.csv"), encoding = "UTF-8")
  s <- fread(file.path(path_csv, "synthesis-debrief-used-strategy-shares.csv"), encoding = "UTF-8")

  items <- c(perceived_win_probability = "Q1 Expected wins in 100 rounds",
             perceived_risk            = "Q2 Riskiness of the betting",
             belief_independence       = "Q3 Belief in independence of tosses",
             reliance_on_sequence      = "Q4 Use of the six-toss history",
             action_seeking            = "Q5 Betting because acting feels better",
             fatigue                   = "Q6 Boredom or tiredness",
             self_risk_tolerance       = "Q7 Willingness to take financial risks")

  ms <- function(it, sm) { r <- d[item == it & sample == sm]; sprintf("%.1f (%.1f)", r$mean, r$sd) }

  rows <- vapply(names(items), function(it) {
    a <- d[item == it & sample == "all"]
    paste0("\t\t\t", items[[it]], " & ", ms(it, "m25"), " & ", ms(it, "m19"), " & ", ms(it, "all"),
           " & ", format(a$median, trim = TRUE), " & ", a$min, "--", a$max, " \\\\")
  }, "")

  sh <- function(ans, sm) { r <- s[answer == ans & sample == sm]; sprintf("%d (%.0f\\%%)", r$n, 100 * r$share) }
  q8 <- vapply(c(yes = "Yes", no = "No", unsure = "Not sure"), function(a) "", "")
  q8 <- vapply(names(q8), function(a) {
    paste0("\t\t\t\\hspace{1em}", c(yes = "Yes", no = "No", unsure = "Not sure")[[a]], " & ",
           sh(a, "m25"), " & ", sh(a, "m19"), " & ", sh(a, "all"), " & & \\\\")
  }, "")

  n <- function(sm) d[item == "perceived_risk" & sample == sm, n]

  out <- c(
    "% Source: synthesis-debrief-descriptives.csv, synthesis-debrief-used-strategy-shares.csv",
    "% Produced by: 212_synthesis-debrief-descriptives-table.R",
    "\\begin{table}[h]", "\t\\centering", "\t\\footnotesize", "\t\\begin{tabular}{lccccc}", "\t\t\\toprule",
    "\t\t\tItem & FN & FP & All & Median & Range \\\\",
    paste0("\t\t\t$n$ & ", n("m25"), " & ", n("m19"), " & ", n("all"), " & & \\\\"),
    "\t\t\t\\midrule", rows,
    "\t\t\t\\midrule",
    "\t\t\tQ8 Used a rule or strategy & & & & & \\\\", q8,
    "\t\t\t\\bottomrule", "\t\t\\end{tabular}",
    "\t\\caption{\\textbf{Descriptive statistics of the debriefing questions.} Cells for Q1 to Q7 show the mean with the standard deviation in parentheses for each treatment and for all participants, followed by the median and the range for all participants; Q8 shows counts with shares. Item numbers and wording are those of the post-task questionnaire (Section~\\ref{App:PostTaskQuestionsEn}). Q1 is a number between 0 and 100, Q2 to Q7 are answered on a 0 to 10 scale.}",
    "\t\\label{Tab:debrief-descriptives}", "\\end{table}")

  f <- file.path(path_tab, "debrief-descriptives.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
