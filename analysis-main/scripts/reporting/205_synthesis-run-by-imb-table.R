# ============================================================
# 205_synthesis-run-by-imb-table.R
#
# PURPOSE
#   Table of the online appendix on the personal run slope by imbalance class.
#   Participants are grouped by the classification of their personal side-choice
#   imbalance slope (against, uncertain, with). For each group the table reports
#   the number of participants, the median personal run slope, and the share of
#   participants whose run slope is credibly negative (against the run). It also
#   reports Mann-Whitney tests of the run slopes between the groups and the
#   Spearman correlation between the continuous run and imbalance slopes.
#
# SAMPLE
#   Full sample of each treatment (FN = m25, FP = m19), as fitted in
#   91_synthesis-stan-joint.R.
#
# INPUT (data/output/csv/)
#   synthesis-<tr>-joint-run-slope-participants.csv   personal side run slope (posterior median), its 95% interval, and class
#   synthesis-<tr>-joint-imb-slope-participants.csv   personal side imbalance slope, its 95% interval, and class
#   (both written by 95_synthesis-joint-heterogeneity.R)
#
# OUTPUT (data/output/tables/)
#   run-by-imb-class.tex   table with the groups, the tests, and the correlation
#
# NOTES
#   - Negative run slopes mean that the participant bets against the run
#   - Tests are two-sided Mann-Whitney tests (wilcox.test with default settings)
#   - The Spearman correlation is computed on the posterior medians of the personal slopes
#     and does not carry the posterior uncertainty of the slopes
#   - Output files are overwritten at every run
# ============================================================

run_by_imb_table <- function(cfg) {

  trs <- c(FN = "m25", FP = "m19")
  rd  <- function(tr, name) fread(file.path(path_csv, paste0("synthesis-", tr, "-joint-", name, "-slope-participants.csv")),
                                  encoding = "UTF-8", colClasses = c(pid = "character"))
  dat <- lapply(trs, function(tr) merge(rd(tr, "run")[, .(pid, run = median, run_class = class)],
                                        rd(tr, "imb")[, .(pid, imb = median, imb_class = class)], by = "pid"))

  grp <- c("Against Imbalance" = "against_imb", "Uncertain" = "uncertain", "With Imbalance" = "with_imb")
  fp  <- function(p) if (p < 0.001) "$<.001$" else sub("^0", "", sprintf("%.3f", p))
  num <- function(x) if (round(x, 2) < 0) paste0("$-", sprintf("%.2f", abs(x)), "$") else sprintf("%.2f", x)

  # One row per imbalance class: n, median run slope, share credibly against the run
  cells <- function(t, g) {
    x <- dat[[t]][imb_class == g]
    paste0(nrow(x), " & ", num(median(x$run)), " & ", sprintf("%.1f", 100 * mean(x$run_class == "against_run")), "\\%")
  }
  body <- vapply(names(grp), function(l) paste0("\t\t\t", l, " & ", cells("FN", grp[[l]]), " & ", cells("FP", grp[[l]]), " \\\\"), "")

  # Mann-Whitney tests of the run slopes between the classes
  mw <- function(t, a, b) wilcox.test(dat[[t]][imb_class == a, run], dat[[t]][imb_class == b, run])$p.value
  tst <- function(label, a, b) {
    paste0("\t\t\t", label, " & \\multicolumn{3}{c|}{", fp(mw("FN", a, b)), "} & \\multicolumn{3}{c}{", fp(mw("FP", a, b)), "} \\\\")
  }
  rho <- function(t) {
    r <- cor.test(dat[[t]]$run, dat[[t]]$imb, method = "spearman", exact = FALSE)
    p <- if (r$p.value < 0.001) "p<.001" else paste0("p=", sub("^0", "", sprintf("%.3f", r$p.value)))
    paste0(sub("^0", "", sprintf("%.2f", r$estimate)), " ($", p, "$)")
  }

  out <- c(
    "% Source: synthesis-<m25|m19>-joint-run-slope-participants.csv, synthesis-<m25|m19>-joint-imb-slope-participants.csv",
    "% Produced by: 205_synthesis-run-by-imb-table.R",
    "\\begin{table}[h]", "\t\\centering", "\t\\footnotesize", "\t\\begin{tabular}{lrrr|rrr}", "\t\t\\toprule",
    "\t\t\t& \\multicolumn{3}{c|}{FN Treatment} & \\multicolumn{3}{c}{FP Treatment} \\\\",
    "\t\t\t\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    "\t\t\tImbalance class & $n$ & Median run slope & Against the run & $n$ & Median run slope & Against the run \\\\",
    "\t\t\t\\midrule", body,
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{7}{l}{\\textit{Mann--Whitney test of the run slopes, $p$-value}} \\\\",
    tst("Against vs.\\ With Imbalance", "against_imb", "with_imb"),
    tst("Against Imbalance vs.\\ Uncertain", "against_imb", "uncertain"),
    tst("With Imbalance vs.\\ Uncertain", "with_imb", "uncertain"),
    "\t\t\t\\midrule",
    paste0("\t\t\tSpearman $\\rho$, run and imbalance slopes & \\multicolumn{3}{c|}{", rho("FN"), "} & \\multicolumn{3}{c}{", rho("FP"), "} \\\\"),
    "\t\t\t\\bottomrule", "\t\t\\end{tabular}",
    "\t\\caption{\\textbf{Personal run slope by imbalance class, by treatment.}",
    "\t\tParticipants are grouped by the 95\\% credible interval of their personal side-choice imbalance slope (Panel B of Table~\\ref{M-tab:synthesis-heterogeneity}). \\textit{Median run slope} is the median across participants of the posterior median of the personal side-choice run slope, and negative values mean that the participant bets against the run. \\textit{Against the run} is the share of participants whose run slope is credibly negative at the 95\\% level. The tests are two-sided Mann--Whitney tests of the run slopes between the classes. The last row is the Spearman correlation between the posterior medians of the personal run and imbalance slopes across all participants of the treatment.}",
    "\t\\label{Tab:run-by-imb-class}", "\\end{table}")

  f <- file.path(path_tab, "run-by-imb-class.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
