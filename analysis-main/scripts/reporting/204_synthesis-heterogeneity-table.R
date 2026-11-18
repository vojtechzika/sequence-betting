# ============================================================
# 204_synthesis-heterogeneity-table.R
#
# PURPOSE
#   Main-text table on the participant-level heterogeneity of the
#   joint model of the Synthesis. Panel A: posterior standard deviation (tau) of
#   the participant-level random slopes and its ratio to the absolute population
#   slope. Panel B: classification of the participants by the 95% credible
#   interval of their personal side-choice slope (against, with, uncertain).
#   Panel C: correlations between the participant effects of the betting and of the
#   side-choice equation.
#
# SAMPLE
#   Full sample of each treatment (FN = m25, FP = m19), as fitted in
#   91_synthesis-stan-joint.R.
#
# INPUT (data/output/csv/)
#   synthesis-<tr>-joint-tau.csv           tau, its 95% credible interval, and tau / |fixed effect|
#                                          (94_synthesis-joint-effects.R)
#   synthesis-<tr>-joint-dissociation.csv  cross-equation correlations of the participant effects
#   synthesis-<tr>-joint-class-counts.csv  number and share of participants by class
#                                          (95_synthesis-joint-heterogeneity.R)
#
# OUTPUT (data/output/tables/)
#   synthesis-heterogeneity.tex   table with panels A and B
#
# NOTES
#   - The ratio in bold is the largest of each treatment
#   - The number of participants in the column labels is the sum of the class counts
#   - Percentages are rounded to integers and may not add to 100
#   - Output files are overwritten at every run
# ============================================================

synthesis_heterogeneity_table <- function(cfg) {

  trs <- c(FN = "m25", FP = "m19")
  rd  <- function(tr, name) fread(file.path(path_csv, paste0("synthesis-", tr, "-joint-", name, ".csv")), encoding = "UTF-8")
  tau <- lapply(trs, rd, name = "tau")
  cls <- lapply(trs, rd, name = "class-counts")
  rho <- lapply(trs, rd, name = "dissociation")
  n_part <- vapply(cls, function(x) sum(x[regressor == "side_run_slope", n]), numeric(1))

  # Betting-equation classes quoted in the caption
  share_unc <- vapply(cls, function(x) x[regressor == "bet_run_slope" & class == "uncertain", share], numeric(1))
  neg_imb   <- vapply(cls, function(x) { r <- x[regressor == "bet_imb_slope" & class == "disengaged", n]; if (length(r)) r else 0L }, numeric(1))

  # Panel A: one row per random slope
  rows <- c("(Bet) Has Terminal Run"    = "bet_run_sd",
            "(Bet) Imbalance Magnitude" = "bet_imb_sd",
            "(Heads) Terminal Run"      = "side_run_sd",
            "(Heads) Imbalance"         = "side_imb_sd")
  cell_a <- function(t, p) {
    r <- tau[[t]][param == p]; stopifnot(nrow(r) == 1L)
    ratio <- sprintf("%.2f", r$heterogeneity_ratio)
    if (r$heterogeneity_ratio == max(tau[[t]][param %in% rows, heterogeneity_ratio])) ratio <- paste0("\\textbf{", ratio, "}")
    paste0(sprintf("%.2f [%.2f, %.2f]", r$median, r$q025, r$q975), " & ", ratio)
  }
  body_a <- vapply(seq_along(rows), function(i) {
    paste0("\t\t\t", sprintf("%-25s", names(rows)[i]), " & ", cell_a("FN", rows[[i]]), " & ", cell_a("FP", rows[[i]]), " \\\\")
  }, "")

  # Panel C: correlations of the participant effects, median and interval on one line
  fr <- function(x) {
    s <- sub("^0", "", sprintf("%.2f", abs(x)))
    if (x < 0 && s != ".00") paste0("$-", s, "$") else s
  }
  cell_c <- function(t, k) {
    r <- rho[[t]][correlation == k]; stopifnot(nrow(r) == 1L)
    paste0(fr(r$median), " [", fr(r$q025), ", ", fr(r$q975), "]")
  }
  line_c <- function(label, k) {
    paste0("\t\t\t", label, " & \\multicolumn{2}{c|}{", cell_c("FN", k), "} & \\multicolumn{2}{c}{", cell_c("FP", k), "} \\\\")
  }

  # Panel B: counts and shares of the side-choice classes
  cell_b <- function(t, reg, k) {
    x <- cls[[t]][regressor == reg]
    one <- function(cl) { r <- x[class == cl]; if (nrow(r) == 0L) "0 (0\\%)" else sprintf("%d (%.0f\\%%)", r$n, 100 * r$share) }
    paste(vapply(k, one, ""), collapse = " / ")
  }
  line_b <- function(label, reg, k) {
    paste0("\t\t\t", label, " & \\multicolumn{2}{c|}{", cell_b("FN", reg, k), "} & \\multicolumn{2}{c}{", cell_b("FP", reg, k), "} \\\\")
  }

  out <- c(
    "% Source: synthesis-<m25|m19>-joint-tau.csv, -class-counts.csv, -dissociation.csv",
    "% Produced by: 204_synthesis-heterogeneity-table.R",
    "\\begin{table}[h]", "\t\\centering", "\t\\footnotesize", "\t\\begin{tabular}{lcc|cc}", "\t\t\\toprule",
    paste0("\t\t\t& \\multicolumn{2}{c|}{FN Treatment ($n=", n_part[["FN"]], "$)} & \\multicolumn{2}{c}{FP Treatment ($n=", n_part[["FP"]], "$)} \\\\"),
    "\t\t\tRegressor & $\\tau$ [95\\% CrI] & $\\tau/|\\hat\\beta|$ & $\\tau$ [95\\% CrI] & $\\tau/|\\hat\\beta|$ \\\\",
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{5}{l}{\\textit{Panel A: Participant-level heterogeneity}} \\\\",
    body_a,
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{5}{l}{\\textit{Panel B: Participants with a credible personal side-choice slope}} \\\\",
    line_b("Run  (against / with / uncertain)",       "side_run_slope", c("against_run", "with_run", "uncertain")),
    line_b("Imbalance  (against / with / uncertain)", "side_imb_slope", c("against_imb", "with_imb", "uncertain")),
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{5}{l}{\\textit{Panel C: Correlation between the betting and side-choice participant effects ($\\rho$)}} \\\\",
    line_c("Intercept",    "bet_intercept <-> side_intercept"),
    line_c("Terminal run", "bet_run_slope <-> side_run_slope"),
    line_c("Imbalance",    "bet_imb_slope <-> side_imb_slope"),
    "\t\t\t\\bottomrule", "\t\\end{tabular}",
    paste0("\t\\caption{\\textbf{Participant-level heterogeneity in run- and imbalance-sensitivity, by treatment.} Panel A: posterior SD ($\\tau$) of each participant-level random slope and its ratio to the absolute population slope (Table~\\ref{tab:synthesis-fixed}). Panel B: each participant's personal side-choice slope (population slope plus their own random-effect deviation), classified as credibly positive (with the cue), credibly negative (against the cue), or indistinguishable from zero (uncertain) at the 95\\% level, separately by treatment. Panel C: posterior median and 95\\% credible interval of the correlation between the participant-level random intercepts, run slopes and imbalance slopes of the betting and the side-choice equation. The classification of the betting-equation slopes is omitted: ",
           sprintf("%.0f\\%% of the participants in FN and %.0f\\%% in FP are uncertain on the run slope, and only %.0f of %d in FN and %.0f of %d in FP are credibly negative on the imbalance slope, so direction is not the informative axis for the decision to bet as it is for side choice.}",
                   100 * share_unc[["FN"]], 100 * share_unc[["FP"]], neg_imb[["FN"]], as.integer(n_part[["FN"]]), neg_imb[["FP"]], as.integer(n_part[["FP"]]))),
    "\t\\label{tab:synthesis-heterogeneity}", "\\end{table}")

  f <- file.path(path_tab, "synthesis-heterogeneity.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
