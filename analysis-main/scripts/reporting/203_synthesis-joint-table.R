# ============================================================
# 203_synthesis-joint-table.R
#
# PURPOSE
#   Main-text table of the joint model of the Synthesis: population-average
#   effects of the terminal run and of the imbalance in the betting and in the
#   side-choice equation (log-odds slope, average marginal effect, Cohen's h)
#   for FN and FP.
#
# SAMPLE
#   Full sample of each treatment (FN = m25, FP = m19), as fitted in
#   91_synthesis-stan-joint.R.
#
# INPUT (data/output/csv/)
#   synthesis-<tr>-joint-fixed-effects.csv   log-odds slopes and intercepts (94_synthesis-joint-effects.R)
#   synthesis-<tr>-joint-ame.csv             average marginal effects and baseline probabilities
#   synthesis-<tr>-joint-cohens-h.csv        Cohen's h and the discrete change in probability
#   synthesis-<tr>-joint-class-counts.csv    number of participants (sum of the class counts)
#
# OUTPUT (data/output/tables/)
#   synthesis-fixed-fn.tex   table with the betting and side-choice panels
#
# NOTES
#   - Every cell is the posterior median with the 95% credible interval below it
#   - Average marginal effects are in percentage points, baseline probabilities in percent
#   - Cohen's h of the imbalance in the side equation is the one for the fully imbalanced
#     sequences against the balanced ones (quantity side_imb_full)
#   - Output files are overwritten at every run
# ============================================================

synthesis_joint_table <- function(cfg) {

  trs <- c(FN = "m25", FP = "m19")
  rd  <- function(tr, name) fread(file.path(path_csv, paste0("synthesis-", tr, "-joint-", name, ".csv")), encoding = "UTF-8")
  dat <- lapply(trs, function(tr) list(fe = rd(tr, "fixed-effects"), ame = rd(tr, "ame"), h = rd(tr, "cohens-h"),
                                       cls = rd(tr, "class-counts")))

  # Formatting: negative numbers in math mode, intervals below the median
  num <- function(x, d = 2) {
    s <- sprintf(paste0("%.", d, "f"), x)
    if (as.numeric(s) == 0) return(sub("^-", "", s))
    if (as.numeric(s) < 0) paste0("$-", sub("^-", "", s), "$") else s
  }
  pp  <- function(x) { r <- round(100 * x, 1); paste0("$", ifelse(r < 0, "-", "+"), sprintf("%.1f", abs(r)), "$") }
  pct <- function(x) paste0(sprintf("%.1f", 100 * x), "\\%")
  cell <- function(med, lo, hi) paste0("\\makecell{", med, "\\\\\\scriptsize[", lo, ", ", hi, "]}")

  # Cells of one estimate in one treatment (t = "FN" or "FP")
  fe_cell  <- function(t, p)   { r <- dat[[t]]$fe[param == p];                        stopifnot(nrow(r) == 1L); cell(num(r$median), num(r$q025), num(r$q975)) }
  ame_cell <- function(t, q)   { r <- dat[[t]]$ame[startsWith(quantity, paste0(q, " "))]; stopifnot(nrow(r) == 1L); cell(pp(r$median), num(100 * r$q025, 1), num(100 * r$q975, 1)) }
  base_cell <- function(t, q)  { r <- dat[[t]]$ame[startsWith(quantity, paste0(q, " "))]; stopifnot(nrow(r) == 1L); cell(pct(r$median), num(100 * r$q025, 1), num(100 * r$q975, 1)) }
  h_cell   <- function(t, q)   { r <- dat[[t]]$h[quantity == q];                      stopifnot(nrow(r) == 1L); cell(num(r$median), num(r$q025), num(r$q975)) }
  line_int <- function(label, fe, ame) {
    dash <- "---"
    paste0("\t\t\t", label, " & ", fe_cell("FN", fe), " & ", base_cell("FN", ame), " & ", dash, " & ",
           fe_cell("FP", fe), " & ", base_cell("FP", ame), " & ", dash, " \\\\[6pt]")
  }
  line_est <- function(label, fe, ame, h) {
    paste0("\t\t\t", label, " & ",
           fe_cell("FN", fe), " & ", ame_cell("FN", ame), " & ", h_cell("FN", h), " & ",
           fe_cell("FP", fe), " & ", ame_cell("FP", ame), " & ", h_cell("FP", h), " \\\\[6pt]")
  }
  n_part <- vapply(dat, function(x) sum(x$cls[regressor == "bet_run_slope", n]), numeric(1))
  dpp <- vapply(dat, function(x) x$h[quantity == "bet_imb", diff_pp_median], numeric(1))

  out <- c(
    "% Source: synthesis-<m25|m19>-joint-fixed-effects.csv, -ame.csv, -cohens-h.csv, -class-counts.csv",
    "% Produced by: 203_synthesis-joint-table.R",
    "\\begin{table}[h]", "\t\\centering", "\t\\footnotesize", "\t\\setlength{\\tabcolsep}{4pt}",
    "\t{\\renewcommand{\\arraystretch}{1.3}",
    "\t\t\\begin{tabular}{lccc|ccc}", "\t\t\t\\toprule",
    "\t\t\t& \\multicolumn{3}{c|}{FN Treatment} & \\multicolumn{3}{c}{FP Treatment} \\\\",
    "\t\t\t\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    "\t\t\t& Slope & AME (pp) & Cohen's $h$ & Slope & AME (pp) & Cohen's $h$ \\\\",
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{7}{l}{\\textit{Betting equation (Bet = 1)}} \\\\[2pt]",
    line_int("Intercept", "bet_intercept", "bet_baseline"),
    line_est("\\makecell[l]{Terminal run\\\\\\scriptsize\\{$=1$\\}}", "bet_run_slope", "bet_run", "bet_run"),
    line_est("\\makecell[l]{Imbalance magnitude\\\\\\scriptsize\\{$0,\\frac13,\\frac23,1$\\}}", "bet_imb_slope", "bet_imb", "bet_imb"),
    "\t\t\t\\midrule",
    "\t\t\t\\multicolumn{7}{l}{\\textit{Side-choice equation (Heads = 1)}} \\\\[2pt]",
    line_int("Intercept", "side_intercept", "side_baseline"),
    line_est("\\makecell[l]{Terminal run\\\\\\scriptsize\\{$-1=TT$, 0, $1=HH$\\}}", "side_run_slope", "side_run", "side_run"),
    line_est("\\makecell[l]{Imbalance\\\\\\scriptsize\\{$-1$ (T) to $1$ (H), by $\\frac13$\\}}", "side_imb_slope", "side_imb", "side_imb_full"),
    "\t\t\t\\bottomrule", "\t\t\\end{tabular}}",
    "\t\\caption{",
    "\t\t\\textbf{Population-average effects of the terminal run and imbalance, by treatment.}",
    paste0("\t\tEach cell reports the posterior median with its 95\\% credible interval in brackets below. Curly braces $\\{\\}$ beneath each regressor give its coding or admissible value set. FN: ",
           n_part[["FN"]], " participants; FP: ", n_part[["FP"]], " participants."),
    "\t\t\\textit{Slope} is the log-odds coefficient of the fixed effect. \\textit{AME} is the population-averaged marginal effect on the probability scale, integrating over the fitted participant- and sequence-level random effects rather than assuming them at zero. For the Intercept rows it is the RE-integrated baseline predicted probability (no terminal run, zero imbalance). For the other rows it is the instantaneous per-unit marginal effect ($dp/dx$) evaluated at each observation's actual data: because the betting-equation regressors are coded 0 to 1, ``per unit'' spans their full observed range; because the side-choice regressors are coded $-1$ to $1$, ``per unit'' corresponds to moving from neutral (no cue) to the cue fully present in one direction, half of the total coded range.",
    paste0("\t\t\\textit{Cohen's $h$}, $2\\arcsin\\sqrt{p_1}-2\\arcsin\\sqrt{p_0}$, is computed per posterior draw from the discrete change in probability when the cue moves from its reference to its treated level on the actual trials (run present vs.\\ absent, and fully imbalanced vs.\\ balanced, for betting; the probability of choosing the run's symbol and the majority symbol, respectively, for side choice). For imbalance and betting this discrete change is ",
           sprintf("%.1f", dpp[["FN"]]), "pp (FN) and ", sprintf("%.1f", dpp[["FP"]]),
           "pp (FP), below the AME, which is a per-unit slope. AME and $h$ are evaluated at the estimated random effects of the participants and sequences in the sample, so their intervals need not coincide with those of the slope."),
    "\t}", "\t\\label{tab:synthesis-fixed}", "\\end{table}")

  f <- file.path(path_tab, "synthesis-fixed-fn.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
