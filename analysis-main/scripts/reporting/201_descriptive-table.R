# ============================================================
# 201_descriptive-table.R
#
# PURPOSE
#   Table of descriptive statistics by treatment in the main text: sample
#   composition, Holt-Laury risk preferences, the sequence betting task and
#   covariates, for the FN treatment (full and confirmatory sample) and the
#   FP treatment (full sample).
#
# SAMPLES
#   FN full          all participants of the FN treatment (m25)
#   FN confirmatory  FN participants with a consistent Holt-Laury choice
#                    pattern (m25/consistent), all of whom are normative betters
#   FP full          all participants of the FP treatment (m19)
#
# INPUT
#   data/output/csv/descriptive-participants-summary.csv   rows m25, m25/consistent, m19
#   data/output/csv/descriptive-sequences-summary.csv      rows m25, m25/consistent, m19
#   data/output/csv/processing_lotr-alpha.csv              Cronbach's alpha of the LOT-R, by treatment
#   model "mpl": a-star/a-star-flags-<tr>-full.csv         posterior probabilities that a zero (P0) or a
#                                                          positive (Pplus) stake is EU-optimal
#   model "mpl": scored-<tr>.csv                           Holt-Laury consistency flag
#
# OUTPUT
#   data/output/tables/descriptive.tex
#
# NOTES
#   - Entries are means with standard deviations in square brackets, counts or shares
#   - The standard deviation of the share of bets is that of a binary variable,
#     sqrt(p * (1 - p)); the standard deviations of the staking and decision-time
#     rows are taken over the sequence-level observations, the others over participants
#   - Stakes are expressed as a share of the 100-ECU endowment of each round
#   - Normative non-betters and betters are the participants with P0 or Pplus at least
#     as large as the primary threshold cfg$design$a_flags$tau[1]; participants who
#     meet neither criterion are reported as unclassified
#   - The classification of the confirmatory sample is restricted to participants
#     with a consistent Holt-Laury choice pattern
#   - The rows m25/consistent of the two summary files are written by
#     01_descriptive-participants.R and 02_descriptive-sequences.R
#   - The output file is overwritten at every run
# ============================================================

descriptive_table <- function(cfg) {
  
  tau       <- as.numeric(cfg$design$a_flags$tau[1])   # threshold of the normative classification
  endowment <- 100                                     # ECU per round
  
  # ---- Inputs ----
  f_par <- file.path(path_csv, "descriptive-participants-summary.csv")
  f_seq <- file.path(path_csv, "descriptive-sequences-summary.csv")
  f_alp <- file.path(path_csv, "processing_lotr-alpha.csv")
  stopifnot(file.exists(f_par), file.exists(f_seq), file.exists(f_alp))
  
  par <- fread(f_par, encoding = "UTF-8")
  seq <- fread(f_seq, encoding = "UTF-8")
  alp <- fread(f_alp, encoding = "UTF-8")
  
  # Table columns: FN full, FN confirmatory, FP full
  grp <- c("m25", "m25/consistent", "m19")
  if (!all(grp %in% par$group) || !all(grp %in% seq$treatment))
    stop("Rows ", paste(grp, collapse = ", "), " are missing in the summary files; ",
         "run 01_descriptive-participants.R and 02_descriptive-sequences.R first.")
  p <- par[match(grp, group)]
  s <- seq[match(grp, treatment)]
  
  # ---- Normative classification ----
  classify <- function(tr, consistent_only) {
    fl <- fread(model_path("mpl", file.path("a-star", paste0("a-star-flags-", tr, "-full.csv"))),
                encoding = "UTF-8", colClasses = c(pid = "character"))
    if (consistent_only) {
      sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
      fl <- fl[pid %in% as.character(sc[as.integer(inconsistent) == 0L, pid])]
    }
    c(non = sum(fl$P0 >= tau), bet = sum(fl$Pplus >= tau), n = nrow(fl))
  }
  cl  <- list(classify("m25", FALSE), classify("m25", TRUE), classify("m19", FALSE))
  cls <- function(k) vapply(cl, function(x) x[[k]], numeric(1))
  
  # ---- Cell formatting ----
  pct <- function(x) sprintf("%.0f\\%%", 100 * x)
  ms  <- function(m, sd, d = 1) sprintf(paste0("%.", d, "f [%.", d, "f]"), m, sd)
  row <- function(label, x, tail = "") paste0("\t\t", label, " & ", paste(x, collapse = " & "), " \\\\", tail)
  
  alpha_fn <- sprintf("%.2f", alp[sample == "m25", raw_alpha])
  alpha_fp <- sprintf("%.2f", alp[sample == "m19", raw_alpha])
  
  out <- c(
    paste0("% Sources: descriptive-participants-summary.csv, descriptive-sequences-summary.csv (rows m25, m25/consistent, m19), ",
           "processing_lotr-alpha.csv, a-star-flags-<m25|m19>-full.csv, scored-<m25|m19>.csv"),
    "% Produced by: 201_descriptive-table.R",
    "\\begin{table}[!htbp]",
    "\t\\footnotesize",
    "\t\\centering",
    "\t\\setlength{\\tabcolsep}{4pt}",
    "\t\\begin{tabular}{lccc}",
    "\t\t\\toprule",
    "\t\t& \\multicolumn{2}{c}{FN ($m=2.5$)} & FP ($m=1.9$) \\\\",
    "\t\t\\cmidrule(lr){2-3}\\cmidrule(lr){4-4}",
    "\t\t& Full & Confirmatory & Full \\\\",
    "\t\t\\midrule",
    row("Sample Size", p$N),
    row("\\quad of which Women, Men, Unspecified", sprintf("%d, %d, %d", p$sex_F_n, p$sex_M_n, p$sex_O_n)),
    row("Age in Years", ms(p$age_mean, p$age_sd)),
    row("Tasks' Earnings in ECU", sprintf("%.0f [%.1f]", p$payoff_mean, p$payoff_sd)),
    "\t\t\\midrule",
    "\t\t\\multicolumn{4}{l}{\\textit{Holt--Laury Task}} \\\\",
    "\t\t\\midrule",
    row("$\\bar{r}$ (Bayesian)", ms(p$r_bayes_mean, p$r_bayes_sd, 2)),
    row("$\\bar{r}$ (Switching Point)", c("--", ms(p$r_classic_mean[2], p$r_classic_sd[2], 2), "--")),
    row("Share of Consistent Participants", pct(1 - p$hl_inconsistent_frac)),
    "\t\t\\midrule",
    "\t\t\\multicolumn{4}{l}{\\textit{Sequence Betting Task}} \\\\",
    "\t\t\\midrule",
    row("Sequence-Level Observations", format(s$n_obs, big.mark = "{,}", trim = TRUE)),
    row("Share of Bets Placed", sprintf("%s [%s]", pct(s$bet_rate_mean), pct(sqrt(s$bet_rate_mean * (1 - s$bet_rate_mean))))),
    row("\\quad Normative Non-Betters", cls("non")),
    row("\\quad Normative Betters", cls("bet")),
    row("\\quad Unclassified", cls("n") - cls("non") - cls("bet")),
    row("Share of Endowment Staked in Bets Placed",
        sprintf("%s [%s]", pct(s$stake_mean_bet / endowment), pct(s$stake_sd_bet / endowment)), " % converted from ECU"),
    row("Share of Heads in Bets Placed", pct(s$side_H_rate_cond_bet)),
    "\t\t\\midrule",
    "\t\t\\multicolumn{4}{l}{\\textit{Covariate Measures}} \\\\",
    "\t\t\\midrule",
    row("Seconds to Decide per Sequence", ms(s$screen_ms_mean / 1000, s$screen_ms_sd / 1000)),
    row("Life Orientation Test Revised (0--24 pts.)", ms(p$lotr_mean, p$lotr_sd)),
    row("Life Orientation Test Revised (Cronbach's $\\alpha)$", c(alpha_fn, "--", alpha_fp)),
    "\t\t\\bottomrule",
    "\t\\end{tabular}",
    "\t\\caption{",
    "\t\t\t\\textbf{Descriptive statistics by treatment.}",
    paste0("\t\t\tEntries are means with standard deviations in square brackets, or counts, or shares. ",
           "The confirmatory sample, as per Section~\\ref{Sec:ModellingStrategyEcon} (\\nameref{Sec:SampleEcon}), ",
           "consists of the FN participants with a consistent Holt--Laury choice pattern, all of whom are normative betters; ",
           "the FP treatment has no confirmatory sample. ",
           "Standard deviations of the betting, staking and decision-time rows are across sequences, and the others across participants. ",
           "The switching-point estimate exists only for participants with a consistent Holt--Laury choice pattern. ",
           "Cronbach's $\\alpha$ is computed by treatment."),
    paste0("\t\t\tNormative betters and non-betters are model-derived classifications of participants by the posterior probability ",
           "that a positive or a zero stake is EU-optimal, defined in \\OAsec{App:econ-sample}. ",
           "Participants for whom neither probability reaches ", sprintf("%.2f", tau), " are unclassified."),
    "\t}",
    "\t\\label{tab:descriptive}",
    "\\end{table}")
  
  f <- file.path(path_tab, "descriptive.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
