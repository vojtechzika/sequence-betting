# ============================================================
# 202_model-summary-table.R
#
# PURPOSE
#   Table of model summaries in the main text: betting rate (RQ1), stake
#   deviation with its all-in component (RQ2) and welfare loss (RQ3) for the
#   confirmatory sample of the FN treatment. Rows: EU benchmark, grand mean
#   with its 95% credible interval, median participant, standard deviations
#   of the participant and sequence effects and their ratio, classification
#   of participants and of sequences, and the sample sizes.
#
# INPUT (data/output/csv/, FN treatment m25, confirmatory sample)
#   rq1-m25-confirmatory-model-summary.csv      grand_mean_bet; sigma_u, sigma_s (betting rate)
#   rq2-m25-confirmatory-model-summary.csv      grand_mean_a; sigma_u, sigma_s (Beta-Binomial stake component);
#                                               sigma_v, sigma_c (all-in component)
#   rq3-m25-confirmatory-model-summary.csv     grand_mean_c; sup, sbp (positive-loss component of the hurdle-Gamma model)
#   rq1-3-participants-distribution-summary.csv median participant
#   rq1-m25-confirmatory-participants-summary.csv, rq2-..., rq3-...   classification of participants
#   rq1-m25-confirmatory-sequences-summary.csv, rq2-..., rq3-...      classification of sequences
#   The model summaries are written by the RQ tables scripts 11_rq1-tables.R,
#   21_rq2-tables.R and 31_rq3-tables.R, the participant distribution by
#   33_rq1-3-figures.R.
#
# OUTPUT
#   data/output/tables/rq1-3_model_summary.tex
#
# NOTES
#   - Entries are posterior medians; the ratio is that of the posterior medians
#     of the participant and the sequence standard deviation, on the scale of
#     the respective linear predictor
#   - Participant classes show the solid and the neutral classes (and, for
#     stakes, solid over- and under-staking) as a share of the participants of
#     the model; sequences are shown as counts of the 64 sequences
#   - The script stops if the sequence classification contains labels other
#     than those shown, or if the footnote about the three excluded normative
#     betters no longer matches the participant counts of RQ1 and RQ2
#   - The output file is overwritten at every run
# ============================================================

model_summary_table <- function(cfg) {
  
  tr <- "m25"; smp <- "confirmatory"   # FN treatment, confirmatory sample
  
  # ---- Inputs ----
  rd <- function(rq, kind) {
    f <- file.path(path_csv, paste0(rq, "-", tr, "-", smp, "-", kind, ".csv"))
    if (!file.exists(f)) stop("Missing input file: ", f, call. = FALSE)
    fread(f, encoding = "UTF-8")
  }
  ms  <- list(rd("rq1", "model-summary"), rd("rq2", "model-summary"), rd("rq3", "model-summary"))
  pr <- lapply(1:3, function(i) rd(paste0("rq", i), "participants-summary"))
  sq <- lapply(1:3, function(i) rd(paste0("rq", i), "sequences-summary"))
  f_dist <- file.path(path_csv, "rq1-3-participants-distribution-summary.csv")
  if (!file.exists(f_dist)) stop("Missing input file: ", f_dist, call. = FALSE)
  dd <- fread(f_dist, encoding = "UTF-8")[sample == "FN, confirmatory"]
  
  # one posterior summary (median, 95% interval) of a named parameter
  pm <- function(i, name) {
    x <- ms[[i]][parameter == name]
    stopifnot(nrow(x) == 1L)
    x
  }
  med <- function(q) {
    x <- dd[quantity == q]
    stopifnot(nrow(x) == 1L)
    x$median
  }
  
  # ---- Cell formatting ----
  m    <- function(s) paste0("$", s, "$")                                  # math mode
  pct1 <- function(x) sprintf("%.1f\\%%", 100 * x)                         # percent, one decimal
  sg   <- function(x) paste0(ifelse(x >= 0, "+", "-"), sprintf("%.1f", abs(100 * x)), "\\%")   # signed percent
  ci   <- function(f, lo, hi) m(paste0("[", f(lo), ", ", f(hi), "]"))      # interval
  sd2  <- function(x) paste(m(sprintf("%.2f", x$median)), m(sprintf("[%.2f, %.2f]", x$q025, x$q975)))
  mc2  <- function(s) paste0("\\multicolumn{2}{c}{", s, "}")               # cell across the two stake columns
  row  <- function(label, x, tail = "") paste0("\t\t", label, " & ", paste(x, collapse = " & "), " \\\\", tail)
  
  # ---- Model quantities ----
  gb <- pm(1, "grand_mean_bet"); ga <- pm(2, "grand_mean_a"); gc <- pm(3, "grand_mean_c")
  su <- list(pm(1, "sigma_u"), pm(2, "sigma_u"), pm(2, "sigma_v"), pm(3, "sup"))
  ss <- list(pm(1, "sigma_s"), pm(2, "sigma_s"), pm(2, "sigma_c"), pm(3, "sbp"))
  ratio <- vapply(1:4, function(k) su[[k]]$median / ss[[k]]$median, numeric(1))
  
  # ---- Classification ----
  n_par <- vapply(ms, function(x) x$n_participants[1], numeric(1))
  cnt <- function(d, col, lab) { x <- d[get(col) == lab, N]; if (length(x) == 0L) 0L else x[1] }
  shr <- function(N, n) sprintf("%.0f\\%%", 100 * N / n)
  
  seq_pre <- lapply(sq, function(d) d[classification == "preregistered"])
  stopifnot(all(seq_pre[[1]]$label == "strong"),
            all(seq_pre[[3]]$label == "strong"),
            all(seq_pre[[2]]$label %in% c("neutral", "over_moderate", "over_strong")),
            all(vapply(seq_pre, function(d) sum(d$N), numeric(1)) == ms[[1]]$n_sequences[1]))
  n_excluded <- n_par[1] - n_par[2]
  stopifnot(n_excluded == 3)   # the footnote names three excluded normative betters
  
  out <- c(
    "% Sources (data/output/csv/): rq1-m25-confirmatory-model-summary.csv, rq2-m25-confirmatory-model-summary.csv, rq3-m25-confirmatory-model-summary.csv,",
    "%   rq1-3-participants-distribution-summary.csv, rq1-m25-confirmatory-participants-summary.csv, rq2-m25-confirmatory-participants-summary.csv,",
    "%   rq3-m25-confirmatory-participants-summary.csv, rq1-m25-confirmatory-sequences-summary.csv, rq2-m25-confirmatory-sequences-summary.csv,",
    "%   rq3-m25-confirmatory-sequences-summary.csv",
    "% Produced by: 202_model-summary-table.R",
    "\\begin{table}[h]",
    "\t\\footnotesize",
    "\t\\centering",
    "\t\\begin{tabular}{lcccc}",
    "\t\t\\hline",
    "\t\t& Betting rate & \\multicolumn{2}{c}{Stake} & Welfare loss \\\\",
    "\t\t\\cline{3-4}",
    "\t\t& $\\hat{\\mu}^b$ & Deviation, $\\hat{\\mu}^a$ & All-in & $\\hat{\\mu}^c$ \\\\",
    "\t\t\\hline",
    row("EU-optimal benchmark", c(m("100\\%"), mc2(m("0\\%")), m("0\\%")), "[2pt]"),
    row("Grand mean", c(m(pct1(gb$median)), mc2(m(sg(ga$median))), m(pct1(gc$median)))),
    row("", c(ci(pct1, gb$q025, gb$q975), mc2(ci(sg, ga$q025, ga$q975)), ci(pct1, gc$q025, gc$q975))),
    row("Median participant", c(m(pct1(med("betting probability"))), mc2(m(sg(med("stake deviation")))),
                                m(pct1(med("welfare loss")))), "[4pt]"),
    row("Participant SD", vapply(su, sd2, character(1))[c(1, 2, 3, 4)]),
    row("Sequence SD", vapply(ss, sd2, character(1))[c(1, 2, 3, 4)]),
    row("Ratio", m(sprintf("%.2f", ratio)), "[4pt]"),
    "\t\t\\multirow{3}{*}{Participants}",
    paste0("\t\t& solid under: ", m(shr(cnt(pr[[1]], "underbet_label", "solid"), n_par[1])),
           " & ", mc2(paste0("solid over: ", m(shr(cnt(pr[[2]], "label", "over_solid"), n_par[2])))),
           " & solid: ", m(shr(cnt(pr[[3]], "label", "solid"), n_par[3])), " \\\\"),
    paste0("\t\t& neutral: ", m(shr(cnt(pr[[1]], "underbet_label", "neutral"), n_par[1])),
           " & ", mc2(paste0("solid under: ", m(shr(cnt(pr[[2]], "label", "under_solid"), n_par[2])))),
           " & neutral: ", m(shr(cnt(pr[[3]], "label", "neutral"), n_par[3])), " \\\\"),
    paste0("\t\t&  & ", mc2(paste0("neutral: ", m(shr(cnt(pr[[2]], "label", "neutral"), n_par[2])))), " & \\\\[4pt]"),
    "\t\t\\multirow{2}{*}{Sequences}",
    paste0("\t\t& strong under: ", m(sum(seq_pre[[1]]$N)), " & ",
           mc2(paste0("neutral: ", m(cnt(seq_pre[[2]], "label", "neutral")))),
           " & strong: ", m(sum(seq_pre[[3]]$N)), " \\\\"),
    paste0("\t\t&  & ", mc2(paste0("over-staking: ", m(sum(seq_pre[[2]][grepl("^over", label), N])))), " & \\\\[4pt]"),
    "\t\t\\hline",
    row("Likelihood", c("Bernoulli", mc2("One-inflated Beta--Binomial"), "Hurdle--Gamma")),
    row("Participants", c(m(n_par[1]), mc2(m(paste0(n_par[2], "^\\ddagger"))), m(n_par[3]))),
    row("Observations", c(m(format(ms[[1]]$n_trials[1], big.mark = "{,}", trim = TRUE)),
                          mc2(m(format(ms[[2]]$n_trials[1], big.mark = "{,}", trim = TRUE))),
                          m(format(ms[[3]]$n_trials[1], big.mark = "{,}", trim = TRUE))), "[4pt]"),
    "\t\t\\hline",
    "\t\\end{tabular}",
    "\t\\caption{Model summaries for the three behavioral quantities, FN",
    "\t\ttreatment, restricted to consistent Holt--Laury participants.",
    "\t\tPosterior medians with 95\\% credible intervals in brackets.",
    "\t\tThe median participant is the median across participants of their",
    "\t\tposterior medians, and describes the distribution of individual",
    "\t\testimates rather than a model parameter.",
    "\t\tThe participant and sequence SDs ($\\sigma_u$ and $\\sigma_s$;",
    "\t\t$\\sigma_v$ and $\\sigma_c$ in the all-in column) and their ratio, that of the posterior medians, reflect the relative",
    "\t\tcontribution of individual heterogeneity versus sequence-level",
    "\t\tvariation; values above 1 indicate that who you are matters more",
    "\t\tthan which sequence you face. The stake deviation column refers to the",
    "\t\tBeta--Binomial component of stakes below the endowment, the all-in column",
    "\t\tto the probability of staking the whole endowment, and the welfare loss",
    "\t\tcolumn to the positive-loss component of the hurdle model.",
    "\t\tClassification percentages",
    "\t\tare based on posterior concentration outside the preregistered",
    "\t\tpractical equivalence region. Sequences are classified against the",
    paste0("\t\tsame benchmark, and their counts refer to the $", ms[[1]]$n_sequences[1], "$ sequences; the two"),
    "\t\tover-staking sequences are \\texttt{HHHHHH} and \\texttt{TTTTTT}. $^\\ddagger$RQ2 excludes three normative betters",
    "\t\t(EUT prescribes a positive stake given their risk preferences)",
    "\t\twho placed too few bets in practice to identify within-participant",
    "\t\tstake variation.}",
    "\t\\label{tab:rq13_modelSummaryFN}",
    "\\end{table}")
  
  f <- file.path(path_tab, "rq1-3_model_summary.tex")
  writeLines(out, f)
  msg("Saved: ", f)
  invisible(TRUE)
}
