# ============================================================
# 207_sequence-estimates-table.R
#
# PURPOSE
#   Table of the sequence-level estimates of RQ1 to RQ4 for all 64 sequences,
#   one table per treatment, for the online appendix. Each cell is the
#   posterior mean of the sequence-level quantity with its 95% credible
#   interval: the bet probability (RQ1), the staking deviation (RQ2), the
#   welfare loss (RQ3), and the Heads probability (RQ4).
#
# SAMPLE
#   FN: confirmatory sample for RQ1 to RQ3, full sample for RQ4.
#   FP: full sample for all four (the FP treatment has no confirmatory sample).
#
# INPUT (data/output/csv/)
#   rq<k>-<treatment>-<sample>-sequences.csv   k = 1, ..., 4
#   columns mu_b_*, mu_a_*, mu_c_*, mu_h_*: mean, q025, q975
#   rq<k>-<treatment>-<sample>-model-summary.csv   grand_mean_bet, grand_mean_a,
#   grand_mean_c, hbar (mean, q025, q975)
#
# OUTPUT (data/output/tables/)
#   sequence-estimates-fn.tex, sequence-estimates-fp.tex   longtable, one row per sequence
#
# NOTES
#   - The first row is the grand mean (RQ1 to RQ3) and the baseline h-bar (RQ4)
#   - Sequences are listed alphabetically (H before T)
#   - Output files are overwritten at every run
# ============================================================

sequence_estimates_table <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  style  <- treatment_style(cfg)

  # research question, parameter prefix in the sequence files, and the
  # parameter of the model summary used as the reference row (grand mean of
  # RQ1 to RQ3, baseline Heads probability of RQ4)
  qs <- list(
    list(rq = "rq1", par = "mu_b", ref = "grand_mean_bet"),
    list(rq = "rq2", par = "mu_a", ref = "grand_mean_a"),
    list(rq = "rq3", par = "mu_c", ref = "grand_mean_c"),
    list(rq = "rq4", par = "mu_h", ref = "hbar")
  )

  fmt <- function(m, lo, hi) {
    f <- function(x) { x <- round(x, 2); x[x == 0] <- 0; gsub("-", "$-$", sprintf("%.2f", x)) }
    paste0(f(m), " [", f(lo), ", ", f(hi), "]")
  }

  for (k in seq_along(tr_vec)) {
    tr <- tr_vec[k]
    # confirmatory sample of RQ1 to RQ3 exists for the first treatment (FN) only
    smp <- function(rq) if (k == 1L && rq != "rq4") "confirmatory" else "full"

    cols <- lapply(qs, function(q) {
      f <- file.path(path_csv, paste0(q$rq, "-", tr, "-", smp(q$rq), "-sequences.csv"))
      stopifnot(file.exists(f))
      d <- fread(f, encoding = "UTF-8")
      d[, .(sequence, cell = fmt(get(paste0(q$par, "_mean")),
                                 get(paste0(q$par, "_q025")),
                                 get(paste0(q$par, "_q975"))))]
    })
    tbl <- Reduce(function(a, b) merge(a, b, by = "sequence", sort = FALSE), cols)
    stopifnot(nrow(tbl) == 64L)
    setorder(tbl, sequence)   # alphabetical, H before T
    setnames(tbl, c("sequence", vapply(qs, `[[`, "", "rq")))

    lab  <- if (k == 1L) "fn" else "fp"
    head <- paste(
      "\\begin{longtable}{@{}>{\\raggedright\\arraybackslash}p{0.105\\textwidth}*{4}{>{\\centering\\arraybackslash}p{0.2055\\textwidth}}@{}}",
      paste0("\\caption{Sequence-level estimates, ", style$type[k], " treatment. Posterior mean with 95\\% credible interval of the bet probability (RQ1), the staking deviation (RQ2), the welfare loss (RQ3), and the \\textit{Heads} probability (RQ4) for each sequence. ",
             if (k == 1L) "RQ1 to RQ3 use the confirmatory sample and RQ4 the full sample. " else "All estimates use the full sample. ",
             "The first row gives the grand mean of RQ1 to RQ3 and the baseline Heads probability $\\bar h$ of RQ4. The sequences are ordered with H before T.}"),
      paste0("\\label{Tab:sequence-estimates-", lab, "} \\\\"),
      "\\toprule",
      "Sequence & Bet Decision (RQ1) & Staking Deviation (RQ2) & Welfare Loss (RQ3) & Side Choice (RQ4) \\\\",
      "\\midrule",
      "\\endfirsthead",
      "\\toprule",
      "Sequence & Bet Decision (RQ1) & Staking Deviation (RQ2) & Welfare Loss (RQ3) & Side Choice (RQ4) \\\\",
      "\\midrule",
      "\\endhead",
      "\\bottomrule",
      "\\endfoot",
      "\\bottomrule",
      "\\endlastfoot",
      sep = "\n")
    # reference row from the model summaries
    ref <- vapply(qs, function(q) {
      m <- fread(file.path(path_csv, paste0(q$rq, "-", tr, "-", smp(q$rq), "-model-summary.csv")),
                 encoding = "UTF-8")[parameter == q$ref]
      stopifnot(nrow(m) == 1L)
      fmt(m$mean, m$q025, m$q975)
    }, character(1))
    ref_row <- c(paste0("\\emph{Grand mean} & ", paste(ref, collapse = " & "), " \\\\"), "\\midrule")

    rows <- apply(tbl, 1L, function(r) paste0("\\texttt{", r[1L], "} & ", paste(r[-1L], collapse = " & "), " \\\\"))

    out <- file.path(path_tab, paste0("sequence-estimates-", lab, ".tex"))
    writeLines(c("\\begingroup", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}", "\\setlength{\\LTcapwidth}{\\textwidth}", head, ref_row, rows, "\\end{longtable}", "\\endgroup"), out)
    msg("Saved: ", out)
  }
  invisible(TRUE)
}
