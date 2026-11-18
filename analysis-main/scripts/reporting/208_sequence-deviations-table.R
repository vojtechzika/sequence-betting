# ============================================================
# 208_sequence-deviations-table.R
#
# PURPOSE
#   Table of the deviations of the sequence-level estimates from their grand
#   means for all 64 sequences, one table per treatment, for the online
#   appendix. Each cell is the posterior median of the paired difference
#   d_s = mu_s - grand mean, in percentage points, with its 95% credible
#   interval, for the bet probability (RQ1), the staking deviation (RQ2),
#   the welfare loss (RQ3), and the Heads probability (RQ4). Cells whose
#   interval excludes zero are set in bold.
#
# SAMPLE
#   FN: confirmatory sample for RQ1 to RQ3, full sample for RQ4.
#   FP: full sample for all four (the FP treatment has no confirmatory sample).
#
# INPUT (data/output/csv/)
#   synthesis-<treatment>-<sample>-sequence-diffs.csv   (91_synthesis-sequences.R)
#   columns <e>_d_median, <e>_d_lo, <e>_d_hi, <e>_sig for e = b, a, c, h
#   The Heads probability (e = h) of the confirmatory run comes from the full
#   RQ4 fit, as set in 91_synthesis-sequences.R.
#
# OUTPUT (data/output/tables/)
#   sequence-deviations-fn.tex, sequence-deviations-fp.tex   longtable, one row per sequence
#
# NOTES
#   - The grand means are the first rows of sequence-estimates-fn.tex and
#     sequence-estimates-fp.tex (207_sequence-estimates-table.R)
#   - Sequences are listed alphabetically (H before T)
#   - Output files are overwritten at every run
# ============================================================

sequence_deviations_table <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  style  <- treatment_style(cfg)

  # estimand prefix in the synthesis files, in the order of the columns (RQ1 to RQ4)
  est <- c("b", "a", "c", "h")

  fmt <- function(m, lo, hi, sig) {
    f <- function(x) { x <- round(100 * x, 1); x[x == 0] <- 0; gsub("-", "$-$", sprintf("%.1f", x)) }
    cell <- paste0(f(m), " [", f(lo), ", ", f(hi), "]")
    ifelse(sig, paste0("\\textbf{", cell, "}"), cell)
  }

  for (k in seq_along(tr_vec)) {
    tr  <- tr_vec[k]
    smp <- if (k == 1L && isTRUE(cfg$design$a_flags$betting_normative[[tr]])) "confirmatory" else "full"
    f   <- file.path(path_csv, paste0("synthesis-", tr, "-", smp, "-sequence-diffs.csv"))
    stopifnot(file.exists(f))
    d <- fread(f, encoding = "UTF-8")
    stopifnot(nrow(d) == 64L)
    setorder(d, sequence)   # alphabetical, H before T

    tbl <- data.table(sequence = d$sequence)
    for (e in est) {
      tbl[[e]] <- fmt(d[[paste0(e, "_d_median")]], d[[paste0(e, "_d_lo")]],
                      d[[paste0(e, "_d_hi")]], as.logical(d[[paste0(e, "_sig")]]))
    }

    lab    <- if (k == 1L) "fn" else "fp"
    header <- paste("Sequence & Bet Decision (RQ1) & Staking Deviation (RQ2) & Welfare Loss (RQ3) & Side Choice (RQ4) \\\\",
                    " & $\\Delta^b_s$ & $\\Delta^a_s$ & $\\Delta^c_s$ & $\\Delta^h_s$ \\\\", sep = "\n")
    head <- paste(
      "\\begin{longtable}{@{}>{\\raggedright\\arraybackslash}p{0.105\\textwidth}*{4}{>{\\centering\\arraybackslash}p{0.2055\\textwidth}}@{}}",
      paste0("\\caption{Deviations of the sequence-level estimates from the grand mean, ", style$type[k], " treatment. The deviation $\\Delta_s$ is the posterior median of the difference between the estimate of sequence $s$ and the grand mean of the same posterior draw, in percentage points, with the 95\\% credible interval, for the bet probability (RQ1), the staking deviation (RQ2), the welfare loss (RQ3), and the \\textit{Heads} probability (RQ4). Intervals that exclude zero are set in bold. ",
             if (k == 1L) "RQ1 to RQ3 use the confirmatory sample and RQ4 the full sample. " else "All estimates use the full sample. ",
             "The grand means are given in Table~\\ref{Tab:sequence-estimates-", lab, "}. The sequences are ordered with H before T.}"),
      paste0("\\label{Tab:sequence-deviations-", lab, "} \\\\"),
      "\\toprule", header, "\\midrule", "\\endfirsthead",
      "\\toprule", header, "\\midrule", "\\endhead",
      "\\bottomrule", "\\endfoot", "\\bottomrule", "\\endlastfoot",
      sep = "\n")

    rows <- apply(tbl, 1L, function(r) paste0("\\texttt{", r[1L], "} & ", paste(r[-1L], collapse = " & "), " \\\\"))

    out <- file.path(path_tab, paste0("sequence-deviations-", lab, ".tex"))
    writeLines(c("\\begingroup", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}", "\\setlength{\\LTcapwidth}{\\textwidth}", head, rows, "\\end{longtable}", "\\endgroup"), out)
    msg("Saved: ", out)
  }
  invisible(TRUE)
}
