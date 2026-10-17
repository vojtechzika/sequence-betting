# ============================================================
# 82_ex4-canonical-rules.R
#
# PURPOSE
#   EX4: correspondence between five canonical rule-based classifications
#   of the side choice and the sequence-level directional estimates of RQ4.
#   For every betting trial (i, s) and rule k, the chosen side c_is
#   (Heads = +1, Tails = -1) is classified as hot-handish (+1) or
#   gamblerish (-1), q_is(k); trials on which a rule is undefined are
#   omitted. With x_s1,...,x_s6 the tosses of sequence s (Heads = +1),
#   the rules are
#     last   q = c_is * x_s6
#     run2   q = c_is * x_s6 if the terminal run has length >= 2
#     run3   q = c_is * x_s6 if the terminal run has length >= 3
#     imb    q = -c_is * sign(B_s), B_s = sum_t x_st; undefined if B_s = 0
#     alt    q = c_is * x_s6 if A_s <= 2, -c_is * x_s6 if A_s >= 4,
#            undefined if A_s = 3, A_s = number of alternations
#   The sequence score is Q_s(k) = mean of q_is(k) over the included betting
#   trials on s. Association with RQ4 is the Spearman correlation across
#   the 64 sequences between Q_s(k) and (a) the posterior median of mu_h_s,
#   (b) the posterior probability of Head-bias H_s(rho), (c) the posterior
#   probability of Tail-bias T_s(rho), with p-values. No posterior
#   uncertainty is propagated (81_ex4-similarity.R gives the preregistered
#   comparison with the EX1 indices, with propagation).
#
# SAMPLE
#   The full sample of RQ4 (all betting trials of a treatment), the sample
#   of the RQ4 sequence table.
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, through analysis_trials()
#   data/output/csv/rq4-<tr>-full-sequences.csv   (from 41_rq4-tables.R)
#   cfg$design$rq4$rho (first value is the main threshold)
#
# OUTPUT (data/output/csv/)
#   ex4-<tr>-canonical-seq-scores.csv    RQ4 sequence estimates with Q_s(k) and
#                                        the number of included trials n_<rule>
#   ex4-<tr>-canonical-correlations.csv  rule, outcome, Spearman rho, p-value, n
#
# NOTES
#   - The Spearman correlation is undefined (NA) for a rule whose scores are
#     constant across sequences
#   - Output files are overwritten at every run
# ============================================================

ex4_canonical_rules <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  heads  <- as.character(design$seq$side_labels$heads)
  tails  <- as.character(design$seq$side_labels$tails)
  rho_main <- as.numeric(design$rq4$rho)[1L]
  nm_main  <- gsub("\\.", "", sprintf("%.2f", rho_main))

  rules  <- c("last", "run2", "run3", "imb", "alt")

  # Sequence features: x_s1..x_s6 (Heads = +1), terminal run, imbalance, alternations
  seq_features <- function(s) {
    ch <- strsplit(s, "", fixed = TRUE)[[1L]]
    stopifnot(length(ch) == 6L, all(ch %in% c(heads, tails)))
    x  <- ifelse(ch == heads, 1L, -1L)
    list(last_out = x[6L], term_run = tail(rle(x)$lengths, 1L),
         imbalance = sum(x), alternations = sum(x[-1L] != x[-6L]))
  }

  for (tr in tr_vec) {

    f_rq4 <- file.path(path_csv, paste0("rq4-", tr, "-full-sequences.csv"))
    if (!file.exists(f_rq4)) stop("Missing ", f_rq4, "; run 41_rq4-tables.R first.", call. = FALSE)
    rq4 <- fread(f_rq4, encoding = "UTF-8")
    col_mu <- "mu_h_median"
    col_H  <- paste0("H_rho_", nm_main)
    col_T  <- paste0("T_rho_", nm_main)
    stopifnot(all(c("sequence", col_mu, col_H, col_T) %in% names(rq4)))

    d <- analysis_trials(tr, "rq4-full")
    d <- d[is.finite(stake) & stake > 0 & side %in% c(heads, tails)]
    d[, c_is := fifelse(side == heads, 1L, -1L)]

    sd_tbl <- data.table(seq = unique(as.character(d$seq)))
    ft <- lapply(sd_tbl$seq, seq_features)
    sd_tbl[, `:=`(last_out = vapply(ft, `[[`, integer(1L), "last_out"),
                  term_run = vapply(ft, `[[`, integer(1L), "term_run"),
                  imbalance = vapply(ft, `[[`, integer(1L), "imbalance"),
                  alternations = vapply(ft, `[[`, integer(1L), "alternations"))]
    d <- merge(d[, seq := as.character(seq)], sd_tbl, by = "seq")

    # Trial-level classification: +1 hot-handish, -1 gamblerish, NA undefined
    d[, q_last := c_is * last_out]
    d[, q_run2 := fifelse(term_run >= 2L, as.numeric(c_is * last_out), NA_real_)]
    d[, q_run3 := fifelse(term_run >= 3L, as.numeric(c_is * last_out), NA_real_)]
    d[, q_imb  := fifelse(imbalance != 0L, as.numeric(-c_is * sign(imbalance)), NA_real_)]
    d[, q_alt  := fcase(alternations <= 2L,  as.numeric(c_is * last_out),
                        alternations >= 4L, -as.numeric(c_is * last_out),
                        default = NA_real_)]

    # Sequence scores Q_s(k) and numbers of included trials
    qcols <- paste0("q_", rules)
    Q <- d[, lapply(.SD, function(x) if (any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_),
           .SDcols = qcols, by = seq]
    setnames(Q, qcols, paste0("Q_", rules))
    N <- d[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols = qcols, by = seq]
    setnames(N, qcols, paste0("n_", rules))
    scores <- merge(rq4, merge(Q, N, by = "seq"), by.x = "sequence", by.y = "seq", all.x = TRUE)

    outcomes <- c(col_mu, col_H, col_T)
    cors <- rbindlist(lapply(paste0("Q_", rules), function(rv) {
      rbindlist(lapply(outcomes, function(ov) {
        x <- scores[[rv]]; y <- scores[[ov]]
        ok <- is.finite(x) & is.finite(y)
        ct <- if (sum(ok) >= 3L && sd(x[ok]) > 0 && sd(y[ok]) > 0)
          suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)) else NULL
        data.table(rule = sub("^Q_", "", rv), outcome = ov,
                   rho = if (is.null(ct)) NA_real_ else unname(ct$estimate),
                   p_value = if (is.null(ct)) NA_real_ else ct$p.value,
                   n = sum(ok))
      }))
    }))
    cors[, `:=`(treatment = tr, sample = "full", rho_threshold = rho_main)]
    setcolorder(cors, c("treatment", "sample", "rule", "outcome", "rho", "p_value", "n", "rho_threshold"))

    f1 <- file.path(path_csv, paste0("ex4-", tr, "-canonical-seq-scores.csv"))
    f2 <- file.path(path_csv, paste0("ex4-", tr, "-canonical-correlations.csv"))
    fwrite(scores, f1);  msg("Saved: ", f1)
    fwrite(cors,   f2);  msg("Saved: ", f2)
  }

  invisible(TRUE)
}
