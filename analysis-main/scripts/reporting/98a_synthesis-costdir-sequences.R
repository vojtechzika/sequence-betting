# ============================================================
# 98a_synthesis-costdir-sequences.R
#
# PURPOSE
#   Synthesis, sequence level: are the sequences that are economically
#   costly the same as those whose side choice is steered by the sequence?
#   The 64 sequences are the population of interest, so posterior
#   uncertainty is the only uncertainty. Every quantity is computed within
#   posterior draws and summarised by the posterior median and 95% credible
#   interval.
#     A1  Link between side bias |mu_h_s - hbar| (RQ4) and welfare loss mu_c_s
#         (RQ3): per draw, the Spearman correlation (all 64 sequences, and
#         excluding the two pure sequences) and the OLS slope of loss on bias,
#         alone, with imbalance-class dummies, and with imbalance class and
#         terminal run.
#     A2  Feature regressions: side bias, betting rate (RQ1), stake deviation
#         (RQ2) and welfare loss each regressed on a terminal-run indicator
#         and imbalance-class dummies (class dummies rather than a linear
#         term, because loss and stake effects are concentrated in the most
#         imbalanced sequences).
#
# SAMPLE
#   Treatments with a confirmatory sample (cfg$design$a_flags$betting_normative):
#   RQ1-RQ3 from the confirmatory fits; other treatments: full fits. RQ4 is
#   always from the full fit.
#
# INPUT (models selected under the keys rq1, rq2, rq3, rq4)
#   fit-, pid-levels-, seq-levels-<tr>-<sample>.rds: mu_b, mu_a, mu_c (sequence
#   means), mu_h and hbar (RQ4)
#
# OUTPUT
#   data/output/csv/synthesis-<tr>-costdir-seq-link.csv
#   data/output/csv/synthesis-<tr>-costdir-seq-features.csv
#   data/output/tables/synthesis-<tr>-costdir-features.tex
#
# NOTES
#   - Each model is estimated separately, so the K draws (cfg$synthesis$costdir)
#     are subsampled independently from each fit and paired by row index;
#     within a model the draws stay joint across sequences
#   - Imbalance class: |number of heads - 3|, labelled 3-3, 2-4, 1-5, 0-6;
#     the 0-6 class contains only the two pure sequences
#   - Seeded by cfg$run$seed
#   - Output files are overwritten at every run
# ============================================================

synthesis_costdir_sequences <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  set.seed(as.integer(cfg$run$seed))
  tl <- synthesis_draw_tools(cfg)
  K  <- tl$K
  
  n_H <- function(s) nchar(gsub("T", "", s))
  seq_features <- function(s) data.table(
    sequence  = s,
    imb_class = factor(abs(n_H(s) - 3), levels = 0:3, labels = c("3-3", "2-4", "1-5", "0-6")),
    run       = as.integer(substr(s, 5, 5) == substr(s, 6, 6)))
  
  out <- list()
  
  for (tr in tr_vec) {
    
    conf_ok <- isTRUE(cfg$design$a_flags$betting_normative[[tr]])
    tg_main <- if (conf_ok) "confirmatory" else "full"
    
    R4 <- tl$load_rq4(tr)
    M  <- list(rq1 = tl$load_rq("rq1", tr, tg_main, "mu_b", "mu_b_i"),
               rq2 = tl$load_rq("rq2", tr, tg_main, "mu_a", "mu_a_i"),
               rq3 = tl$load_rq("rq3", tr, tg_main, "mu_c", "mu_c_i"))
    if (is.null(R4) || any(vapply(M, is.null, logical(1L)))) {
      warning("Synthesis cost vs direction (sequences): fits missing for ", tr, "; skipped.",
              call. = FALSE)
      next
    }
    
    # ---- A1: side bias versus welfare loss ----
    seqs   <- sort(colnames(R4$seq))
    Fe     <- seq_features(seqs)
    bias   <- 100 * abs(R4$seq[, seqs] - R4$hbar)          # K x 64, percentage points
    bet    <- 100 * M$rq1$seq[, seqs]
    stake  <- 100 * M$rq2$seq[, seqs]
    loss   <- 100 * M$rq3$seq[, seqs]
    nopure <- Fe$imb_class != "0-6"
    
    X_alone  <- matrix(1, length(seqs), 1L, dimnames = list(NULL, "(Intercept)"))
    X_imb    <- model.matrix(~ imb_class, Fe)
    X_imbrun <- model.matrix(~ imb_class + run, Fe)
    slope <- function(Xbase, k) {
      X <- cbind(Xbase, bias = bias[k, ])
      qr.coef(qr(X), loss[k, ])[["bias"]]
    }
    a1_draws <- data.table(
      rho_all    = sapply(1:K, function(k) cor(bias[k, ], loss[k, ], method = "spearman")),
      rho_nopure = sapply(1:K, function(k) cor(bias[k, nopure], loss[k, nopure], method = "spearman")),
      b_alone    = sapply(1:K, function(k) slope(X_alone, k)),
      b_imb      = sapply(1:K, function(k) slope(X_imb, k)),
      b_imbrun   = sapply(1:K, function(k) slope(X_imbrun, k)))
    a1 <- rbindlist(list(
      tl$sumrow(a1_draws$rho_all,    treatment = tr, statistic = "spearman_rho", spec = "all 64"),
      tl$sumrow(a1_draws$rho_nopure, treatment = tr, statistic = "spearman_rho", spec = "excl. pure"),
      tl$sumrow(a1_draws$b_alone,    treatment = tr, statistic = "ols_slope_loss_on_bias", spec = "alone"),
      tl$sumrow(a1_draws$b_imb,      treatment = tr, statistic = "ols_slope_loss_on_bias", spec = "+ imbalance class"),
      tl$sumrow(a1_draws$b_imbrun,   treatment = tr, statistic = "ols_slope_loss_on_bias", spec = "+ imbalance class + run"),
      tl$sumrow(10 * a1_draws$b_imb, treatment = tr, statistic = "loss_pp_per_10pp_bias", spec = "+ imbalance class")))
    a1[, `:=`(sample_bet_stake_loss = tg_main, n_draws = K)]
    f1 <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-seq-link.csv"))
    fwrite(a1, f1);  msg("Saved: ", f1)
    
    # ---- A2: feature regressions ----
    X_feat <- model.matrix(~ run + imb_class, Fe)
    qrF    <- qr(X_feat)
    Y      <- list(bias = bias, bet = bet, stake = stake, loss = loss)
    a2 <- rbindlist(lapply(names(Y), function(y) {
      B <- qr.coef(qrF, t(Y[[y]]))                          # p x K
      rbindlist(lapply(rownames(B), function(t)
        tl$sumrow(B[t, ], treatment = tr, outcome = y, term = t)))
    }))
    a2[, sample_bet_stake_loss := tg_main]
    f2 <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-seq-features.csv"))
    fwrite(a2, f2);  msg("Saved: ", f2)
    
    # ---- LaTeX table of A2 ----
    terms  <- c("run", "imb_class2-4", "imb_class1-5", "imb_class0-6", "(Intercept)")
    labels <- c("Terminal run", "Imbalance 2--4", "Imbalance 1--5", "Imbalance 0--6",
                "Constant (balanced, no run)")
    fmt1 <- function(x) sprintf("%.1f", ifelse(abs(x) < 0.05, 0, x))
    body <- unlist(lapply(seq_along(terms), function(i) {
      cs <- lapply(names(Y), function(y) tl$cell(a2[outcome == y & term == terms[i]], fmt1))
      c(paste(labels[i], "&", paste(sapply(cs, `[`, 1), collapse = " & "), "\\\\"),
        paste(" &", paste(sapply(cs, `[`, 2), collapse = " & "),
              if (i < length(terms)) "\\\\[3pt]" else "\\\\"))
    }))
    tl$write_tex(
      file.path(path_tab, paste0("synthesis-", tr, "-costdir-features.tex")),
      caption = "Sequence features and outcomes across the four margins",
      label   = paste0("Tab:synthesis-features", if (tr != tr_vec[1L]) paste0("-", tr)),
      colspec = "l*{4}{c}",
      header  = c(" & Direction & \\multicolumn{3}{c}{Engagement and welfare} \\\\",
                  "\\cmidrule(lr){2-2}\\cmidrule(lr){3-5}",
                  " & Side bias & Betting rate & Stake deviation & Welfare loss \\\\",
                  " & $|\\hat\\mu^h_s-\\bar h|$ & $\\hat\\mu^b_s$ & $\\hat\\mu^a_s$ & $\\hat\\mu^c_s$ \\\\"),
      body    = body,
      footer  = c("Sequences & 64 & 64 & 64 & 64 \\\\",
                  paste("Posterior draws &", paste(rep(K, 4), collapse = " & "), "\\\\")),
      note    = paste0("Each column regresses a sequence-level quantity on a terminal-run indicator ",
                       "(last two symbols identical) and imbalance-class indicators (reference: balanced 3--3), ",
                       "separately within each posterior draw. Cells report the posterior median of the coefficient ",
                       "(percentage points) with its 95\\% credible interval in brackets; an asterisk denotes an ",
                       "interval excluding zero. Side bias uses the full sample; betting rate, stake deviation, and ",
                       "welfare loss (shares of the endowment) use the ",
                       if (conf_ok) "confirmatory" else "full", " sample. ",
                       "The 0--6 class contains only the two pure sequences."))
    
    out[[tr]] <- list(seq_link = a1, seq_features = a2)
  }
  
  invisible(out)
}
