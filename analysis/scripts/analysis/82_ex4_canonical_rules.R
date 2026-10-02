# ============================================================
# 82_ex4_canonical_rules.R
#
# PURPOSE
#   For each sequence, compute the five canonical rule-based
#   scores from the preregistration (EX4):
#     (i)   last-outcome rule
#     (ii)  run-length-conditioned rule (terminal run >= 2)
#     (iii) run-length-conditioned rule (terminal run >= 3)
#     (iv)  imbalance rule
#     (v)   alternation rule
#
#   Q_s^(k) = mean_i(q_is^(k)) across betting observations
#   for sequence s under rule k.
#
#   Merges with rq4_<tr>_sequences.csv and correlates
#   Q_s^(k) with mu_h_median, H_delta_005, T_delta_005.
#
# INPUTS
#   path_src/master_sequences.csv
#   path_out/rq4_<tr>_sequences.csv
#
# OUTPUTS
#   path_out/ex4_<tr>_canonical_seq_scores.csv
#   path_out/ex4_<tr>_canonical_correlations.csv
#
# CALL ORDER
#   rq4_tables(cfg)          -- 43_
#   ex4_canonical_rules(cfg) -- this script
# ============================================================

library(data.table)

ex4_canonical_rules <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(nzchar(lab_heads), nzchar(lab_tails))
  
  delta_vec  <- as.numeric(design$rq4$delta %||% c(0.05, 0.03, 0.08))
  delta_main <- delta_vec[1]
  nm_main    <- gsub("\\.", "", sprintf("%.2f", delta_main))
  
  # ----------------------------
  # Load master
  # ----------------------------
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  master <- fread(infile, encoding = "UTF-8")
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  master[, side  := as.character(side)]
  master[is.na(stake), stake := 0]
  
  # ----------------------------
  # Sequence feature extraction
  # ----------------------------
  seq_features <- function(s) {
    chars <- strsplit(s, "")[[1]]
    x     <- ifelse(chars == "H", 1L, -1L)   # H=+1, T=-1
    n     <- length(x)
    
    rle_out      <- rle(x)
    term_run     <- tail(rle_out$lengths, 1)
    last_out     <- x[n]
    imbalance    <- sum(x)                     # B_s
    alternations <- sum(x[-1] != x[-n])        # A_s
    
    list(last_out     = last_out,
         term_run     = term_run,
         imbalance    = imbalance,
         alternations = alternations)
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ----------------------------
    # Load RQ4 sequence estimates
    # ----------------------------
    f_rq4 <- file.path(path_out, paste0("rq4_", tr, "_sequences.csv"))
    if (!file.exists(f_rq4)) {
      warning("ex4_canonical_rules: missing RQ4 sequences for tr='", tr,
              "'. Run rq4_tables() first. Skipping.")
      next
    }
    rq4 <- fread(f_rq4)
    
    # Outcome columns for correlations
    col_mu  <- "mu_h_median"
    col_H   <- paste0("H_delta_", nm_main)
    col_T   <- paste0("T_delta_", nm_main)
    stopifnot(all(c(col_mu, col_H, col_T) %in% names(rq4)))
    
    # ----------------------------
    # Subset to treatment betting trials
    # ----------------------------
    dt <- master[
      treat == tr &
        is.finite(stake) & stake > 0 &
        side %in% c(lab_heads, lab_tails)
    ]
    
    # Encode chosen side: H = +1, T = -1
    dt[, c_is := fifelse(side == lab_heads, 1L, -1L)]
    
    # ----------------------------
    # Attach sequence features
    # ----------------------------
    seq_dict <- unique(dt[, .(seq)])
    feat_list <- lapply(seq_dict$seq, seq_features)
    seq_dict[, last_out     := sapply(feat_list, `[[`, "last_out")]
    seq_dict[, term_run     := sapply(feat_list, `[[`, "term_run")]
    seq_dict[, imbalance    := sapply(feat_list, `[[`, "imbalance")]
    seq_dict[, alternations := sapply(feat_list, `[[`, "alternations")]
    
    dt <- seq_dict[dt, on = "seq"]
    
    # ----------------------------
    # Trial-level rule scores q_is^(k)
    # +1 = hot-handish, -1 = gamblerish, NA = undefined
    # ----------------------------
    
    # (i) Last-outcome: q = c_is * x_s6
    dt[, q_last := c_is * last_out]
    
    # (ii) Run >= 2
    dt[, q_run2 := fifelse(term_run >= 2L, c_is * last_out, NA_real_)]
    
    # (iii) Run >= 3
    dt[, q_run3 := fifelse(term_run >= 3L, c_is * last_out, NA_real_)]
    
    # (iv) Imbalance: q = -c_is * sign(B_s); undefined if B_s == 0
    dt[, q_imb := fifelse(imbalance != 0L,
                          -c_is * sign(imbalance),
                          NA_real_)]
    
    # (v) Alternation rule
    dt[, q_alt := fcase(
      alternations <= 2L,  as.numeric(c_is * last_out),
      alternations >= 4L, -as.numeric(c_is * last_out),
      default = NA_real_
    )]
    
    # ----------------------------
    # Sequence-level scores Q_s^(k)
    # Q_s^(k) = mean_i(q_is^(k)), NA removed
    # ----------------------------
    rules     <- c("q_last", "q_run2", "q_run3", "q_imb", "q_alt")
    rule_lbls <- c("Q_last", "Q_run2", "Q_run3", "Q_imb", "Q_alt")
    
    seq_scores <- dt[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                     .SDcols = rules, by = seq]
    setnames(seq_scores, rules, rule_lbls)
    
    # N available per rule per sequence
    seq_n <- dt[, lapply(.SD, function(x) sum(!is.na(x))),
                .SDcols = rules, by = seq]
    setnames(seq_n, rules,
             paste0("n_", c("last", "run2", "run3", "imb", "alt")))
    
    seq_scores <- seq_n[seq_scores, on = "seq"]
    
    # ----------------------------
    # Merge with RQ4 estimates
    # ----------------------------
    merged <- rq4[seq_scores, on = c(sequence = "seq")]
    
    # ----------------------------
    # Correlations (Spearman)
    # ----------------------------
    outcome_vars <- c(col_mu, col_H, col_T)
    
    cor_results <- rbindlist(lapply(rule_lbls, function(rv) {
      rbindlist(lapply(outcome_vars, function(ov) {
        ct <- cor.test(merged[[rv]], merged[[ov]],
                       method = "spearman", exact = FALSE)
        data.table(
          rule      = rv,
          outcome   = ov,
          rho       = round(ct$estimate, 3),
          p_value   = round(ct$p.value,  4),
          n         = sum(!is.na(merged[[rv]]) & !is.na(merged[[ov]]))
        )
      }))
    }))
    
    # ----------------------------
    # Print
    # ----------------------------
    cat("\n=== EX4 canonical rule correlations (", tr, ") ===\n\n", sep = "")
    
    mat <- matrix(
      NA_real_,
      nrow = length(rule_lbls),
      ncol = length(outcome_vars),
      dimnames = list(rule_lbls, outcome_vars)
    )
    pmat <- mat
    for (rv in rule_lbls)
      for (ov in outcome_vars) {
        mat[rv,  ov] <- cor_results[rule == rv & outcome == ov, rho]
        pmat[rv, ov] <- cor_results[rule == rv & outcome == ov, p_value]
      }
    
    cat("Spearman rho:\n");  print(round(mat,  3))
    cat("\np-values:\n");    print(round(pmat, 4))
    
    # ----------------------------
    # Save
    # ----------------------------
    f_scores <- file.path(path_out,
                          paste0("ex4_", tr, "_canonical_seq_scores.csv"))
    f_cor    <- file.path(path_out,
                          paste0("ex4_", tr, "_canonical_correlations.csv"))
    
    if (!should_skip(f_scores, cfg, "output",
                     paste0("EX4 canonical seq scores (", tr, ")"))) {
      fwrite(merged,      f_scores)
      msg("Saved: ", f_scores)
    }
    if (!should_skip(f_cor, cfg, "output",
                     paste0("EX4 canonical correlations (", tr, ")"))) {
      fwrite(cor_results, f_cor)
      msg("Saved: ", f_cor)
    }
    
    outputs[[tr]] <- list(
      canonical_seq_scores_csv  = f_scores,
      canonical_correlations_csv = f_cor
    )
  }
  
  invisible(outputs)
}