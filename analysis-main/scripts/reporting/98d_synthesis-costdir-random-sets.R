# ============================================================
# 98d_synthesis-costdir-random-sets.R
#
# PURPOSE
#   Robustness check for the contrast "at least one strategy vs. no strategy"
#   (98c_synthesis-costdir-participants.R). Every participant faces the same
#   exhaustive set of 64 six-toss sequences. This script asks whether the
#   welfare-loss gap between participants with and without a credible cue
#   slope depends on that set. For random subsets of k sequences (the same
#   subset for all participants), it computes the gap in mean welfare loss
#   and summarizes its distribution over subsets.
#
# WHAT IT DOES
#   1. Builds the participant x sequence matrix of trial-level welfare loss
#      of RQ3 (percentage points of the endowment) for the confirmatory sample.
#   2. Defines "strategy" as a credible personal slope (95% interval excluding
#      zero) on the terminal run or on imbalance in the side-choice model;
#      "none" as uncertain on both.
#   3. For k in {8, 16, 32} draws n_sets random subsets of k sequences and
#      records the gap (strategy minus none) in mean loss; k = 64 is the gap
#      on all sequences.
#
# INPUT
#   data/output/csv/synthesis-<tr>-joint-run-slope-participants.csv  (95_synthesis-joint-heterogeneity.R)
#   data/output/csv/synthesis-<tr>-joint-imb-slope-participants.csv  (95_synthesis-joint-heterogeneity.R)
#   data/output/csv/rq2-<tr>-confirmatory-participants.csv           (participants of the contrast in 98c)
#   scripts/modelling/_rq3-welfare-loss.R                            (trial-level welfare loss, as in RQ3)
#
# OUTPUT
#   data/output/csv/synthesis-<tr>-costdir-random-sets.csv
#     columns: treatment, sample, k, n_sets, n_strategy, n_none,
#              mean_gap, q025, q975, share_gt0
#
# NOTES
#   - Only treatments in which betting is normative (cfg$design$a_flags$betting_normative)
#     have a confirmatory sample; others are skipped.
#   - The gap is a plain difference of participant means, not the
#     posterior/Bayesian-bootstrap contrast of 98c; its value on all 64
#     sequences therefore differs somewhat from the 98c estimate.
#   - The slope classes are estimated on all 64 sequences, so classification
#     and outcome share trials, as in 98c.
#   - Seeded by cfg$run$seed; the output file is overwritten at every run.
# ============================================================

synthesis_costdir_random_sets <- function(cfg) {
  
  source(here::here("scripts", "modelling", "_rq3-welfare-loss.R"))
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  set.seed(as.integer(cfg$run$seed))
  ks     <- c(8L, 16L, 32L)    # subset sizes
  n_sets <- 5000L              # random subsets per size
  
  for (tr in tr_vec) {
    
    if (!isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
    
    f_run <- file.path(path_csv, paste0("synthesis-", tr, "-joint-run-slope-participants.csv"))
    f_imb <- file.path(path_csv, paste0("synthesis-", tr, "-joint-imb-slope-participants.csv"))
    f_rq2 <- file.path(path_csv, paste0("rq2-", tr, "-confirmatory-participants.csv"))
    if (!all(file.exists(f_run, f_imb, f_rq2))) {
      warning("Random sets (", tr, "): input files missing; skipped.", call. = FALSE)
      next
    }
    
    # ---- Slope classes and participants of the contrast ----
    run <- fread(f_run, encoding = "UTF-8", colClasses = c(pid = "character"))
    imb <- fread(f_imb, encoding = "UTF-8", colClasses = c(pid = "character"))
    rq2 <- fread(f_rq2, encoding = "UTF-8", colClasses = c(pid = "character"))
    
    # ---- Trial-level welfare loss, one value per participant and sequence ----
    wl <- rq3_welfare_loss(cfg, tr, "confirmatory")
    L  <- data.table(pid = wl$d$pid, seq = wl$d$seq, loss = 100 * wl$y)
    stopifnot(!anyDuplicated(L[, .(pid, seq)]))
    M  <- dcast(L, pid ~ seq, value.var = "loss")
    
    P <- Reduce(intersect, list(M$pid, run$pid, imb$pid, rq2$pid))
    M <- M[match(P, pid)]
    X <- as.matrix(M[, -"pid"])
    stopifnot(!anyNA(X))
    
    strat <- (run$class[match(P, run$pid)] != "uncertain") |
             (imb$class[match(P, imb$pid)] != "uncertain")
    
    # ---- Gap in mean loss (strategy minus none) over random subsets ----
    gap <- function(cols) {
      v <- rowMeans(X[, cols, drop = FALSE])
      mean(v[strat]) - mean(v[!strat])
    }
    rows <- list()
    rows[["64"]] <- data.table(k = ncol(X), n_sets = 1L, mean_gap = gap(seq_len(ncol(X))),
                               q025 = NA_real_, q975 = NA_real_, share_gt0 = NA_real_)
    for (k in ks) {
      g <- replicate(n_sets, gap(sample.int(ncol(X), k)))
      rows[[as.character(k)]] <- data.table(k = k, n_sets = n_sets, mean_gap = mean(g),
                                            q025 = unname(quantile(g, .025)),
                                            q975 = unname(quantile(g, .975)),
                                            share_gt0 = mean(g > 0))
    }
    res <- rbindlist(rows)
    res[, `:=`(treatment = tr, sample = "confirmatory",
               n_strategy = sum(strat), n_none = sum(!strat))]
    setcolorder(res, c("treatment", "sample", "k", "n_sets", "n_strategy", "n_none"))
    
    f_out <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-random-sets.csv"))
    fwrite(res, f_out)
    cat("Saved:", f_out, "\n")
  }
  
  invisible(NULL)
}
