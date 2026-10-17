# ============================================================
# 101_discussion-split-sample.R
#
# PURPOSE
#   Supporting material for the Discussion: do majority followers bet against
#   the terminal run less as the imbalance grows, when the classification and
#   the outcome rest on disjoint trials? Participants are classified from side
#   choices on imbalanced sequences without a terminal run (binomial test of
#   majority-symbol choice against 0.5 at level cfg$discussion$alpha):
#   follower, restorer or uncertain. The outcome, the share of bets against
#   the terminal run, is computed on run-ending sequences by |H - T| = 0, 2,
#   4, 6.
#
# SAMPLE
#   All bets of each treatment in cfg$run$treatment.
#
# DEFINITIONS (from the six-symbol sequence string)
#   imbalance  (nH - nT) / 6, signed (positive = more heads)
#   run        last two symbols identical; run_sign = +1 (heads), -1 (tails), 0
#   abs_imb    |nH - nT| in {0, 2, 4, 6}
#   against    chosen side differs from the symbol of the run
#
# INPUT
#   data/derived/master-sequences.csv   pid, treat, seq, side
#
# OUTPUT (data/output/csv/)
#   discussion-split-sample.csv   one row per treatment, class and |H - T|:
#     participants, trials, share of bets against the run, class size
#
# NOTES
#   - Participants without bets on imbalanced no-run sequences cannot be
#     classified and are treated as uncertain
#   - Output files are overwritten at every run
# ============================================================

discussion_split_sample <- function(cfg) {
  
  alpha  <- cfg$discussion$alpha
  heads  <- as.character(cfg$design$seq$side_labels$heads)
  tails  <- as.character(cfg$design$seq$side_labels$tails)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  f_in <- file.path(path_der, "master-sequences.csv")
  stopifnot(file.exists(f_in))
  m <- fread(f_in, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "side") %in% names(m)))
  
  m[, `:=`(pid = as.character(pid), treat = as.character(treat),
           seq = as.character(seq), side = as.character(side))]
  m <- m[side %in% c(heads, tails)]                   # bets only
  
  m[, nH        := nchar(gsub("T", "", seq))]
  m[, last      := substr(seq, 6, 6)]
  m[, run       := substr(seq, 5, 5) == last]
  m[, run_sign  := fifelse(run, fifelse(last == heads, 1L, -1L), 0L)]
  m[, imbalance := (2 * nH - 6) / 6]
  m[, abs_imb   := abs(2 * nH - 6)]
  m[, chose_H   := as.integer(side == heads)]
  
  out <- list()
  for (tr in tr_vec) {
    d_all <- m[treat == tr]
    if (nrow(d_all) == 0L) next
    
    # classification on imbalanced sequences without a terminal run
    nr <- d_all[run_sign == 0L & abs_imb > 0L]
    nr[, maj := as.integer((chose_H == 1L & imbalance > 0) | (chose_H == 0L & imbalance < 0))]
    g <- nr[, .(k = sum(maj), n = .N), by = pid]
    g[, p   := mapply(function(k, n) stats::binom.test(k, n, 0.5)$p.value, k, n)]
    g[, cls := fifelse(p < alpha, fifelse(k / n > 0.5, "follower", "restorer"), "uncertain")]
    
    # outcome on run-ending sequences
    r <- d_all[run_sign != 0L]
    r[, against := as.integer((chose_H == 0L & run_sign == 1L) | (chose_H == 1L & run_sign == -1L))]
    r <- merge(r, g[, .(pid, cls)], by = "pid", all.x = TRUE)
    r[is.na(cls), cls := "uncertain"]
    
    size <- unique(r[, .(pid, cls)])[, .(class_size = .N), by = cls]
    res  <- r[, .(n_participants = uniqueN(pid), n_trials = .N,
                  share_against_run = mean(against)), by = .(cls, abs_imb)]
    res  <- merge(res, size, by = "cls")
    res[, treatment := tr]
    out[[tr]] <- res
  }
  
  tbl <- rbindlist(out)
  setorder(tbl, treatment, cls, abs_imb)
  
  f <- file.path(path_csv, "discussion-split-sample.csv")
  fwrite(tbl, f);  msg("Saved: ", f)
  invisible(tbl)
}