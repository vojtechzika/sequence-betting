# ============================================================
# 101_discussion.R
#
# PURPOSE
#   Supporting material for the Discussion (Sec:Dis): do majority-followers
#   bet against the terminal run less as imbalance grows, when the
#   classification and the outcome rest on DISJOINT trials?
#
# METHOD
#   Participants are classified from side choices on imbalanced sequences
#   WITHOUT a terminal run (binomial test of majority-symbol choice
#   against 0.5, p < .05): follower / restorer / uncertain. The outcome,
#   the share of bets against the terminal run, is computed on run-ending
#   sequences by |H-T| = 0, 2, 4, 6.
#
# DEFINITIONS (from the 6-symbol sequence string itself)
#   imbalance : (nH - nT) / 6, signed (positive = more heads)
#   run       : last two symbols identical; run_sign = +1 (H), -1 (T), 0
#   abs_imb   : |nH - nT| in {0, 2, 4, 6}
#   against   : chosen side differs from the run's symbol
#
# INPUT
#   path_src/master_sequences.csv  (pid, treat, seq, side)
#     side: "H", "T", or "NB" (no bet)
#
# OUTPUT
#   path_out/discussion_split_sample.csv -- one row per treatment x class
#     x |H-T|: participants, trials, share of bets against the run, and
#     class size
#
# NOTES
#   - Participants without bets on imbalanced no-run sequences cannot be
#     classified and are treated as uncertain.
#   - Only console output: "Saved: <path>".
# ============================================================
discussion <- function(cfg) {
  
  alpha  <- 0.05
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  f_in <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_in))
  m <- fread(f_in, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "side") %in% names(m)))
  
  m[, `:=`(pid = as.character(pid), treat = as.character(treat),
           seq = as.character(seq), side = as.character(side))]
  m <- m[side %in% c("H", "T")]                       # bets only
  
  m[, nH        := nchar(gsub("T", "", seq))]
  m[, last      := substr(seq, 6, 6)]
  m[, run       := substr(seq, 5, 5) == last]
  m[, run_sign  := fifelse(run, fifelse(last == "H", 1L, -1L), 0L)]
  m[, imbalance := (2 * nH - 6) / 6]
  m[, abs_imb   := abs(2 * nH - 6)]
  m[, chose_H   := as.integer(side == "H")]
  
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
  
  f <- file.path(path_out, "discussion_split_sample.csv")
  if (!should_skip(f, cfg, "output", "Discussion, raw split-sample")) {
    fwrite(tbl, f); msg("Saved: ", f)
  }
  invisible(tbl)
}