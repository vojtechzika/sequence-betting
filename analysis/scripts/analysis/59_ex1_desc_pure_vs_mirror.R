# ============================================================
# 59_ex1_desc_pure_vs_mirror.R
#
# Cross-tabulates participant behavior on pure sequences
# (HHHHHH, TTTTTT) against four mirror-pair sequences to
# test whether pure-sequence hot-handers are gamblerish on
# sequences with strong terminal runs.
#
# Mirror pairs (H<->T flips of each other):
#   HHTHTT <-> TTHTHH   (end TT / end HH)
#   HTHHTT <-> THTTHH   (end TT / end HH)
#   TTHHTT <-> HHTTHH   (end TT / end HH)
#   TTTHHH <-> HHHTTT   (end HH / end TT)
#
# Extended: reports raw p(Heads) on HHHHHH and TTTTTT
# broken down by pure_type, to assess whether hot-hand
# consistency reflects genuine streak-following or an
# unconditional Heads preference.
#
# OUTPUTS
#   path_out/ex1_pure_vs_mirror.csv
#   path_out/ex1_pure_vs_mirror_xtab.csv
#   path_out/ex1_pure_vs_mirror_ph.csv     <- new
# ============================================================

ex1_pure_vs_mirror <- function(cfg) {
  
  lab_heads <- as.character(cfg$design$seq$side_labels$heads)
  lab_tails <- as.character(cfg$design$seq$side_labels$tails)
  pure_H    <- as.character(cfg$design$seq$anchor_labels$pure_heads)
  pure_T    <- as.character(cfg$design$seq$anchor_labels$pure_tails)
  tr_vec    <- unique(as.character(cfg$run$treatment))
  
  mirror_cfg <- data.table(
    seq       = c("HHTHTT", "TTHTHH",
                  "HTHHTT", "THTTHH",
                  "TTHHTT", "HHTTHH",
                  "HHHTTT", "TTTHHH"),
    gambler_h = c(1L, 0L,
                  1L, 0L,
                  1L, 0L,
                  1L, 0L)
  )
  
  dt <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  dt <- dt[
    treat %in% tr_vec &
      is.finite(stake) & stake > 0 &
      side %in% c(lab_heads, lab_tails)
  ]
  dt[, h := as.integer(side == lab_heads)]
  
  # ---- Pure sequence type ----
  d_pure <- dt[seq %in% c(pure_H, pure_T),
               .(h = as.integer(mean(h) >= 0.5)), by = .(treat, pid, seq)]
  
  d_pH <- d_pure[seq == pure_H, .(treat, pid, h_pureH = h)]
  d_pT <- d_pure[seq == pure_T, .(treat, pid, h_pureT = h)]
  
  pure_tbl <- merge(d_pH, d_pT, by = c("treat", "pid"), all = TRUE)
  pure_tbl[, pure_type := fcase(
    h_pureH == 1L & h_pureT == 0L,  "hothand",
    h_pureH == 0L & h_pureT == 1L,  "gambler",
    h_pureH == 1L & h_pureT == 1L,  "inconsistent",
    h_pureH == 0L & h_pureT == 0L,  "inconsistent",
    default = "missing"
  )]
  
  # ---- Raw p(Heads) by pure_type ----
  # Merge pure_type back onto raw trial data for pure sequences
  ph_raw <- dt[seq %in% c(pure_H, pure_T)]
  ph_raw <- merge(ph_raw,
                  pure_tbl[, .(treat, pid, pure_type)],
                  by = c("treat", "pid"), all.x = TRUE)
  
  ph_summary <- ph_raw[
    !is.na(pure_type),
    .(
      n_participants = uniqueN(pid),
      n_trials       = .N,
      p_heads_mean   = mean(h),
      p_heads_sd     = sd(h),
      p_heads_se     = sd(h) / sqrt(.N),
      ci_lo          = mean(h) - qt(0.975, df = .N - 1) * sd(h) / sqrt(.N),
      ci_hi          = mean(h) + qt(0.975, df = .N - 1) * sd(h) / sqrt(.N)
    ),
    by = .(treat, seq, pure_type)
  ]
  
  ph_summary[, pure_type := factor(pure_type,
                                   levels = c("hothand", "gambler",
                                              "inconsistent", "missing"))]
  setorder(ph_summary, treat, seq, pure_type)
  
  msg("\n=== Raw p(Heads) on pure sequences by pure_type ===")
  for (tr in tr_vec) {
    msg("\nTreatment: ", tr)
    print(ph_summary[treat == tr,
                     .(seq, pure_type, n_participants, n_trials,
                       p_heads_mean = round(p_heads_mean, 3),
                       ci_lo        = round(ci_lo,        3),
                       ci_hi        = round(ci_hi,        3))],
          nrows = Inf)
  }
  
  # ---- Mirror sequence responses ----
  d_mirror <- dt[seq %in% mirror_cfg$seq]
  d_mirror <- merge(d_mirror, mirror_cfg, by = "seq")
  
  d_mirror[, gambler_score := fifelse(h == gambler_h, 1L, -1L)]
  
  mirror_tbl <- d_mirror[, .(
    n_mirror           = .N,
    n_seqs_observed    = uniqueN(seq),
    gambler_score_mean = mean(gambler_score),
    mirror_type        = fcase(
      mean(gambler_score) >  0, "gambler",
      mean(gambler_score) <  0, "hothand",
      mean(gambler_score) == 0, "split"
    )
  ), by = .(treat, pid)]
  
  # ---- Merge ----
  combined <- merge(pure_tbl, mirror_tbl, by = c("treat", "pid"), all = TRUE)
  
  # ---- Cross-tabulation ----
  xtab <- combined[
    pure_type %in% c("hothand", "gambler", "inconsistent") &
      !is.na(mirror_type),
    .(n = .N),
    by = .(treat, pure_type, mirror_type)
  ]
  
  xtab[, pct_within_pure := round(100 * n / sum(n), 1), by = .(treat, pure_type)]
  xtab[, pure_type   := factor(pure_type,   levels = c("hothand", "gambler", "inconsistent"))]
  xtab[, mirror_type := factor(mirror_type, levels = c("gambler", "hothand", "split"))]
  setorder(xtab, treat, pure_type, mirror_type)
  
  msg("\n=== Pure vs Mirror cross-tabulation ===")
  msg("(mirror_type: gambler = consistent with terminal-run fallacy)")
  for (tr in tr_vec) {
    msg("\nTreatment: ", tr)
    print(xtab[treat == tr])
  }
  
  # ---- Save ----
  f_pid  <- file.path(path_out, "ex1_pure_vs_mirror.csv")
  f_xtab <- file.path(path_out, "ex1_pure_vs_mirror_xtab.csv")
  f_ph   <- file.path(path_out, "ex1_pure_vs_mirror_ph.csv")
  
  fwrite(combined,   f_pid)
  fwrite(xtab,       f_xtab)
  fwrite(ph_summary, f_ph)
  
  msg("\nSaved: ", f_pid)
  msg("Saved: ", f_xtab)
  msg("Saved: ", f_ph)
  
  invisible(list(participants = combined, xtab = xtab, ph_summary = ph_summary))
}