# ============================================================
# 02_descriptive-sequences.R
#
# PURPOSE
#   Descriptive statistics of betting behaviour by sequence and for each
#   treatment overall: betting rate, stake size, side choice and response
#   time.
#
# INPUT
#   path_der/master-sequences.csv
#
# OUTPUT (data/output/csv/)
#   descriptive-sequences-<tr>.csv      one row per sequence, by treatment
#   descriptive-sequences-summary.csv   pooled and by treatment, all trials
#
# NOTES
#   - A bet is a trial with a positive stake. Stakes and the share of Heads
#     are computed over betting trials
#   - Heads and Tails are the side labels of cfg$design$seq$side_labels
#   - Standard deviations are across trials
# ============================================================

descriptive_sequences <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  infile <- file.path(path_der, "master-sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "seq", "treat", "stake", "side", "screen_ms") %in% names(dt)))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, side  := as.character(side)]
  dt[, stake := as.numeric(stake)]
  
  heads_label <- as.character(cfg$design$seq$side_labels$heads)
  tails_label <- as.character(cfg$design$seq$side_labels$tails)
  
  dt[, bet   := !is.na(stake) & stake > 0]
  dt[, bet_H := bet & side == heads_label]
  dt[, bet_T := bet & side == tails_label]
  
  # ---- By sequence ----
  seq_desc <- function(x) {
    out <- x[, .(
      n_pid          = uniqueN(pid),
      bet_rate_mean  = mean(bet),
      stake_mean     = mean(stake[bet]),
      stake_sd       = sd(stake[bet]),
      screen_ms_mean = mean(screen_ms, na.rm = TRUE),
      screen_ms_sd   = sd(screen_ms,   na.rm = TRUE),
      side_H_rate    = mean(side[bet] == heads_label)
    ), by = seq]
    setorder(out, seq)
    out
  }
  
  for (tr in tr_vec) {
    x <- dt[treat == tr]
    if (nrow(x) == 0) next
    f_out <- file.path(path_csv, paste0("descriptive-sequences-", tr, ".csv"))
    fwrite(seq_desc(x), f_out)
    msg("Saved: ", f_out)
  }
  
  # ---- Overall ----
  global_desc <- function(x, group) {
    x_bet <- x[bet == TRUE]
    data.table(
      treatment                    = group,
      n_obs                        = nrow(x),
      n_pid                        = uniqueN(x$pid),
      bet_rate_mean                = mean(x$bet),
      stake_mean_all               = mean(x$stake, na.rm = TRUE),
      stake_sd_all                 = sd(x$stake,   na.rm = TRUE),
      stake_mean_bet               = mean(x_bet$stake),
      stake_sd_bet                 = sd(x_bet$stake),
      side_H_rate_cond_bet         = mean(x_bet$side == heads_label),
      side_T_rate_cond_bet         = mean(x_bet$side == tails_label),
      side_bias_H_minus_T_cond_bet = mean(x_bet$side == heads_label) -
        mean(x_bet$side == tails_label),
      bet_H_rate                   = mean(x$bet_H),
      bet_T_rate                   = mean(x$bet_T),
      screen_ms_mean               = mean(x$screen_ms, na.rm = TRUE),
      screen_ms_sd                 = sd(x$screen_ms,   na.rm = TRUE)
    )
  }
  
  rows <- list(pooled = global_desc(dt, "pooled"))
  for (tr in tr_vec) {
    x <- dt[treat == tr]
    if (nrow(x) > 0) rows[[tr]] <- global_desc(x, tr)
  }
  f_sum <- file.path(path_csv, "descriptive-sequences-summary.csv")
  fwrite(rbindlist(rows, use.names = TRUE, fill = TRUE), f_sum)
  msg("Saved: ", f_sum)
  
  invisible(TRUE)
}
