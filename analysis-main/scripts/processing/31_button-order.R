# ============================================================
# 31_button-order.R
#
# PURPOSE
#   Checks whether button order (HT vs TH, i.e. heads first vs tails
#   first) systematically affects choices. This is a randomization
#   check: if randomization worked, neither the betting rate, the Heads
#   rate nor the stake should differ across button-order conditions.
#
# INPUT
#   path_der/master-sequences.csv  -- trial-level analysis table
#
# OUTPUT
#   path_csv/processing_button-order.csv
#
# NOTES
#   - Outcomes per treatment (one row each, with the tested statistic):
#       bet   (mean)         : share of trials with a bet (stake > 0), all trials
#       heads (mean)         : share of Heads among bets, betting trials only
#       stake (mean, median) : stake among bets, betting trials only
#   - Within-participant permutation test: button-order labels are
#     permuted across each participant's trials; the test statistic is
#     the difference in the statistic (mean or median) between the two
#     conditions
#   - Two-sided p-value (1 + #{|perm diff| >= |obs diff|}) / (B + 1)
#   - Heads and no-bet labels are taken from cfg$design$seq$side_labels
#   - Unexpected data (missing side or button order, a treatment
#     without trials, a single button-order level) stops the script
# ============================================================

button_order_check <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  seed   <- as.integer(cfg$run$seed)
  B      <- 2000L
  
  heads_label <- as.character(cfg$design$seq$side_labels$heads)
  nobet_label <- as.character(cfg$design$seq$side_labels$nobet)
  
  infile <- file.path(path_der, "master-sequences.csv")
  stopifnot(file.exists(infile))
  
  dt0 <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "side", "btn_order") %in% names(dt0)))
  
  dt0[, pid       := as.character(pid)]
  dt0[, treat     := as.character(treat)]
  dt0[, stake     := as.numeric(stake)]
  dt0[, side      := as.character(side)]
  dt0[, btn_order := as.character(btn_order)]
  dt0[, y_bet     := as.integer(stake > 0)]
  dt0[, y_heads   := as.integer(side == heads_label)]
  
  # A trial is a bet if and only if a side was chosen
  if (anyNA(dt0$stake) || anyNA(dt0$side) || anyNA(dt0$btn_order) ||
      any((dt0$side == nobet_label) != (dt0$y_bet == 0L))) {
    stop("Missing values, or stake and side are inconsistent.")
  }
  
  f_csv <- file.path(path_csv, "processing_button-order.csv")
  
  # Permutation test of one outcome (column y) on the trials in d
  perm_test <- function(d, y, tr, outcome, stat = mean, stat_name = "mean") {
    
    if (nrow(d) == 0L) stop("No trials for ", outcome, " (tr = ", tr, ").")
    
    levels_btn <- sort(unique(d$btn_order))
    if (length(levels_btn) != 2L) {
      stop("Expected two button-order levels for tr = ", tr, ", found ",
           length(levels_btn), ".")
    }
    
    in1      <- d$btn_order == levels_btn[1]
    obs_diff <- stat(y[in1]) - stat(y[!in1])
    idx_pid  <- split(seq_len(nrow(d)), d$pid)
    
    # Permute button-order labels within participant
    set.seed(seed)
    diff_rep <- numeric(B)
    for (b in seq_len(B)) {
      btn_perm <- d$btn_order
      for (rows in idx_pid) {
        if (length(rows) > 1L) btn_perm[rows] <- sample(btn_perm[rows])
      }
      p1 <- btn_perm == levels_btn[1]
      diff_rep[b] <- stat(y[p1]) - stat(y[!p1])
    }
    
    data.table(
      treatment   = tr,
      outcome     = outcome,
      statistic   = stat_name,
      btn_order_1 = levels_btn[1],
      btn_order_2 = levels_btn[2],
      n_1         = sum(in1),
      n_2         = sum(!in1),
      value_1     = stat(y[in1]),
      value_2     = stat(y[!in1]),
      diff        = obs_diff,
      n_perm      = B,
      p_two_sided = (1 + sum(abs(diff_rep) >= abs(obs_diff) - 1e-12)) / (B + 1)
    )
  }
  
  all_rows <- list()
  
  for (tr in tr_vec) {
    
    d_all <- dt0[treat == tr]
    d_bet <- d_all[y_bet == 1L]
    
    all_rows[[paste(tr, "bet")]]   <- perm_test(d_all, d_all$y_bet,   tr, "bet")
    all_rows[[paste(tr, "heads")]] <- perm_test(d_bet, d_bet$y_heads, tr, "heads")
    all_rows[[paste(tr, "stake")]] <- perm_test(d_bet, d_bet$stake,   tr, "stake")
    all_rows[[paste(tr, "stake median")]] <-
      perm_test(d_bet, d_bet$stake, tr, "stake", stat = median, stat_name = "median")
  }
  
  fwrite(rbindlist(all_rows), f_csv)
  msg("Saved: ", f_csv)
  
  invisible(TRUE)
}
