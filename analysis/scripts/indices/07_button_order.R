# ============================================================
# 07_button_order.R
#
# PURPOSE
#   Checks whether button order (HT vs TH) systematically affects
#   the Heads choice rate. This is a randomization check — if
#   randomization worked correctly, Heads rates should be identical
#   across button order conditions.
#
# INPUT
#   path_src/master_sequences.csv
#
# OUTPUT
#   path_out/button_order_check.csv
#
# NOTES
#   - Restricted to betting trials only
#   - Within-participant permutation test on Heads choice rate
#   - A significant result indicates positional response bias
# ============================================================

button_order_check <- function(cfg) {
  
  tr_vec      <- unique(as.character(cfg$run$treatment))
  seed        <- as.integer(cfg$run$seed)
  B           <- 2000L
  heads_label <- as.character(cfg$design$seq$side_labels$heads)
  
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt0 <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "stake", "side", "btn_order") %in% names(dt0)))
  
  dt0[, pid       := as.character(pid)]
  dt0[, treat     := as.character(treat)]
  dt0[, stake     := as.numeric(stake)]
  dt0[, y_bet     := as.integer(stake > 0)]
  dt0[, y_heads   := as.integer(side == heads_label)]
  dt0[, btn_order := as.character(btn_order)]
  
  f_csv <- file.path(path_out, "button_order_check.csv")
  
  if (should_skip(
    paths = f_csv,
    cfg   = cfg,
    type  = "output",
    label = "Button order check"
  )) return(invisible(NULL))
  
  all_rows <- list()
  
  for (tr in tr_vec) {
    
    d <- dt0[treat == tr & y_bet == 1L & !is.na(side) & !is.na(btn_order)]
    
    if (nrow(d) == 0) {
      warning("No betting trials for button order check (tr=", tr, ")")
      next
    }
    
    obs        <- d[, .(n = .N, n_heads = sum(y_heads), heads_rate = mean(y_heads)),
                    by = btn_order]
    rates      <- setNames(obs$heads_rate, obs$btn_order)
    btn_orders <- sort(unique(d$btn_order))
    
    if (length(btn_orders) < 2) {
      warning("Only one button order level found for tr=", tr)
      next
    }
    
    obs_diff   <- as.numeric(rates[btn_orders[1]] - rates[btn_orders[2]])
    idx_by_pid <- split(seq_len(nrow(d)), d$pid)
    
    set.seed(seed)
    diff_rep <- numeric(B)
    for (b in seq_len(B)) {
      btn_perm <- d$btn_order
      for (rows in idx_by_pid) {
        if (length(rows) <= 1L) next
        btn_perm[rows] <- sample(btn_perm[rows], size = length(rows), replace = FALSE)
      }
      diff_rep[b] <- mean(d$y_heads[btn_perm == btn_orders[1]]) -
        mean(d$y_heads[btn_perm == btn_orders[2]])
    }
    
    p_two_sided <- mean(abs(diff_rep) >= abs(obs_diff))
    
    msg("Button order check (", tr, "):",
        " ", btn_orders[1], "=", sprintf("%.3f", rates[btn_orders[1]]),
        " vs ", btn_orders[2], "=", sprintf("%.3f", rates[btn_orders[2]]),
        " | diff=", sprintf("%.3f", obs_diff),
        " | p=", sprintf("%.4f", p_two_sided),
        if (p_two_sided < 0.05) " *** BIAS DETECTED" else " (no bias)")
    
    all_rows[[tr]] <- data.table(
      treatment    = tr,
      btn_order_1  = btn_orders[1],
      btn_order_2  = btn_orders[2],
      n_1          = obs[btn_order == btn_orders[1], n],
      n_2          = obs[btn_order == btn_orders[2], n],
      heads_rate_1 = obs[btn_order == btn_orders[1], heads_rate],
      heads_rate_2 = obs[btn_order == btn_orders[2], heads_rate],
      diff         = obs_diff,
      p_two_sided  = p_two_sided,
      bias_detected = p_two_sided < 0.05
    )
  }
  
  fwrite(rbindlist(all_rows), f_csv)
  msg("Saved: ", f_csv)
  
  invisible(TRUE)
}