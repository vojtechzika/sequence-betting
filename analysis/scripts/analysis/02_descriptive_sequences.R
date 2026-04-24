# ============================================================
# 02_descriptive_sequences.R
#
# PURPOSE
#   Computes sequence-level and global descriptive statistics
#   for betting behavior, stake size, side choice, and response
#   time. Produces outputs pooled and by treatment.
#
# INPUT
#   path_src/master_sequences.csv
#
# OUTPUT
#   path_out/descriptive_sequences_pooled.csv  -- sequence-level
#   path_out/descriptive_sequences_<tr>.csv    -- sequence-level by treatment
#   path_out/descriptive_sequences.csv         -- global summary by treatment
# ============================================================

descriptive_sequences <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "seq", "treat", "block", "btn_order",
                "sex", "stake", "side", "screen_ms")
  stopifnot(all(required %in% names(dt)))
  
  dt[, pid       := as.character(pid)]
  dt[, treat     := as.character(treat)]
  dt[, seq       := as.character(seq)]
  dt[, side      := as.character(side)]
  dt[, btn_order := as.character(btn_order)]
  dt[, sex       := fifelse(sex %in% c("F", "M"), sex, "O")]
  dt[is.na(sex) | !nzchar(sex), sex := "O"]
  
  dt[, bet   := !is.na(stake) & stake > 0]
  dt[, bet_H := bet & side == "H"]
  dt[, bet_T := bet & side == "T"]
  
  heads_label <- as.character(cfg$design$seq$side_labels$heads)
  
  # ---- Sequence-level descriptives ----
  seq_desc <- function(x) {
    out <- x[, .(
      n_pid          = uniqueN(pid),
      bet_rate_mean  = mean(bet,                      na.rm = TRUE),
      bet_rate_sd    = sd(bet,                         na.rm = TRUE),
      stake_mean     = mean(stake[bet],                na.rm = TRUE),
      stake_sd       = sd(stake[bet],                  na.rm = TRUE),
      screen_ms_mean = mean(screen_ms,                 na.rm = TRUE),
      screen_ms_sd   = sd(screen_ms,                   na.rm = TRUE),
      side_H_rate    = mean(side[bet] == heads_label,  na.rm = TRUE)
    ), by = seq]
    setorder(out, seq)
    out
  }
  
  # ---- Global descriptives ----
  global_desc <- function(x, tag) {
    x_bet <- x[bet == TRUE]
    data.table(
      treatment                    = tag,
      n_obs                        = nrow(x),
      n_pid                        = uniqueN(x$pid),
      bet_rate_mean                = mean(x$bet,          na.rm = TRUE),
      bet_rate_sd                  = sd(as.numeric(x$bet), na.rm = TRUE),
      stake_mean_all               = mean(x$stake,         na.rm = TRUE),
      stake_sd_all                 = sd(x$stake,           na.rm = TRUE),
      stake_mean_bet               = mean(x_bet$stake,     na.rm = TRUE),
      stake_sd_bet                 = sd(x_bet$stake,       na.rm = TRUE),
      side_H_rate_cond_bet         = mean(x_bet$side == "H", na.rm = TRUE),
      side_T_rate_cond_bet         = mean(x_bet$side == "T", na.rm = TRUE),
      side_bias_H_minus_T_cond_bet = mean(x_bet$side == "H", na.rm = TRUE) -
        mean(x_bet$side == "T", na.rm = TRUE),
      bet_H_rate                   = mean(x$bet_H,         na.rm = TRUE),
      bet_T_rate                   = mean(x$bet_T,         na.rm = TRUE),
      screen_ms_mean               = mean(x$screen_ms,     na.rm = TRUE),
      screen_ms_sd                 = sd(x$screen_ms,       na.rm = TRUE)
    )
  }
  
  # ---- Run sequence-level ----
  run_seq <- function(x, tag) {
    f_out <- file.path(path_out, paste0("descriptive_sequences_", tag, ".csv"))
    if (!should_skip(f_out, cfg, "output",
                     paste0("Descriptive sequences (", tag, ")"))) {
      out <- seq_desc(x)
      fwrite(out, f_out)
      msg("Saved: ", f_out)
    }
  }
  
  # by treatment only
  for (tr in tr_vec) {
    x <- dt[treat == tr]
    if (nrow(x) == 0) next
    run_seq(x, tr)
  }
  
  # ---- summary ----
  f_global <- file.path(path_out, "descriptive_sequences_summary.csv")
  if (!should_skip(f_global, cfg, "output", "Descriptive sequences (global)")) {
    global_rows        <- list()
    global_rows[["pooled"]] <- global_desc(dt, "pooled")
    for (tr in tr_vec) {
      x <- dt[treat == tr]
      if (nrow(x) == 0) next
      global_rows[[tr]] <- global_desc(x, tr)
    }
    out_global <- rbindlist(global_rows, use.names = TRUE, fill = TRUE)
    fwrite(out_global, f_global)
    msg("Saved: ", f_global)
  }
  
  invisible(TRUE)
}