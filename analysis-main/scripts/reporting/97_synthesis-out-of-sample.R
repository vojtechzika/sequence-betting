# ============================================================
# 97_synthesis-out-of-sample.R
#
# PURPOSE
#   Out-of-sample check of three parameter-free decision rules on the
#   independent pilot sample. The direction of the effects is established by
#   the joint model on the main sample; this script asks whether each rule,
#   applied mechanically to a sample it was not fitted to, forecasts
#   correctly.
#     Betting: predict a bet whenever the sequence is not perfectly balanced
#       and 0.5 on balanced sequences.
#     Side choice against the run: on bets after a terminal run of at least
#       two identical symbols, predict the opposite symbol.
#     Side choice with the imbalance majority: on bets without a terminal
#       run and with an unbalanced sequence, predict the majority symbol.
#   Side choice is evaluated by the accuracy of the rule against the chance
#   level 0.5 (symmetric sides, randomised button position), one mean per
#   participant. Betting is evaluated by the rate predicted from the pilot
#   sequences alone against the realised mean bet rate. Intervals are
#   percentile bootstraps over participants, bounded in [0, 1]. The reported
#   check is whether the benchmark lies inside the interval; it is not a
#   null-hypothesis test.
#
# SAMPLE
#   Pilot sequences of the treatments in cfg$run$treatment.
#
# INPUT
#   raw-data/manual/manual_out_of_sample_sequences.csv
#     columns: pid, seq, side (NB for no bet), treatment
#
# OUTPUT (data/output/csv/)
#   synthesis-out-of-sample-participant-rates.csv
#   synthesis-out-of-sample-summary.csv   rule, treat_label, n_participants,
#     predicted_rate, mean_rate, ci_low, ci_high, ci_method, predicted_in_ci
#
# NOTES
#   - No-cue trials contribute an analytic 0.5 to the predicted bet rate
#   - With very few participants in a treatment the bootstrap interval is
#     fragile; n_participants is reported alongside
#   - Seeded by cfg$run$seed; the number of bootstrap draws is set in
#     settings.R (cfg$synthesis$oos$n_boot)
#   - No model object is read
#   - Output files are overwritten at every run
# ============================================================

oos_seq_features <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  imbalance_signed <- (sum(chars == "H") - sum(chars == "T")) / length(chars)
  last_two <- paste(utils::tail(chars, 2), collapse = "")
  run_dir  <- data.table::fifelse(last_two == "HH", 1L,
                                  data.table::fifelse(last_two == "TT", -1L, 0L))
  list(imbalance_abs = abs(imbalance_signed), imbalance_signed = imbalance_signed,
       run_dir = run_dir)
}

# Across-participant mean and percentile-bootstrap 95% interval, resampling
# participants; bounded in [0, 1] by construction.
oos_rate_ci <- function(x, n_boot) {
  n  <- length(x)
  bs <- replicate(n_boot, mean(x[sample.int(n, n, replace = TRUE)]))
  q  <- stats::quantile(bs, c(.025, .975), names = FALSE)
  data.table::data.table(n_participants = n, mean_rate = mean(x),
                         ci_low = q[1], ci_high = q[2],
                         ci_method = "percentile bootstrap over participants")
}

synthesis_out_of_sample <- function(cfg) {
  
  set.seed(as.integer(cfg$run$seed))
  n_boot <- as.integer(cfg$synthesis$oos$n_boot)
  
  # treatment codes of the pilot file and their labels (FN, FP)
  style        <- treatment_style(cfg)
  TREAT_LABELS <- setNames(style$type, style$treat)
  
  f_pilot <- file.path(path_manual, "manual_out_of_sample_sequences.csv")
  if (!file.exists(f_pilot)) stop("Pilot file not found: ", f_pilot, call. = FALSE)
  dt <- data.table::fread(f_pilot, encoding = "UTF-8")
  stopifnot(all(c("pid", "seq", "side", "treatment") %in% names(dt)))
  data.table::setnames(dt, "treatment", "treat")
  
  dt[, treat := as.character(treat)]
  dt[, pid   := as.character(pid)]
  dt[, seq   := as.character(seq)]
  
  missing_treat <- setdiff(unique(dt$treat), names(TREAT_LABELS))
  if (length(missing_treat) > 0L)
    warning("Dropping pilot rows with treatment code(s) not in the run: ",
            paste(missing_treat, collapse = ", "), call. = FALSE)
  dt <- dt[treat %in% names(TREAT_LABELS)]
  dt[, treat_label := TREAT_LABELS[treat]]
  
  dt <- cbind(dt, data.table::rbindlist(lapply(dt$seq, oos_seq_features)))
  dt[, bet := as.integer(side != "NB")]
  
  # ---- Rules ----
  dt[, predicted_bet := data.table::fifelse(imbalance_abs > 0, 1, 0.5)]
  dt[, predicted_side_run := data.table::fifelse(
    run_dir == 1, "T", data.table::fifelse(run_dir == -1, "H", NA_character_))]
  dt[, majority_side := data.table::fifelse(
    imbalance_signed > 0, "H", data.table::fifelse(imbalance_signed < 0, "T", NA_character_))]
  
  # ---- Betting: aggregate rate ----
  bet_predicted   <- dt[, .(predicted_rate = mean(predicted_bet)), by = treat_label]
  bet_participant <- dt[, .(rule = "bet_vs_imbalance", rate = mean(bet), n_trials = .N),
                        by = .(treat_label, pid)]
  bet_summary <- bet_participant[, oos_rate_ci(rate, n_boot), by = treat_label]
  bet_summary[, rule := "bet_vs_imbalance"]
  bet_summary <- merge(bet_summary, bet_predicted, by = "treat_label")
  
  # ---- Side choice: against the run ----
  side_participant <- dt[bet == 1 & run_dir != 0,
                         .(rule = "side_vs_terminal_run",
                           rate = mean(predicted_side_run == side), n_trials = .N),
                         by = .(treat_label, pid)]
  side_summary <- side_participant[, oos_rate_ci(rate, n_boot), by = treat_label]
  side_summary[, rule := "side_vs_terminal_run"]
  side_summary[, predicted_rate := 0.5]
  
  # ---- Side choice: with the imbalance majority (no-run trials only) ----
  side_imb_participant <- dt[bet == 1 & run_dir == 0 & !is.na(majority_side),
                             .(rule = "side_vs_imbalance_majority",
                               rate = mean(majority_side == side), n_trials = .N),
                             by = .(treat_label, pid)]
  side_imb_summary <- side_imb_participant[, oos_rate_ci(rate, n_boot), by = treat_label]
  side_imb_summary[, rule := "side_vs_imbalance_majority"]
  side_imb_summary[, predicted_rate := 0.5]
  
  # ---- Combine, save ----
  participant_rates <- data.table::rbindlist(
    list(bet_participant, side_participant, side_imb_participant), use.names = TRUE)
  data.table::setorder(participant_rates, rule, treat_label, pid)
  
  f_part <- file.path(path_csv, "synthesis-out-of-sample-participant-rates.csv")
  data.table::fwrite(participant_rates, f_part);  msg("Saved: ", f_part)
  
  summary_tbl <- data.table::rbindlist(list(bet_summary, side_summary, side_imb_summary),
                                       use.names = TRUE)
  summary_tbl[, predicted_in_ci := (predicted_rate >= ci_low) & (predicted_rate <= ci_high)]
  data.table::setcolorder(summary_tbl, c("rule", "treat_label", "n_participants",
                                         "predicted_rate", "mean_rate", "ci_low", "ci_high",
                                         "ci_method", "predicted_in_ci"))
  data.table::setorder(summary_tbl, rule, treat_label)
  
  f_summary <- file.path(path_csv, "synthesis-out-of-sample-summary.csv")
  data.table::fwrite(summary_tbl, f_summary);  msg("Saved: ", f_summary)
  
  invisible(list(participant_rates = participant_rates, summary = summary_tbl))
}