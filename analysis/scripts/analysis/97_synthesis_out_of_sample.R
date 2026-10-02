# ============================================================
# 97_synthesis_out_of_sample.R
#
# PURPOSE
#   Produces the "final, lightweight check" reported in the Synthesis
#   section: whether minimal, PARAMETER-FREE decision rules (no fitted
#   coefficients anywhere) forecast correctly in the independent pilot
#   sample. The direction of the effects is established by the fitted
#   slopes on the main sample; this script only asks whether each rule,
#   applied mechanically to a sample it was never fit to, gets the right
#   answer -- including the prediction that one of them should NOT.
#
# RULES
#   Betting (extensive margin): predict "bet" (1) whenever the sequence is
#   not perfectly balanced; on perfectly balanced (3H-3T) sequences the
#   rule has no signal and predicts the coin-flip expectation (0.5).
#   Side choice, against the run (conditional on betting): against the
#   terminal run, on bettor trials with a run of >= 2 identical symbols.
#   Expected to beat chance (the run effect is homogeneous across
#   participants in the fitted model).
#   Side choice, with the imbalance majority (conditional on betting):
#   predict the majority symbol, on bettor trials that are NOT balanced and
#   have NO terminal run. The fitted model implies only a weak,
#   heterogeneous population pull toward the majority, so this rule is
#   expected to land close to chance.
#
# METHOD
#   Side choice: per-trial accuracy vs. the fixed chance level 0.5 (H/T
#     symmetric by design; button position randomized). One mean per
#     participant on rule-eligible trials.
#   Betting: aggregate rate -- predicted_rate (mechanical, from the pilot
#     sequences alone) vs. realized_rate (mean of participants' own bet
#     rates).
#   Uncertainty: PERCENTILE BOOTSTRAP OVER PARTICIPANTS (the independent
#     unit), n_boot resamples. Unlike the t-interval used previously, the
#     bootstrap interval is bounded in [0, 1] by construction (a t-interval
#     produced an upper bound of 105.7% in the FP pilot). With very small
#     samples (FP: n = 5) any interval is fragile and the bootstrap tends
#     to under-cover; n is reported alongside and the text should treat FP
#     as uninformative rather than as a test.
#   The reported check is whether the benchmark (0.5 for side choice; the
#   mechanical predicted_rate for betting) lies inside the interval
#   (predicted_in_ci). Deliberately not a null-hypothesis test.
#
# INPUT
#   path_src/manual_out_of_sample_sequences.csv
#     columns: pid, seq, side, treatment (side == "NB" for no bet)
#
# OUTPUT
#   path_out/synthesis_out_of_sample_participant_rates.csv
#   path_out/synthesis_out_of_sample_summary.csv
#     columns: rule, treat_label, n_participants, predicted_rate,
#     mean_rate, ci_low, ci_high, ci_method, predicted_in_ci
#
# NOTES
#   - "No cue" trials contribute an analytic 0.5 to predicted_bet (no
#     simulated coin flips).
#   - Seeded from cfg$run$seed; n_boot from cfg$synthesis$n_boot_oos
#     (default 10000).
#   - Fits nothing and reads no .rds model object.
# ============================================================

TREAT_LABELS <- c(m25 = "FN", m19 = "FP")

seq_features <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  imbalance_signed <- (sum(chars == "H") - sum(chars == "T")) / length(chars)
  last_two <- paste(utils::tail(chars, 2), collapse = "")
  run_dir  <- data.table::fifelse(last_two == "HH", 1L,
                                  data.table::fifelse(last_two == "TT", -1L, 0L))
  list(imbalance_abs = abs(imbalance_signed), imbalance_signed = imbalance_signed,
       run_dir = run_dir)
}

# ---- Helper: across-participant mean + percentile-bootstrap 95% CI,
# resampling participants. Bounded in [0, 1] by construction. ----
rate_ci <- function(x, n_boot) {
  n  <- length(x)
  bs <- replicate(n_boot, mean(x[sample.int(n, n, replace = TRUE)]))
  q  <- stats::quantile(bs, c(.025, .975), names = FALSE)
  data.table::data.table(n_participants = n, mean_rate = mean(x),
                         ci_low = q[1], ci_high = q[2],
                         ci_method = "percentile bootstrap over participants")
}

synthesis_out_of_sample <- function(cfg) {
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  set.seed(as.integer(cfg$run$seed %||% 12345L))
  n_boot <- as.integer(cfg$synthesis$n_boot_oos %||% 10000L)
  
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  
  f_pilot <- file.path(path_src, "manual_out_of_sample_sequences.csv")
  stopifnot(file.exists(f_pilot))
  dt <- data.table::fread(f_pilot, encoding = "UTF-8")
  stopifnot(all(c("pid", "seq", "side", "treatment") %in% names(dt)))
  data.table::setnames(dt, "treatment", "treat")
  
  dt[, treat := as.character(treat)]
  dt[, pid   := as.character(pid)]
  dt[, seq   := as.character(seq)]
  
  missing_treat <- setdiff(unique(dt$treat), names(TREAT_LABELS))
  if (length(missing_treat) > 0L) {
    msg("synthesis_out_of_sample: dropping row(s) with unrecognized treatment code(s): ",
        paste(missing_treat, collapse = ", "))
  }
  dt <- dt[treat %in% names(TREAT_LABELS)]
  dt[, treat_label := TREAT_LABELS[treat]]
  
  dt <- cbind(dt, data.table::rbindlist(lapply(dt$seq, seq_features)))
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
  bet_summary <- bet_participant[, rate_ci(rate, n_boot), by = treat_label]
  bet_summary[, rule := "bet_vs_imbalance"]
  bet_summary <- merge(bet_summary, bet_predicted, by = "treat_label")
  
  # ---- Side choice: against the run ----
  side_participant <- dt[bet == 1 & run_dir != 0,
                         .(rule = "side_vs_terminal_run",
                           rate = mean(predicted_side_run == side), n_trials = .N),
                         by = .(treat_label, pid)]
  side_summary <- side_participant[, rate_ci(rate, n_boot), by = treat_label]
  side_summary[, rule := "side_vs_terminal_run"]
  side_summary[, predicted_rate := 0.5]
  
  # ---- Side choice: with the imbalance majority (no-run trials only) ----
  side_imb_participant <- dt[bet == 1 & run_dir == 0 & !is.na(majority_side),
                             .(rule = "side_vs_imbalance_majority",
                               rate = mean(majority_side == side), n_trials = .N),
                             by = .(treat_label, pid)]
  side_imb_summary <- side_imb_participant[, rate_ci(rate, n_boot), by = treat_label]
  side_imb_summary[, rule := "side_vs_imbalance_majority"]
  side_imb_summary[, predicted_rate := 0.5]
  
  # ---- Combine, save ----
  participant_rates <- data.table::rbindlist(
    list(bet_participant, side_participant, side_imb_participant), use.names = TRUE)
  data.table::setorder(participant_rates, rule, treat_label, pid)
  
  f_part <- file.path(path_out, "synthesis_out_of_sample_participant_rates.csv")
  if (!should_skip(f_part, cfg, "output", "Out-of-sample participant-level rates")) {
    data.table::fwrite(participant_rates, f_part); msg("Saved: ", f_part)
  }
  
  summary_tbl <- data.table::rbindlist(list(bet_summary, side_summary, side_imb_summary),
                                       use.names = TRUE)
  summary_tbl[, predicted_in_ci := (predicted_rate >= ci_low) & (predicted_rate <= ci_high)]
  data.table::setcolorder(summary_tbl, c("rule", "treat_label", "n_participants",
                                         "predicted_rate", "mean_rate", "ci_low", "ci_high",
                                         "ci_method", "predicted_in_ci"))
  data.table::setorder(summary_tbl, rule, treat_label)
  
  f_summary <- file.path(path_out, "synthesis_out_of_sample_summary.csv")
  if (!should_skip(f_summary, cfg, "output", "Out-of-sample predicted-vs-realized rate summary")) {
    data.table::fwrite(summary_tbl, f_summary); msg("Saved: ", f_summary)
  }
  
  invisible(list(participant_rates = participant_rates, summary = summary_tbl))
}