# ============================================================
# 01_descriptive-participants.R
#
# PURPOSE
#   Participant-level descriptive statistics: demographics, risk
#   preferences, optimism and behavioural aggregates from the sequence
#   task. Produces one file with a row per participant and a summary table
#   by group.
#
# INPUT
#   path_src/participants.csv           demographics, treatment, payoff
#   path_der/master-sequences.csv       trials of the sequence task
#   path_der/lotr-scored.csv            LOT-R scores
#   path_der/switching-point-r.csv      classical (switching-point) r
#   model "mpl": scored-<tr>.csv        posterior mean of r and the Holt-Laury
#                                       consistency flag, per treatment fit
#
# OUTPUT
#   data/output/csv/descriptive-participants.csv        one row per participant
#   data/output/csv/descriptive-participants-summary.csv   group summaries
#
# GROUPS (rows of the summary)
#   pooled | pooled/Women | pooled/Men
#   pooled/consistent | pooled/consistent/Women | pooled/consistent/Men
#   <tr> | <tr>/Women | <tr>/Men
#   <tr>/consistent | <tr>/consistent/Women | <tr>/consistent/Men
#   "consistent" are the participants with a consistent Holt-Laury choice
#   pattern. Participants whose sex is neither F nor M appear in the
#   groups without a sex split.
#
# NOTES
#   - Behavioural variables are participant means. In the participant file
#     the *_sd_within columns are the standard deviations across a
#     participant's own trials
#   - In the summary, <variable>_mean is the mean and <variable>_sd the
#     standard deviation across participants of the participant-level means
#   - Stakes are those of betting trials; the share of Heads is computed
#     over betting trials
# ============================================================

descriptive_participants <- function(cfg) {
  
  tr_vec      <- unique(as.character(cfg$run$treatment))
  heads_label <- as.character(cfg$design$seq$side_labels$heads)
  
  f_par     <- file.path(path_src, "participants.csv")
  f_master  <- file.path(path_der, "master-sequences.csv")
  f_lotr    <- file.path(path_der, "lotr-scored.csv")
  f_classic <- file.path(path_der, "switching-point-r.csv")
  stopifnot(file.exists(f_par), file.exists(f_master), file.exists(f_lotr),
            file.exists(f_classic))
  
  # ---- Participant-level data ----
  par <- fread(f_par, encoding = "UTF-8")
  stopifnot(all(c("pid", "age", "sex", "treat", "payoff") %in% names(par)))
  par[, pid   := as.character(pid)]
  par[, treat := as.character(treat)]
  par[, sex   := fifelse(sex %in% c("F", "M"), sex, "O")]
  par[is.na(sex) | !nzchar(sex), sex := "O"]
  
  # posterior mean of r and consistency flag, from the fit of the participant's treatment
  mpl <- rbindlist(lapply(tr_vec, function(tr)
    fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")[
      , .(pid = as.character(pid), r_mean, inconsistent = as.integer(inconsistent))]))
  stopifnot(uniqueN(mpl$pid) == nrow(mpl))
  
  lotr <- fread(f_lotr, encoding = "UTF-8")
  lotr[, pid := as.character(pid)]
  stopifnot(uniqueN(lotr$pid) == nrow(lotr))
  
  classic <- fread(f_classic, encoding = "UTF-8")
  classic[, pid := as.character(pid)]
  
  dt <- merge(par[, .(pid, age, sex, treat, payoff)], mpl, by = "pid", all.x = TRUE)
  dt <- merge(dt, lotr[, .(pid, lotr_score)], by = "pid", all.x = TRUE)
  dt <- merge(dt, classic[, .(pid, r_hl)],    by = "pid", all.x = TRUE)
  
  # ---- Behavioural aggregates from the trials ----
  master <- fread(f_master, encoding = "UTF-8")
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[, side  := as.character(side)]
  master[, bet   := !is.na(stake) & stake > 0]
  
  behav <- master[, .(
    n_trials           = .N,
    bet_rate_mean      = mean(bet),
    bet_rate_sd_within = sd(as.numeric(bet)),
    stake_mean         = mean(stake[bet]),
    stake_sd_within    = sd(stake[bet]),
    screen_ms_mean      = mean(screen_ms, na.rm = TRUE),
    screen_ms_sd_within = sd(screen_ms, na.rm = TRUE),
    side_H_rate        = mean(side[bet] == heads_label)
  ), by = .(pid, treat)]
  
  dt <- merge(dt, behav, by = c("pid", "treat"), all.x = TRUE)
  
  f_pid <- file.path(path_csv, "descriptive-participants.csv")
  fwrite(dt, f_pid)
  msg("Saved: ", f_pid)
  
  # ---- Group summaries ----
  one_row <- function(x, group) {
    if (nrow(x) == 0) return(NULL)
    x[, .(
      group                = group,
      N                    = uniqueN(pid),
      age_mean             = mean(age,        na.rm = TRUE),
      age_sd               = sd(age,          na.rm = TRUE),
      sex_F_n              = sum(sex == "F"),
      sex_M_n              = sum(sex == "M"),
      sex_O_n              = sum(sex == "O"),
      payoff_mean          = mean(payoff,     na.rm = TRUE),
      payoff_sd            = sd(payoff,       na.rm = TRUE),
      lotr_mean            = mean(lotr_score, na.rm = TRUE),
      lotr_sd              = sd(lotr_score,   na.rm = TRUE),
      r_bayes_mean         = mean(r_mean,     na.rm = TRUE),
      r_bayes_sd           = sd(r_mean,       na.rm = TRUE),
      r_classic_mean       = mean(r_hl,       na.rm = TRUE),
      r_classic_sd         = sd(r_hl,         na.rm = TRUE),
      n_r_bayes            = sum(!is.na(r_mean)),
      n_r_classic          = sum(!is.na(r_hl)),
      hl_inconsistent_n    = sum(inconsistent == 1L,  na.rm = TRUE),
      hl_inconsistent_frac = mean(inconsistent == 1L, na.rm = TRUE),
      bet_rate_mean        = mean(bet_rate_mean,  na.rm = TRUE),
      bet_rate_sd_between  = sd(bet_rate_mean,    na.rm = TRUE),
      stake_mean_bet       = mean(stake_mean,     na.rm = TRUE),   # mean over participants of their mean stake when betting
      stake_sd_bet_between = sd(stake_mean,       na.rm = TRUE),
      screen_ms_mean       = mean(screen_ms_mean, na.rm = TRUE),
      screen_ms_sd         = sd(screen_ms_mean,   na.rm = TRUE),
      side_H_rate_mean     = mean(side_H_rate,    na.rm = TRUE),
      side_H_rate_sd       = sd(side_H_rate,      na.rm = TRUE)
    )]
  }
  
  rows <- list()
  add_group <- function(d, prefix) {
    rows[[prefix]]                   <<- one_row(d,             prefix)
    rows[[paste0(prefix, "/Women")]] <<- one_row(d[sex == "F"], paste0(prefix, "/Women"))
    rows[[paste0(prefix, "/Men")]]   <<- one_row(d[sex == "M"], paste0(prefix, "/Men"))
    d_cons <- d[inconsistent == 0L]
    rows[[paste0(prefix, "/consistent")]]       <<- one_row(d_cons,             paste0(prefix, "/consistent"))
    rows[[paste0(prefix, "/consistent/Women")]] <<- one_row(d_cons[sex == "F"], paste0(prefix, "/consistent/Women"))
    rows[[paste0(prefix, "/consistent/Men")]]   <<- one_row(d_cons[sex == "M"], paste0(prefix, "/consistent/Men"))
  }
  
  add_group(dt, "pooled")
  for (tr in tr_vec) add_group(dt[treat == tr], tr)
  
  f_sum <- file.path(path_csv, "descriptive-participants-summary.csv")
  fwrite(rbindlist(Filter(Negate(is.null), rows), use.names = TRUE, fill = TRUE), f_sum)
  msg("Saved: ", f_sum)
  
  invisible(TRUE)
}
