# ============================================================
# 01_descriptive_participants.R
#
# PURPOSE
#   Computes participant-level descriptive statistics including
#   demographics, risk preferences, optimism scores, and
#   behavioral aggregates from the sequence task.
#   Produces both a per-participant file and group-level summaries.
#
# INPUT
#   path_src/participants.csv
#   path_src/master_sequences.csv
#   path_out/lotr_scored.csv
#   path_out/mpl_scored_pooled.csv
#   path_out/mpl_switching.csv        -- classical r estimates
#
# OUTPUT
#   path_out/descriptive_participants.csv         -- one row per participant
#   path_out/descriptive_participants_summary.csv -- group-level summaries
#
# ROW STRUCTURE (summary)
#   pooled | pooled/Women | pooled/Men
#   pooled/consistent | pooled/consistent/Women | pooled/consistent/Men
#   <tr> | <tr>/Women | <tr>/Men
#   <tr>/consistent | <tr>/consistent/Women | <tr>/consistent/Men
#   (consistent rows only if cfg$run$consistent_only == TRUE)
# ============================================================

descriptive_participants <- function(cfg) {
  
  tr_vec          <- unique(as.character(cfg$run$treatment))
  consistent_only <- isTRUE(cfg$run$consistent_only)
  heads_label     <- as.character(cfg$design$seq$side_labels$heads)
  
  f_par     <- file.path(path_src, "participants.csv")
  f_master  <- file.path(path_src, "master_sequences.csv")
  f_lotr    <- file.path(path_out, "lotr_scored.csv")
  f_mpl     <- file.path(path_out, "mpl_scored_pooled.csv")
  f_classic <- file.path(path_out, "mpl_switching.csv")
  
  stopifnot(file.exists(f_par), file.exists(f_master),
            file.exists(f_mpl), file.exists(f_lotr))
  
  f_out_pid <- file.path(path_out, "descriptive_participants.csv")
  f_out_sum <- file.path(path_out, "descriptive_participants_summary.csv")
  
  skip_pid <- should_skip(f_out_pid, cfg, "output", "Descriptive participants (individual)")
  skip_sum <- should_skip(f_out_sum, cfg, "output", "Descriptive participants (summary)")
  
  if (skip_pid && skip_sum) return(invisible(NULL))
  
  # ---- Load and merge participant-level data ----
  par  <- fread(f_par,    encoding = "UTF-8")
  mpl  <- fread(f_mpl,    encoding = "UTF-8")
  lotr <- fread(f_lotr,   encoding = "UTF-8")
  
  stopifnot(all(c("pid", "age", "sex", "treat", "payoff") %in% names(par)))
  stopifnot(all(c("pid", "r_mean", "inconsistent") %in% names(mpl)))
  stopifnot(all(c("pid", "lotr_score") %in% names(lotr)))
  
  par[,  pid          := as.character(pid)]
  par[,  sex          := fifelse(sex %in% c("F", "M"), sex, "O")]
  par[is.na(sex) | !nzchar(sex), sex := "O"]
  par[,  treat        := as.character(treat)]
  mpl[,  pid          := as.character(pid)]
  mpl[,  inconsistent := as.integer(inconsistent)]
  lotr[, pid          := as.character(pid)]
  
  stopifnot(mpl[,  uniqueN(pid)] == nrow(mpl))
  stopifnot(lotr[, uniqueN(pid)] == nrow(lotr))
  
  dt <- merge(par[, .(pid, age, sex, treat, payoff)],
              mpl[,  .(pid, r_mean, inconsistent)], by = "pid", all.x = TRUE)
  dt <- merge(dt, lotr[, .(pid, lotr_score)],       by = "pid", all.x = TRUE)
  
  if (file.exists(f_classic)) {
    classic <- fread(f_classic, encoding = "UTF-8")
    classic[, pid := as.character(pid)]
    dt <- merge(dt, classic[, .(pid, r_hl)], by = "pid", all.x = TRUE)
  } else {
    dt[, r_hl := NA_real_]
  }
  
  # ---- Compute behavioral aggregates from master_sequences ----
  master <- fread(f_master, encoding = "UTF-8")
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[, side  := as.character(side)]
  master[, bet   := !is.na(stake) & stake > 0]
  
  behav <- master[, .(
    n_trials       = .N,
    bet_rate_mean  = mean(bet,                     na.rm = TRUE),
    bet_rate_sd    = sd(as.numeric(bet),            na.rm = TRUE),
    stake_mean     = mean(stake[bet],               na.rm = TRUE),
    stake_sd       = sd(stake[bet],                 na.rm = TRUE),
    screen_ms_mean = mean(screen_ms,                na.rm = TRUE),
    screen_ms_sd   = sd(screen_ms,                  na.rm = TRUE),
    side_H_rate    = mean(side[bet] == heads_label, na.rm = TRUE)
  ), by = .(pid, treat)]
  
  # Merge behavioral aggregates — use treat-specific merge since
  # participants only appear in one treatment
  dt <- merge(dt, behav, by = c("pid", "treat"), all.x = TRUE)
  
  # ---- Save per-participant file ----
  if (!skip_pid) {
    fwrite(dt, f_out_pid)
    msg("Saved: ", f_out_pid)
  }
  
  # ---- Group-level summary ----
  one_row <- function(x, group) {
    if (nrow(x) == 0) return(NULL)
    x[, .(
      group                = group,
      N                    = uniqueN(pid),
      age_mean             = mean(age,            na.rm = TRUE),
      age_sd               = sd(age,              na.rm = TRUE),
      sex_F_n              = sum(sex == "F"),
      sex_M_n              = sum(sex == "M"),
      sex_O_n              = sum(sex == "O"),
      payoff_mean          = mean(payoff,          na.rm = TRUE),
      payoff_sd            = sd(payoff,            na.rm = TRUE),
      lotr_mean            = mean(lotr_score,      na.rm = TRUE),
      lotr_sd              = sd(lotr_score,        na.rm = TRUE),
      r_bayes_mean         = mean(r_mean,          na.rm = TRUE),
      r_bayes_sd           = sd(r_mean,            na.rm = TRUE),
      r_classic_mean       = mean(r_hl,            na.rm = TRUE),
      r_classic_sd         = sd(r_hl,              na.rm = TRUE),
      hl_inconsistent_n    = sum(inconsistent == 1L, na.rm = TRUE),
      hl_inconsistent_frac = mean(inconsistent == 1L, na.rm = TRUE),
      bet_rate_mean        = mean(bet_rate_mean,   na.rm = TRUE),
      bet_rate_sd          = mean(bet_rate_sd,     na.rm = TRUE),
      stake_mean           = mean(stake_mean,      na.rm = TRUE),
      stake_sd             = mean(stake_sd,        na.rm = TRUE),
      screen_ms_mean       = mean(screen_ms_mean,  na.rm = TRUE),
      screen_ms_sd         = mean(screen_ms_sd,    na.rm = TRUE),
      side_H_rate          = mean(side_H_rate,     na.rm = TRUE)
    )]
  }
  
  build_group <- function(d, prefix) {
    rows[[prefix]]                   <<- one_row(d,             prefix)
    rows[[paste0(prefix, "/Women")]] <<- one_row(d[sex == "F"], paste0(prefix, "/Women"))
    rows[[paste0(prefix, "/Men")]]   <<- one_row(d[sex == "M"], paste0(prefix, "/Men"))
    
    if (consistent_only) {
      d_cons <- d[inconsistent == 0L]
      rows[[paste0(prefix, "/consistent")]]        <<- one_row(d_cons,             paste0(prefix, "/consistent"))
      rows[[paste0(prefix, "/consistent/Women")]]  <<- one_row(d_cons[sex == "F"], paste0(prefix, "/consistent/Women"))
      rows[[paste0(prefix, "/consistent/Men")]]    <<- one_row(d_cons[sex == "M"], paste0(prefix, "/consistent/Men"))
    }
  }
  
  if (!skip_sum) {
    rows <- list()
    build_group(dt, "pooled")
    for (tr in tr_vec) build_group(dt[treat == tr], tr)
    
    out_sum <- rbindlist(Filter(Negate(is.null), rows), use.names = TRUE, fill = TRUE)
    fwrite(out_sum, f_out_sum)
    msg("Saved: ", f_out_sum)
  }
  
  invisible(TRUE)
}