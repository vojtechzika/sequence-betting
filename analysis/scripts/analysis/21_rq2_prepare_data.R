# scripts/analysis/21_rq2_prepare_data.R
source(here::here("scripts", "00_setup.R"))

library(data.table)

rq2_prepare_data <- function(dataset = "pilot",
                             treat_fn = "m25",
                             sd_floor = 2) {
  
  infile <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  req <- c("pid", "treat", "seq", "stake")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("master_sequences.csv missing: ", paste(miss, collapse = ", "))
  
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  
  # FN only
  d <- dt[treat == treat_fn]
  if (nrow(d) == 0) stop("No rows for treat_fn=", treat_fn)
  
  # Betting trials only (RQ2 is intensive margin)
  d <- d[stake > 0]
  if (nrow(d) == 0) stop("No betting trials in FN. Cannot run RQ2.")
  
  # Participant-level summaries across betting trials
  per_i <- d[, .(
    n_bet = .N,
    a_bar = mean(stake),
    a_sd  = sd(stake)
  ), by = pid]
  
  # Exclude <3 betting trials (prereg)
  per_i <- per_i[n_bet >= 3]
  if (nrow(per_i) == 0) stop("No participants with >=3 betting trials in FN.")
  
  per_i[, s_star := pmax(a_sd, sd_floor)]
  
  # Keep only eligible participants
  d <- merge(d, per_i[, .(pid, n_bet, a_bar, a_sd, s_star)], by = "pid", all.x = FALSE)
  
  # Z = (a_is - a_bar_i) / s_star_i
  d[, z := (stake - a_bar) / s_star]
  
  # Levels
  pid_levels <- sort(unique(d$pid))
  seq_levels <- sort(unique(d$seq))
  
  d[, pid_i := match(pid, pid_levels)]
  d[, sid_s := match(seq, seq_levels)]
  
  out <- list(
    d = d,
    pid_levels = pid_levels,
    seq_levels = seq_levels,
    per_i = per_i[pid %in% pid_levels],
    params = list(sd_floor = sd_floor, treat_fn = treat_fn)
  )
  
  f_out <- file.path(path_mod_ds(dataset), paste0("rq2_data_", treat_fn, ".rds"))
  saveRDS(out, f_out)
  msg("Saved RQ2 prepared data: ", f_out)
  
  invisible(out)
}

# Example:
# rq2_prepare_data("pilot")
# rq2_prepare_data("main")