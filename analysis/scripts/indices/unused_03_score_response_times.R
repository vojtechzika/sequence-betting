# ============================================================
# scripts/analysis/03_score_response_times.R
#   Response-time scoring (per dataset)
#
# What this does:
# 1) Loads sequences.csv (must contain: pid, stake, screen_ms)
# 2) Creates two RT aggregates per participant:
#      - rt_bet_*   : betting trials (stake > 0)
#      - rt_nobet_* : non-betting trials (stake == 0)
#    using mean of log(screen_ms) (log-mean), plus mean ms.
# 3) Standardizes (z) log-mean RT across participants separately for bet/nobet.
# 4) Writes: data/clean/<ds>/output/response_times.csv
#
# Notes:
# - Keeps participants even if they only have bet or only nobet trials.
# - screen_ms must be > 0 to be used; non-finite/<=0 are dropped.
# ============================================================

library(data.table)

score_response_times <- function(cfg) {
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$run$dataset))
  ds <- as.character(cfg$run$dataset)
  
  in_file <- file.path(path_clean_ds(ds), "sequences.csv")
  stopifnot(file.exists(in_file))
  
  out_dir <- path_clean_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, "response_times.csv")
  
  dt <- fread(in_file, encoding = "UTF-8")
  
  required <- c("pid", "stake", "screen_ms")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("sequences.csv missing columns: ", paste(missing, collapse = ", "),
         "\nExpected at least: ", paste(required, collapse = ", "))
  }
  
  dt[, pid := as.character(pid)]
  dt[, stake := as.numeric(stake)]
  dt[, screen_ms := as.numeric(screen_ms)]
  
  dt[is.na(stake), stake := 0]
  dt <- dt[is.finite(screen_ms) & screen_ms > 0]
  
  if (nrow(dt) == 0) stop("No valid RT rows after filtering screen_ms > 0 (dataset='", ds, "').")
  
  dt[, is_bet := (is.finite(stake) & stake > 0)]
  dt[, log_rt := log(screen_ms)]
  
  # --- aggregates ---
  # --- aggregates (ALL trials with valid screen_ms) ---
  agg_all <- dt[, .(
    n_all = .N,
    rt_all_mean_ms   = mean(screen_ms),
    rt_all_median_ms = median(screen_ms),
    rt_all_log_mean  = mean(log_rt)
  ), by = pid]
  
  agg_bet <- dt[is_bet == TRUE, .(
    n_bet = .N,
    rt_bet_mean_ms = mean(screen_ms),
    rt_bet_median_ms = median(screen_ms),
    rt_bet_log_mean = mean(log_rt)
  ), by = pid]
  
  agg_nb <- dt[is_bet == FALSE, .(
    n_nobet = .N,
    rt_nobet_mean_ms = mean(screen_ms),
    rt_nobet_median_ms = median(screen_ms),
    rt_nobet_log_mean = mean(log_rt)
  ), by = pid]
  
  out <- merge(agg_all, agg_bet, by = "pid", all = TRUE)
  out <- merge(out, agg_nb, by = "pid", all = TRUE)
  out[, dataset := ds]
  
  # --- z-scoring helper (safe when sd=0 or all NA) ---
  z_safe <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
    (x - m) / s
  }
  
  out[, rt_all_z   := z_safe(rt_all_log_mean)]
  out[, rt_bet_z   := z_safe(rt_bet_log_mean)]
  out[, rt_nobet_z := z_safe(rt_nobet_log_mean)]
  
  # After all merges and z-scoring
  setcolorder(out, c(
    "dataset", "pid",
    
    # ALL first
    "n_all", "rt_all_mean_ms", "rt_all_median_ms", "rt_all_log_mean", "rt_all_z",
    
    # BET second
    "n_bet", "rt_bet_mean_ms", "rt_bet_median_ms", "rt_bet_log_mean", "rt_bet_z",
    
    # NOBET third
    "n_nobet", "rt_nobet_mean_ms", "rt_nobet_median_ms", "rt_nobet_log_mean", "rt_nobet_z"
  ))
  
  setorder(out, pid)
  
  fwrite(out, out_file)
  msg("Saved: ", out_file)
  
  invisible(out)
}

# Example:
# score_response_times(cfg)