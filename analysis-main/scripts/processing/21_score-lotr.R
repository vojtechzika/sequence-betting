# ============================================================
# 21_score-lotr.R
#
# PURPOSE
#   Scores the LOT-R questionnaire: reverse-codes specified items,
#   computes total score, mean, and z-score per participant.
#   Also computes Cronbach alpha as a reliability check, for all
#   participants and separately for each treatment.
#
# INPUT
#   path_src/lotr.csv
#   path_src/participants.csv            treatment of each participant
#
# OUTPUT
#   path_der/lotr-scored.csv     -- pid, lotr_score, lotr_mean, lotr_z
#   path_csv/processing_lotr-alpha.csv -- reliability summary (Cronbach alpha):
#                                        one row for all participants (pooled)
#                                        and one row per treatment
#
# NOTES
#   - Scoring parameters defined in cfg$design$lotr
#   - Existing output files are overwritten on every run
# ============================================================

score_lotr <- function(cfg) {
  
  lotr_cfg <- cfg$design$lotr
  
  LOTR_SCALE_MIN <- as.integer(lotr_cfg$scale_min)
  LOTR_SCALE_MAX <- as.integer(lotr_cfg$scale_max)
  item_cols      <- as.character(lotr_cfg$item_cols)
  scored_items   <- as.character(lotr_cfg$scored_items)
  rev_items      <- as.character(lotr_cfg$rev_items)
  
  stopifnot(
    length(LOTR_SCALE_MIN) == 1L,
    length(LOTR_SCALE_MAX) == 1L,
    LOTR_SCALE_MAX > LOTR_SCALE_MIN,
    length(item_cols) > 0L,
    length(scored_items) > 0L,
    all(rev_items %in% scored_items),
    all(scored_items %in% item_cols)
  )
  
  infile <- file.path(path_src, "lotr.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", item_cols) %in% names(dt)))
  
  outfile_scored    <- file.path(path_der, "lotr-scored.csv")
  outfile_alpha_csv <- file.path(path_csv, "processing_lotr-alpha.csv")
  
  # ---- Cleaning ----
  dt[, pid := as.character(pid)]
  dt[, (item_cols) := lapply(.SD, as.numeric), .SDcols = item_cols]
  
  if (!dt[, all(complete.cases(.SD)), .SDcols = scored_items]) {
    bad <- dt[!complete.cases(dt[, ..scored_items]), .(pid)]
    stop("LOT-R: missing values in scored items for pid(s): ",
         paste(bad$pid, collapse = ", "))
  }
  
  # ---- Scoring ----
  dt[, (rev_items) := lapply(.SD, function(x) LOTR_SCALE_MAX + LOTR_SCALE_MIN - x),
     .SDcols = rev_items]
  
  dt[, lotr_score := rowSums(.SD),  .SDcols = scored_items]
  dt[, lotr_mean  := rowMeans(.SD), .SDcols = scored_items]
  dt[, lotr_z     := as.numeric(scale(lotr_score))]
  
  out_scored <- dt[, .(pid, lotr_score, lotr_mean, lotr_z)]
  fwrite(out_scored, outfile_scored)
  msg("Saved: ", outfile_scored)
  
  # ---- Cronbach alpha: all participants and each treatment ----
  par <- fread(file.path(path_src, "participants.csv"), encoding = "UTF-8")
  par[, pid := as.character(pid)]
  dt  <- merge(dt, par[, .(pid, treat = as.character(treat))], by = "pid", all.x = TRUE)
  
  interpret_alpha <- function(a) {
    if      (a >= 0.9) "Excellent internal consistency (may be inflated in small samples)."
    else if (a >= 0.8) "Good internal consistency."
    else if (a >= 0.7) "Acceptable internal consistency."
    else if (a >= 0.6) "Questionable internal consistency."
    else if (a >= 0.5) "Poor internal consistency."
    else               "Unacceptable internal consistency."
  }
  
  alpha_row <- function(d, label) {
    res <- suppressWarnings(psych::alpha(d[, ..scored_items]))
    data.table(
      sample                          = label,
      n_participants                  = nrow(d),
      n_items                         = length(scored_items),
      raw_alpha                       = round(res$total$raw_alpha, 4),
      standardized_alpha              = round(res$total$std.alpha, 4),
      average_inter_item_correlation  = round(res$total$average_r, 4),
      interpretation                  = interpret_alpha(res$total$raw_alpha)
    )
  }
  
  alpha_table <- rbindlist(c(
    list(alpha_row(dt, "pooled")),
    lapply(unique(na.omit(dt$treat)), function(tr) alpha_row(dt[treat == tr], tr))
  ))
  
  fwrite(alpha_table, outfile_alpha_csv)
  msg("Saved: ", outfile_alpha_csv)
  
  invisible(list(scored = out_scored, alpha_table = alpha_table))
}