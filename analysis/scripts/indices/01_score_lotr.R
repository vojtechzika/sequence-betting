# ============================================================
# 01_score_lotr.R
#
# PURPOSE
#   Scores the LOT-R questionnaire: reverse-codes specified items,
#   computes total score, mean, and z-score per participant.
#   Also computes Cronbach alpha as a reliability check.
#
# INPUT
#   path_src/lotr.csv
#
# OUTPUT
#   path_out/lotr_scored.csv     -- pid, lotr_score, lotr_mean, lotr_z
#   path_out/lotr_alpha.csv      -- reliability summary
#   path_mod/lotr_alpha.rds      -- full psych::alpha object
#
# NOTES
#   - Scoring parameters defined in cfg$design$lotr
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
  
  outfile_scored    <- file.path(path_out, "lotr_scored.csv")
  outfile_alpha_rds <- file.path(path_mod, "lotr_alpha.rds")
  outfile_alpha_csv <- file.path(path_out, "lotr_alpha.csv")
  
  if (should_skip(
    paths = outfile_scored,
    cfg   = cfg,
    type  = "output",
    label = "LOT-R scored"
  )) return(invisible(NULL))
  
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
  msg("Saved: ", outfile_scored, " | rows: ", nrow(out_scored))
  
  # ---- Cronbach alpha ----
  alpha_result <- psych::alpha(dt[, ..scored_items])
  alpha_val    <- alpha_result$total$raw_alpha
  
  msg("Cronbach alpha (LOT-R scored items): ", round(alpha_val, 3))
  
  alpha_interpretation <- if (alpha_val >= 0.9) {
    "Excellent internal consistency (may be inflated in small samples)."
  } else if (alpha_val >= 0.8) {
    "Good internal consistency."
  } else if (alpha_val >= 0.7) {
    "Acceptable internal consistency."
  } else if (alpha_val >= 0.6) {
    "Questionable internal consistency."
  } else if (alpha_val >= 0.5) {
    "Poor internal consistency."
  } else {
    "Unacceptable internal consistency."
  }
  
  if (!should_skip(
    paths = outfile_alpha_rds,
    cfg   = cfg,
    type  = "model",
    label = "LOT-R alpha rds"
  )) {
    saveRDS(alpha_result, outfile_alpha_rds)
    msg("Saved: ", outfile_alpha_rds)
  }
  
  alpha_table <- data.table(
    n_participants                  = nrow(dt),
    n_items                         = length(scored_items),
    raw_alpha                       = round(alpha_val, 4),
    standardized_alpha              = round(alpha_result$total$std.alpha, 4),
    average_inter_item_correlation  = round(alpha_result$total$average_r, 4),
    interpretation                  = alpha_interpretation
  )
  
  if (!should_skip(
    paths = outfile_alpha_csv,
    cfg   = cfg,
    type  = "output",
    label = "LOT-R alpha csv"
  )) {
    fwrite(alpha_table, outfile_alpha_csv)
    msg("Saved: ", outfile_alpha_csv)
  }
  
  invisible(list(scored = out_scored, alpha = alpha_result, alpha_table = alpha_table))
}