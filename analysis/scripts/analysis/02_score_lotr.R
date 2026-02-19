# scripts/analysis/02_score_lotr.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(psych)

score_lotr <- function(dataset = "pilot") {
  
  # =================================================
  # DESIGN PARAMETERS (edit here if LOT-R design changes)
  # =================================================
  
  LOTR_SCALE_MIN <- 0
  LOTR_SCALE_MAX <- 4
  
  # Expected raw item columns in lotr.csv
  item_cols <- paste0("q", 1:10)
  
  # LOT-R key (your screenshot)
  # Scored items: 1, 3, 4, 7, 9, 10
  # Reverse-coded: 3, 7, 9
  scored_items <- c("q1", "q3", "q4", "q7", "q9", "q10")
  rev_items    <- c("q3", "q7", "q9")
  
  # =================================================
  # INPUT
  # =================================================
  
  infile <- file.path(path_clean_ds(dataset), "lotr.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  stopifnot(all(c("pid", item_cols) %in% names(dt)))
  
  # =================================================
  # OUTPUT FOLDERS (dataset-specific artifacts)
  # =================================================
  
  out_dir <- path_out_ds(dataset)
  mod_dir <- path_mod_ds(dataset)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  # =================================================
  # CLEANING
  # =================================================
  
  # Force numeric (guards will catch NAs if something is not numeric)
  dt[, (item_cols) := lapply(.SD, as.numeric), .SDcols = item_cols]
  
  # Guard: no missing in scored items (since responses are forced)
  if (!dt[, all(complete.cases(.SD)), .SDcols = scored_items]) {
    bad <- dt[!complete.cases(dt[, ..scored_items]), .(pid)]
    stop(
      "LOT-R: missing values detected in scored items for pid(s): ",
      paste(bad$pid, collapse = ", ")
    )
  }
  
  # =================================================
  # SCORING
  # =================================================
  
  # Reverse code: x -> max+min-x
  dt[, (rev_items) := lapply(.SD, function(x) LOTR_SCALE_MAX + LOTR_SCALE_MIN - x),
     .SDcols = rev_items]
  
  dt[, lotr_score := rowSums(.SD), .SDcols = scored_items]
  dt[, lotr_mean  := rowMeans(.SD), .SDcols = scored_items]
  dt[, lotr_z     := as.numeric(scale(lotr_score))]
  
  # =================================================
  # RELIABILITY: Cronbach alpha (scored items only)
  # =================================================
  
  alpha_result <- psych::alpha(dt[, ..scored_items])
  alpha_val    <- alpha_result$total$raw_alpha
  
  msg("Cronbach alpha (LOT-R scored items):", round(alpha_val, 3))
  
  # Save full alpha object (models/; should be gitignored)
  saveRDS(alpha_result, file.path(path_mod_ds(dataset), "lotr_alpha.rds"))
  
  # Interpretation (simple)
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
  
  alpha_table <- data.table::data.table(
    dataset = dataset,
    n_participants = nrow(dt),
    n_items = length(scored_items),
    raw_alpha = round(alpha_val, 4),
    standardized_alpha = round(alpha_result$total$std.alpha, 4),
    average_inter_item_correlation = round(alpha_result$total$average_r, 4),
    interpretation = alpha_interpretation
  )
  
  outfile_alpha_csv <- file.path(path_mod_ds(dataset), "lotr_alpha.csv")
  data.table::fwrite(alpha_table, outfile_alpha_csv)
  
  msg("LOT-R alpha table saved:", outfile_alpha_csv)
  
  # =================================================
  # OUTPUT: scored dataset (table lives in data/clean/<ds>/)
  # =================================================
  
  out_scored <- dt[, .(pid, lotr_score, lotr_mean, lotr_z)]
  
  outfile <- file.path(path_clean_ds(dataset), "lotr_scored.csv")
  data.table::fwrite(out_scored, outfile)
  
  msg("LOT-R scored table saved:", outfile, "| rows:", nrow(out_scored))
  
  # Run log (lightweight)
  logfile <- file.path(out_dir, "lotr_run_log.txt")
  cat(
    paste0(
      "dataset=", dataset, "\n",
      "timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "items=", paste(item_cols, collapse = ","), "\n",
      "scored_items=", paste(scored_items, collapse = ","), "\n",
      "rev_items=", paste(rev_items, collapse = ","), "\n",
      "scale_min=", LOTR_SCALE_MIN, " scale_max=", LOTR_SCALE_MAX, "\n",
      "alpha_raw=", round(alpha_val, 6), "\n"
    ),
    file = logfile
  )
  
  invisible(list(scored = out_scored, alpha = alpha_result, alpha_table = alpha_table))
}