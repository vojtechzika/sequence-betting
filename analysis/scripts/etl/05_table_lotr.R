# ============================================================
# 05_table_lotr.R
#
# PURPOSE
#   Extracts LOT-R questionnaire responses from the merged oTree
#   export and saves lotr.csv with one row per participant.
#
# INPUT
#   path_src/merged.csv
#
# OUTPUT
#   path_src/lotr.csv
#
# NOTES
#   - Extracts 10 LOT-R items (q1..q10)
#   - Items are stored as numeric Likert scale values
#   - Scoring (reverse coding, sum) is done in the indices stage
#   - One row per participant (participant.code)
# ============================================================

make_lotr_table <- function(cfg) {
  
  infile <- file.path(path_src, "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  vars <- c(
    "participant.code",
    paste0("lotr.1.player.lotr_", 1:10)
  )
  
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected LOT-R columns: ", paste(missing, collapse = ", "))
  }
  
  lotr <- dt[, ..vars]
  data.table::setnames(lotr, "participant.code", "pid")
  lotr <- unique(lotr, by = "pid")
  
  # Rename item columns to q1..q10
  rename_map <- setNames(
    paste0("q", 1:10),
    paste0("lotr.1.player.lotr_", 1:10)
  )
  old <- names(rename_map)[names(rename_map) %in% names(lotr)]
  data.table::setnames(lotr, old, rename_map[old])
  
  # Coerce to numeric
  for (nm in paste0("q", 1:10)) {
    if (nm %in% names(lotr)) {
      lotr[, (nm) := suppressWarnings(as.numeric(trimws(as.character(get(nm)))))]
    }
  }
  
  outfile <- file.path(path_src, "lotr.csv")
  data.table::fwrite(lotr, outfile)
  
  msg("LOT-R table saved:", outfile, "| rows:", nrow(lotr))
  invisible(lotr)
}