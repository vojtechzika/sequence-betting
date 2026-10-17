# ============================================================
# 14_table-lotr.R
#
# PURPOSE
#   Extracts LOT-R questionnaire responses from the merged oTree
#   export and saves lotr.csv with one row per participant.
#
# INPUT
#   path_src/raw-data.csv
#
# OUTPUT
#   path_src/lotr.csv
#
# NOTES
#   - Extracts 10 LOT-R items (q1..q10)
#   - Items are stored as numeric Likert scale values; a missing or
#     non-numeric response stops the script
#   - Scoring (reverse coding, sum) is done in 21_score-lotr.R
#   - One row per participant (pid = participant.code)
# ============================================================

table_lotr <- function(cfg) {
  
  infile <- file.path(path_src, "raw-data.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # output name = oTree column name
  items <- paste0("q", 1:10)
  vars <- c(
    pid = "participant.code",
    stats::setNames(paste0("lotr.1.player.lotr_", 1:10), items)
  )
  
  missing <- setdiff(unname(vars), names(dt))
  if (length(missing) > 0) {
    stop("Missing expected LOT-R columns: ", paste(missing, collapse = ", "))
  }
  
  lotr <- dt[, unname(vars), with = FALSE]
  data.table::setnames(lotr, unname(vars), names(vars))
  
  if (anyDuplicated(lotr$pid)) stop("Repeated pid in LOT-R table.")
  
  # Coerce items to numeric; any missing or non-numeric entry is an error
  for (nm in items) {
    x <- suppressWarnings(as.numeric(trimws(as.character(lotr[[nm]]))))
    if (anyNA(x)) stop("Missing or non-numeric LOT-R responses in ", nm, ".")
    data.table::set(lotr, j = nm, value = x)
  }
  
  outfile <- file.path(path_src, "lotr.csv")
  data.table::fwrite(lotr, outfile)
  
  msg("Saved: ", outfile)
  invisible(lotr)
}