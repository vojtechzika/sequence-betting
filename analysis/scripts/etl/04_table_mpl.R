# ============================================================
# 04_table_mpl.R
#
# PURPOSE
#   Extracts Holt-Laury MPL choices from the merged oTree export
#   and saves mpl.csv with one row per participant.
#
# INPUT
#   path_src/merged.csv
#
# OUTPUT
#   path_src/mpl.csv
#
# NOTES
#   - Extracts K choice columns as defined in cfg$design$mpl$K
#   - Choices are kept as character (A/B labels from oTree)
#   - One row per participant (participant.code)
# ============================================================

make_mpl_table <- function(cfg) {
  
  infile <- file.path(path_src, "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  K <- cfg$design$mpl$K
  
  keep <- c(
    "participant.code",
    paste0("mpl.1.player.choice_", seq_len(K))
  )
  
  missing <- setdiff(keep, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected MPL columns: ", paste(missing, collapse = ", "))
  }
  
  mpl <- dt[, ..keep]
  data.table::setnames(mpl, "participant.code", "pid")
  mpl <- unique(mpl, by = "pid")
  
  # Rename choice columns to c1..cK
  rename_map <- setNames(
    paste0("c", seq_len(K)),
    paste0("mpl.1.player.choice_", seq_len(K))
  )
  old <- names(rename_map)
  old <- old[old %in% names(mpl)]
  data.table::setnames(mpl, old, rename_map[old])
  
  outfile <- file.path(path_src, "mpl.csv")
  data.table::fwrite(mpl, outfile)
  
  msg("MPL table saved:", outfile, "| rows:", nrow(mpl))
  invisible(mpl)
}