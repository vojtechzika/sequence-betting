# ============================================================
# 13_table-holt-laury.R
#
# PURPOSE
#   Extracts Holt-Laury MPL choices from the merged oTree export
#   and saves mpl.csv with one row per participant.
#
# INPUT
#   path_src/raw-data.csv
#
# OUTPUT
#   path_src/holt-laury.csv
#
# NOTES
#   - Extracts K choice columns as defined in cfg$design$mpl$K
#   - Choices are kept as character (A/B labels from oTree) and
#     validated: any missing or other label stops the script
#   - One row per participant (pid = participant.code)
# ============================================================

table_holt_laury <- function(cfg) {
  
  infile <- file.path(path_src, "raw-data.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  K <- cfg$design$mpl$K
  
  # output name = oTree column name
  vars <- c(
    pid = "participant.code",
    stats::setNames(paste0("mpl.1.player.choice_", seq_len(K)),
                    paste0("c", seq_len(K)))
  )
  
  missing <- setdiff(unname(vars), names(dt))
  if (length(missing) > 0) {
    stop("Missing expected MPL columns: ", paste(missing, collapse = ", "))
  }
  
  mpl <- dt[, unname(vars), with = FALSE]
  data.table::setnames(mpl, unname(vars), names(vars))
  
  if (anyDuplicated(mpl$pid)) stop("Repeated pid in MPL table.")
  
  choices <- unlist(mpl[, paste0("c", seq_len(K)), with = FALSE], use.names = FALSE)
  if (anyNA(choices) || !all(choices %in% c("A", "B"))) {
    stop("MPL choices must be non-missing 'A' or 'B'.")
  }
  
  outfile <- file.path(path_src, "holt-laury.csv")
  data.table::fwrite(mpl, outfile)
  
  msg("Saved: ", outfile)
  invisible(mpl)
}
