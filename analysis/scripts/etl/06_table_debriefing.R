# ============================================================
# 06_table_debriefing.R
#
# PURPOSE
#   Extracts debriefing questionnaire responses from the merged
#   oTree export and saves debriefing.csv with one row per participant.
#
# INPUT
#   path_src/merged.csv
#
# OUTPUT
#   path_src/debriefing.csv
#
# NOTES
#   - Debriefing columns are detected automatically via regex
#   - Numeric Likert-style variables are coerced to numeric
#   - Free-text variables are kept as character
#   - One row per participant (participant.code)
# ============================================================

make_debriefing_table <- function(cfg) {
  
  infile <- file.path(path_src, "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # Detect all debriefing fields automatically
  deb_cols <- grep("^debriefing\\.1\\.player\\.", names(dt), value = TRUE)
  
  vars <- c("participant.code", deb_cols)
  
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected debriefing columns: ", paste(missing, collapse = ", "))
  }
  
  deb <- dt[, ..vars]
  data.table::setnames(deb, "participant.code", "pid")
  deb <- unique(deb, by = "pid")
  
  # Rename: debriefing.1.player.X -> X
  rename_map <- setNames(
    sub("^debriefing\\.1\\.player\\.", "", deb_cols),
    deb_cols
  )
  old <- intersect(names(rename_map), names(deb))
  data.table::setnames(deb, old, rename_map[old])
  
  # Coerce numeric Likert-style variables
  numeric_vars <- intersect(
    c("belief_independence", "reliance_on_sequence", "perceived_realism",
      "action_seeking", "enjoyment", "fatigue", "self_risk_tolerance", "payoff"),
    names(deb)
  )
  for (nm in numeric_vars) {
    deb[, (nm) := suppressWarnings(as.numeric(trimws(as.character(get(nm)))))]
  }
  
  # Keep free-text as character
  text_vars <- intersect(c("strategy_text", "comment"), names(deb))
  for (nm in text_vars) {
    deb[, (nm) := as.character(get(nm))]
  }
  
  # Order columns
  ordered <- c("pid", "belief_independence", "reliance_on_sequence",
               "perceived_realism", "action_seeking", "enjoyment",
               "fatigue", "self_risk_tolerance", "used_strategy",
               "strategy_text", "comment", "payoff")
  ordered <- ordered[ordered %in% names(deb)]
  data.table::setcolorder(deb, ordered)
  
  outfile <- file.path(path_src, "debriefing.csv")
  data.table::fwrite(deb, outfile)
  
  msg("Debriefing table saved:", outfile, "| rows:", nrow(deb))
  invisible(deb)
}