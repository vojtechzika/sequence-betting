# ============================================================
# 15_table-debriefing.R
#
# PURPOSE
#   Extracts debriefing questionnaire responses from the merged
#   oTree export and saves debriefing.csv with one row per participant.
#
# INPUT
#   path_src/raw-data.csv
#
# OUTPUT
#   path_src/debriefing.csv
#
# NOTES
#   - Columns are listed explicitly in `vars`; a missing column stops
#     the script
#   - Likert-style items and the perceived probability are numeric;
#     missing answers are allowed (NA), non-numeric entries are errors
#   - used_strategy must be one of yes / no / unsure
#   - Free-text variables are kept as character
#   - Excluded oTree fields: id_in_group, role, payoff (constant zero),
#     group.id_in_subsession, subsession.round_number
#   - One row per participant (pid = participant.code)
# ============================================================

table_debriefing <- function(cfg) {
  
  infile <- file.path(path_src, "raw-data.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  numeric_vars <- c("belief_independence", "reliance_on_sequence",
                    "action_seeking", "fatigue", "self_risk_tolerance",
                    "perceived_win_probability", "perceived_risk")
  text_vars <- c("used_strategy", "strategy_text", "comment")
  
  # output name = oTree column name
  items <- c(numeric_vars, text_vars)
  vars <- c(
    pid = "participant.code",
    stats::setNames(paste0("debriefing.1.player.", items), items)
  )
  
  missing <- setdiff(unname(vars), names(dt))
  if (length(missing) > 0) {
    stop("Missing expected debriefing columns: ", paste(missing, collapse = ", "))
  }
  
  deb <- dt[, unname(vars), with = FALSE]
  data.table::setnames(deb, unname(vars), names(vars))
  
  if (anyDuplicated(deb$pid)) stop("Repeated pid in debriefing table.")
  
  # Numeric items: NA allowed, other non-numeric entries are errors
  for (nm in numeric_vars) {
    raw <- trimws(as.character(deb[[nm]]))
    raw[raw == ""] <- NA_character_
    x <- suppressWarnings(as.numeric(raw))
    if (any(is.na(x) & !is.na(raw))) {
      stop("Non-numeric debriefing responses in ", nm, ".")
    }
    data.table::set(deb, j = nm, value = x)
  }
  
  # Free text and categorical items stay character
  for (nm in text_vars) {
    data.table::set(deb, j = nm, value = as.character(deb[[nm]]))
  }
  
  if (!all(deb$used_strategy %in% c("yes", "no", "unsure"))) {
    stop("used_strategy must be one of 'yes', 'no', 'unsure'.")
  }
  
  data.table::setcolorder(deb, c("pid", items))
  
  outfile <- file.path(path_src, "debriefing.csv")
  data.table::fwrite(deb, outfile)
  
  msg("Saved: ", outfile)
  invisible(deb)
}
