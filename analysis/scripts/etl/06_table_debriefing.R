# scripts/06_table_debriefing.R
source(here::here("scripts", "00_setup.R"))

make_debriefing_table <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # ---- Define debreifing fields ----
  vars <- c(
    "participant.code",
    "debriefing.1.player.id_in_group",
    "debriefing.1.player.role",
    "debriefing.1.player.payoff",
    "debriefing.1.player.belief_independence",
    "debriefing.1.player.reliance_on_sequence",
    "debriefing.1.player.perceived_realism_history",
    "debriefing.1.player.action_seeking",
    "debriefing.1.player.enjoyment",
    "debriefing.1.player.fatigue",
    "debriefing.1.player.self_risk_tolerance",
    "debriefing.1.player.used_strategy",
    "debriefing.1.player.strategy_text",
    "debriefing.1.player.comment"
  )
  
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected debriefing columns: ", paste(missing, collapse = ", "))
  }
  
  deb <- dt[, ..vars]
  
  # One row per participant
  data.table::setnames(deb, "participant.code", "pid")
  deb <- unique(deb, by = "pid")
  
  # ---- Rename columns ----
  rename_map <- c(
    "debriefing.1.player.id_in_group"          = "id_in_group",
    "debriefing.1.player.role"                 = "role",
    "debriefing.1.player.payoff"               = "payoff",
    "debriefing.1.player.belief_independence"  = "belief_independence",
    "debriefing.1.player.reliance_on_sequence" = "reliance_on_sequence",
    "debriefing.1.player.perceived_realism_history" = "perceived_realism",
    "debriefing.1.player.action_seeking"       = "action_seeking",
    "debriefing.1.player.enjoyment"            = "enjoyment",
    "debriefing.1.player.fatigue"              = "fatigue",
    "debriefing.1.player.self_risk_tolerance"  = "self_risk_tolerance",
    "debriefing.1.player.used_strategy"        = "used_strategy",
    "debriefing.1.player.strategy_text"        = "strategy_text",
    "debriefing.1.player.comment"              = "comment"
  )
  
  old <- names(rename_map)
  old <- old[old %in% names(deb)]
  data.table::setnames(deb, old, rename_map[old])
  
  # ---- Clean numeric Likert-style variables ----
  to_num <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "None", "null", "NULL")] <- NA
    suppressWarnings(as.numeric(x))
  }
  
  numeric_vars <- c(
    "belief_independence",
    "reliance_on_sequence",
    "perceived_realism",
    "action_seeking",
    "enjoyment",
    "fatigue",
    "self_risk_tolerance"
  )
  
  numeric_vars <- numeric_vars[numeric_vars %in% names(deb)]
  
  for (nm in numeric_vars) {
    deb[, (nm) := to_num(get(nm))]
  }
  
  if ("payoff" %in% names(deb)) {
    deb[, payoff := to_num(payoff)]
  }
  
  # Ensure free-text stays character
  text_vars <- c("strategy_text", "comment")
  text_vars <- text_vars[text_vars %in% names(deb)]
  
  for (nm in text_vars) {
    deb[, (nm) := as.character(get(nm))]
  }
  
  # ---- Order columns ----
  ordered <- c(
    "pid",
    numeric_vars,
    "used_strategy",
    "strategy_text",
    "comment",
    "payoff"
  )
  ordered <- ordered[ordered %in% names(deb)]
  data.table::setcolorder(deb, ordered)
  
  # ---- Save ----
  outfile <- file.path(path_clean_ds(dataset), "debriefing.csv")
  data.table::fwrite(deb, outfile)
  
  msg("Debriefing table saved:", outfile, "| rows:", nrow(deb))
  invisible(deb)
}

# Example:
# source("scripts/06_table_debriefing.R")
# make_debriefing_table("pilot")