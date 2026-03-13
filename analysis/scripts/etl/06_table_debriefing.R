# scripts/06_table_debriefing.R
source(here::here("scripts", "00_setup.R"))

make_debriefing_table <- function(cfg) {
  dataset <- as.character(cfg$run$dataset)
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # ---- Detect all debriefing fields automatically ----
  deb_cols <- grep("^debriefing\\.1\\.player\\.", names(dt), value = TRUE)
  
  vars <- c("participant.code", deb_cols)
  
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected debriefing columns: ", paste(missing, collapse = ", "))
  }
  
  deb <- dt[, ..vars]
  
  # One row per participant
  data.table::setnames(deb, "participant.code", "pid")
  deb <- unique(deb, by = "pid")
  
  # ---- Build rename map automatically ----
  # debriefing.1.player.X  ->  X
  rename_map <- setNames(
    sub("^debriefing\\.1\\.player\\.", "", deb_cols),
    deb_cols
  )
  
  old <- intersect(names(rename_map), names(deb))
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