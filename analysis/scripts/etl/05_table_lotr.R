# scripts/05_table_lotr.R
source(here::here("scripts", "00_setup.R"))

make_lotr_table <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # ---- Define LOT-R fields (single round) ----
  vars <- c(
    "participant.code",
    # "lotr.1.player.id_in_group",
    # "lotr.1.player.role",
    # "lotr.1.player.payoff",
    paste0("lotr.1.player.lotr_", 1:10)
  )
  
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected LOT-R columns: ", paste(missing, collapse = ", "))
  }
  
  lotr <- dt[, ..vars]
  
  # One row per participant
  data.table::setnames(lotr, "participant.code", "pid")
  lotr <- unique(lotr, by = "pid")
  
  # ---- Rename columns ----
  rename_map <- c(
    "lotr.1.player.id_in_group" = "id_in_group",
    "lotr.1.player.role"        = "role",
    "lotr.1.player.payoff"      = "payoff"
  )
  
  for (i in 1:10) {
    rename_map[paste0("lotr.1.player.lotr_", i)] <- paste0("q", i)
  }
  
  old <- names(rename_map)
  old <- old[old %in% names(lotr)]
  data.table::setnames(lotr, old, rename_map[old])
  
  # ---- Clean responses (keep numeric Likert scale) ----
  to_num <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "None", "null", "NULL")] <- NA
    suppressWarnings(as.numeric(x))
  }
  
  for (i in 1:10) {
    nm <- paste0("q", i)
    if (nm %in% names(lotr)) {
      lotr[, (nm) := to_num(get(nm))]
    }
  }
  
  if ("payoff" %in% names(lotr)) {
    lotr[, payoff := to_num(payoff)]
  }
  
  # ---- Order columns ----
  ordered <- c("pid", paste0("q", 1:10), "payoff")
  ordered <- ordered[ordered %in% names(lotr)]
  data.table::setcolorder(lotr, ordered)
  
  # ---- Save ----
  outfile <- file.path(path_clean_ds(dataset), "lotr.csv")
  data.table::fwrite(lotr, outfile)
  
  msg("LOT-R table saved:", outfile, "| rows:", nrow(lotr))
  invisible(lotr)
}

# Example:
# source("scripts/05_table_lotr.R")
# make_lotr_table("pilot")