# scripts/04_table_mpl.R
source(here::here("scripts", "00_setup.R"))

make_mpl_table <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  keep <- c(
     "participant.code",
    # "mpl.1.player.id_in_group",
    # "mpl.1.player.role",
    # "mpl.1.player.payoff",
    "mpl.1.player.choice_1",
    "mpl.1.player.choice_2",
    "mpl.1.player.choice_3",
    "mpl.1.player.choice_4",
    "mpl.1.player.choice_5",
    "mpl.1.player.choice_6",
    "mpl.1.player.choice_7",
    "mpl.1.player.choice_8",
    "mpl.1.player.choice_9",
    "mpl.1.player.choice_10"
    # "mpl.1.player.random_draw",
    # "mpl.1.player.choice_to_pay",
    # "mpl.1.player.option_to_pay",
    # "mpl.1.player.inconsistent",
    # "mpl.1.player.switching_row",
    # "mpl.1.player.round_earnings"
  )
  
  missing <- setdiff(keep, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected MPL columns: ", paste(missing, collapse = ", "))
  }
  
  mpl <- dt[, ..keep]
  
  # one row per participant
  data.table::setnames(mpl, "participant.code", "pid")
  mpl <- unique(mpl, by = "pid")
  
  # rename to shorter names
  rename_map <- c(
    "mpl.1.player.id_in_group"     = "id_in_group",
    "mpl.1.player.role"            = "role",
    "mpl.1.player.payoff"          = "payoff",
    "mpl.1.player.random_draw"     = "random_draw",
    "mpl.1.player.choice_to_pay"   = "choice_to_pay",
    "mpl.1.player.option_to_pay"   = "option_to_pay",
    "mpl.1.player.inconsistent"    = "inconsistent",
    "mpl.1.player.switching_row"   = "switch_row",
    "mpl.1.player.round_earnings"  = "earn"
  )
  # choices 1..10
  for (i in 1:10) {
    rename_map[paste0("mpl.1.player.choice_", i)] <- paste0("c", i)
  }
  
  old <- names(rename_map)
  old <- old[old %in% names(mpl)]
  data.table::setnames(mpl, old, rename_map[old])
  
  # Optional: coerce choices to integer (often 0/1 or "A"/"B")
  # Keep as character unless you are sure; here we keep as character and normalize blanks.
  for (nm in c(paste0("c", 1:10), "random_draw", "choice_to_pay", "option_to_pay", "switch_row")) {
    if (nm %in% names(mpl)) {
      mpl[, (nm) := trimws(as.character(get(nm)))]
      mpl[get(nm) == "", (nm) := NA_character_]
    }
  }
  
  # numeric-ish fields (safe)
  to_num <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "None", "null", "NULL")] <- NA
    suppressWarnings(as.numeric(x))
  }
  for (nm in c("payoff", "earn")) {
    if (nm %in% names(mpl)) mpl[, (nm) := to_num(get(nm))]
  }
  for (nm in c("inconsistent")) {
    if (nm %in% names(mpl)) mpl[, (nm) := to_num(get(nm))]
  }
  
  # order columns
  ordered <- c("pid", "id_in_group", "role", paste0("c", 1:10),
               "random_draw", "choice_to_pay", "option_to_pay",
               "inconsistent", "switch_row", "payoff", "earn")
  ordered <- ordered[ordered %in% names(mpl)]
  data.table::setcolorder(mpl, ordered)
  
  outfile <- file.path(path_clean_ds(dataset), "mpl.csv")
  data.table::fwrite(mpl, outfile)
  
  msg("MPL table saved:", outfile, "| rows:", nrow(mpl))
  invisible(mpl)
}

# Example:
# source("scripts/04_table_mpl.R")
# make_mpl_table("pilot")