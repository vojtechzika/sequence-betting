source("scripts/00_setup.R")

make_participants_table <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # Columns you want to keep (raw names)
  keep <- c(
    "participant.code",
    "participant.label",
    "participant.time_started_utc",
    "participant.payoff",
    "participant.treatment",
    "participant.multiplier",
    "session.code",
    "participant.experiment_version",
    "session.config.name",
    "session.config.participation_fee",
    "intro.1.player.age",
    "intro.1.player.sex"
  )
  
  missing <- setdiff(keep, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected columns: ", paste(missing, collapse = ", "))
  }
  
  # Keep only requested columns
  p <- dt[, ..keep]
  
  # Rename to shorter, analysis-friendly names
  data.table::setnames(
    p,
    old = keep,
    new = c(
      "pid",          # participant.code
      "label",        # participant.label
      "started_utc",  # participant.time_started_utc
      "payoff",       # participant.payoff
      "treat",        # participant.treatment
      "mult",         # participant.multiplier
      "session",      # session.code
      "exp_version",  # participant.experiment_version
      "session_name", # session.config.name
      "fee",          # session.config.participation_fee
      "age",          # intro.1.player.age
      "sex"           # intro.1.player.sex
    )
  )
  
  # Ensure 1 row per participant (if duplicates exist, keep first and warn)
  dup_n <- p[, .N, by = pid][N > 1, sum(N)]
  if (!is.na(dup_n) && dup_n > 0) {
    warning("Duplicate pid rows detected. Keeping first occurrence per pid.")
    data.table::setkey(p, pid)
    p <- unique(p, by = "pid")
  }
  
  outfile <- file.path(path_clean_ds(dataset), "participants.csv")
  data.table::fwrite(p, outfile)
  
  msg("Participants table saved:", outfile, "| rows:", nrow(p))
  invisible(p)
}

# Example run (comment out if you prefer calling manually)
# make_participants_table("pilot")