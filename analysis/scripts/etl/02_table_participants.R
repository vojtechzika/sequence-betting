# ============================================================
# 02_table_participants.R
#
# PURPOSE
#   Extracts participant-level variables from the merged oTree export,
#   renames columns to short pipeline names, and saves participants.csv.
#
# INPUT
#   path_src/merged.csv
#
# OUTPUT
#   path_src/participants.csv
#
# NOTES
#   - One row per participant (participant.code)
#   - Stops if any expected column is missing from merged.csv
#   - Warns if duplicate participant codes are found before deduplication
# ============================================================

make_participants_table <- function(cfg) {
  
  infile <- file.path(path_src, "merged.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # ------------------------------------------------------------------
  # All participant-level variables available in raw export.
  # Active variables define schema of participants.csv.
  # Commented variables are intentionally excluded.
  # ------------------------------------------------------------------
  vars <- c(
    
    # --- Core participant info ---
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
    "intro.1.player.sex",
    
    # --- Payout information ---
    "payout.1.player.id_in_group",
    "payout.1.player.role",
    "payout.1.player.payoff",
    "payout.1.player.paid_seq_round",
    "payout.1.player.paid_seq_amount",
    "payout.1.player.paid_mpl_row",
    "payout.1.player.paid_mpl_choice",
    "payout.1.player.paid_mpl_amount",
    "payout.1.player.total_paid"
  )
  
  # Ensure all active vars exist
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected columns: ", paste(missing, collapse = ", "))
  }
  
  p <- dt[, ..vars]
  
  # ---- Warn if duplicates exist before deduplication ----
  dup_codes <- p[duplicated(participant.code), unique(participant.code)]
  if (length(dup_codes) > 0) {
    warning(
      length(dup_codes), " participant code(s) have duplicate rows. ",
      "Keeping first occurrence. Check merged.csv for data issues.\n",
      "Affected codes: ", paste(dup_codes, collapse = ", ")
    )
  }
  p <- unique(p, by = "participant.code")
  
  # ------------------------------------------------------------------
  # Rename
  # ------------------------------------------------------------------
  rename_map <- c(
    "participant.code"                  = "pid",
    "participant.label"                 = "label",
    "participant.time_started_utc"      = "started_utc",
    "participant.payoff"                = "payoff",
    "participant.treatment"             = "treat",
    "participant.multiplier"            = "mult",
    "session.code"                      = "session",
    "participant.experiment_version"    = "exp_version",
    "session.config.name"               = "session_name",
    "session.config.participation_fee"  = "fee",
    "intro.1.player.age"                = "age",
    "intro.1.player.sex"                = "sex",
    "payout.1.player.id_in_group"       = "payout_id",
    "payout.1.player.role"              = "payout_role",
    "payout.1.player.payoff"            = "payout_payoff",
    "payout.1.player.paid_seq_round"    = "paid_seq_round",
    "payout.1.player.paid_seq_amount"   = "paid_seq_amount",
    "payout.1.player.paid_mpl_row"      = "paid_mpl_row",
    "payout.1.player.paid_mpl_choice"   = "paid_mpl_choice",
    "payout.1.player.paid_mpl_amount"   = "paid_mpl_amount",
    "payout.1.player.total_paid"        = "total_paid"
  )
  
  old <- names(rename_map)
  old <- old[old %in% names(p)]
  data.table::setnames(p, old, rename_map[old])
  
  # Keep final columns
  final_vars <- vars
  final_vars <- ifelse(final_vars %in% names(rename_map),
                       rename_map[final_vars],
                       final_vars)
  final_vars <- final_vars[final_vars %in% names(p)]
  p <- p[, ..final_vars]
  
  outfile <- file.path(path_src, "participants.csv")
  data.table::fwrite(p, outfile)
  
  msg("Participants table saved:", outfile, "| rows:", nrow(p))
  invisible(p)
}