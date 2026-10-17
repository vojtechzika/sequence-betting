# ============================================================
# 11_table-participants.R
#
# PURPOSE
#   Extracts participant-level variables from the merged oTree export,
#   renames columns to short pipeline names, and saves participants.csv.
#
# INPUT
#   path_src/raw-data.csv
#
# OUTPUT
#   path_src/participants.csv
#
# NOTES
#   - One row per participant (participant.code); uniqueness is
#     guaranteed by 01_read-merge-filter.R
#   - Stops if any expected column is missing from raw-data.csv
# ============================================================

table_participants <- function(cfg) {
  
  infile <- file.path(path_src, "raw-data.csv")
  stopifnot(file.exists(infile))
  
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # ------------------------------------------------------------------
  # Participant-level variables of the raw export (right) and the short
  # pipeline names used in participants.csv (left). The order defines
  # the column order of the output.
  # ------------------------------------------------------------------
  vars <- c(
    # --- Core participant info ---
    pid              = "participant.code",
    label            = "participant.label",
    started_utc      = "participant.time_started_utc",
    payoff           = "participant.payoff",
    treat            = "participant.treatment",
    mult             = "participant.multiplier",
    session          = "session.code",
    exp_version      = "participant.experiment_version",
    session_name     = "session.config.name",
    fee              = "session.config.participation_fee",
    age              = "intro.1.player.age",
    sex              = "intro.1.player.sex",

    # --- Payout information ---
    payout_id        = "payout.1.player.id_in_group",
    payout_role      = "payout.1.player.role",
    payout_payoff    = "payout.1.player.payoff",
    paid_seq_round   = "payout.1.player.paid_seq_round",
    paid_seq_amount  = "payout.1.player.paid_seq_amount",
    paid_mpl_row     = "payout.1.player.paid_mpl_row",
    paid_mpl_choice  = "payout.1.player.paid_mpl_choice",
    paid_mpl_amount  = "payout.1.player.paid_mpl_amount",
    total_paid       = "payout.1.player.total_paid"
  )
  
  # Ensure all expected columns exist
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected columns: ", paste(missing, collapse = ", "))
  }
  
  p <- dt[, unname(vars), with = FALSE]
  data.table::setnames(p, unname(vars), names(vars))
  
  outfile <- file.path(path_src, "participants.csv")
  data.table::fwrite(p, outfile)
  
  msg("Saved: ", outfile)
  invisible(p)
}