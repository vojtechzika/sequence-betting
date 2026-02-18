# scripts/02_table_participants.R
source(here::here("scripts", "00_setup.R"))

make_participants_table <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
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
    
    # Example of excluding:
    # "some_unused_field"
  )
  
  # Ensure all active vars exist
  missing <- setdiff(vars, names(dt))
  if (length(missing) > 0) {
    stop("Missing expected columns: ", paste(missing, collapse = ", "))
  }
  
  # Keep only selected columns
  p <- dt[, ..vars]
  
  # One row per participant
  p <- unique(p, by = "participant.code")
  
  # ------------------------------------------------------------------
  # Rename (this does NOT define schema; vars does)
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
    
    # payout
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
  
  # ------------------------------------------------------------------
  # Keep final columns derived ONLY from vars
  # ------------------------------------------------------------------
  
  final_vars <- vars
  final_vars <- ifelse(final_vars %in% names(rename_map),
                       rename_map[final_vars],
                       final_vars)
  
  final_vars <- final_vars[final_vars %in% names(p)]
  p <- p[, ..final_vars]
  
  # ------------------------------------------------------------------
  
  outfile <- file.path(path_clean_ds(dataset), "participants.csv")
  data.table::fwrite(p, outfile)
  
  msg("Participants table saved:", outfile, "| rows:", nrow(p))
  invisible(p)
}

# make_participants_table("pilot")