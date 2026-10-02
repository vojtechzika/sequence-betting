# ============================================================
# 09_build_master_sequences.R
#
# PURPOSE
#   Merges sequence-level and participant-level ETL outputs into a
#   single analysis-ready master file. Applies Czech-to-English
#   label recoding for coin outcome labels.
#
# INPUT
#   path_src/sequences.csv
#   path_src/participants.csv
#
# OUTPUT
#   path_src/master_sequences.csv
#
# NOTES
#   - One row per participant x sequence trial
#   - Participant covariates (age, sex, session, treat, seat) are
#     merged from participants.csv
#   - Czech oTree labels (O) are recoded to English (T) here
#   - Two sequence-level rules are computed from `seq` (the 6-toss
#     string of prior outcomes), independent of the participant's bet:
#       imbalance : mean of the 6 tosses coded H=+1, T=-1.
#                   Ranges -1..+1; 0 when heads/tails counts are equal
#                   (e.g. TTTHHH, HHHTTT); +1 for HHHHHH; -1 for TTTTTT.
#       has_terminal_run2      : sign of the terminal run of length >= 2.
#                   +1 if the sequence ends "...HH", -1 if it ends
#                   "...TT", 0 if the last two tosses differ (no
#                   terminal run of 2+) -- coded 0 rather than NA so
#                   it enters linear models without dropping rows.
#   - Stops if duplicate pids found in participants.csv
# ============================================================
build_master_sequences <- function(cfg) {
  
  f_seq <- file.path(path_src, "sequences.csv")
  f_par <- file.path(path_src, "participants.csv")
  stopifnot(file.exists(f_seq), file.exists(f_par))
  
  seq <- fread(f_seq, encoding = "UTF-8")
  par <- fread(f_par, encoding = "UTF-8")
  stopifnot("pid" %in% names(seq), "pid" %in% names(par))
  
  # Keep only required participant-level columns
  par_keep    <- c("pid", "age", "sex", "session", "treat", "label")
  missing_par <- setdiff(par_keep, names(par))
  if (length(missing_par) > 0) {
    stop("participants.csv missing: ", paste(missing_par, collapse = ", "))
  }
  par <- par[, ..par_keep]
  
  # Enforce 1 row per pid
  dup <- par[, .N, by = pid][N > 1]
  if (nrow(dup) > 0) {
    stop(
      "participants.csv has duplicate pid rows. Example pids: ",
      paste(head(dup$pid, 5), collapse = ", "),
      " (n_dups=", nrow(dup), ")"
    )
  }
  
  # Merge
  master <- merge(seq, par, by = "pid", all.x = TRUE)
  
  # Recode Czech coin labels to English
  recode <- cfg$design$seq$label_recode
  recode_labels <- function(x) {
    for (from in names(recode)) x <- gsub(from, recode[[from]], x, fixed = TRUE)
    x
  }
  if ("seq"       %in% names(master)) master[, seq       := recode_labels(seq)]
  if ("realized"  %in% names(master)) master[, realized  := recode_labels(realized)]
  if ("side"      %in% names(master)) master[, side      := recode_labels(side)]
  if ("btn_order" %in% names(master)) master[, btn_order := recode_labels(btn_order)]
  msg("Label recode applied: ", paste(names(recode), unlist(recode), sep = "->", collapse = ", "))
  
  # --------------------------------------------------------
  # Sequence-level rules (computed once per unique `seq`, merged back)
  # --------------------------------------------------------
  seq_to_x <- function(seq_str) {
    toss <- strsplit(seq_str, "")[[1]]
    ifelse(toss == "H", 1L, ifelse(toss == "T", -1L, NA_integer_))
  }
  
  uniq_seq   <- unique(master$seq)
  seq_lookup <- rbindlist(lapply(uniq_seq, function(s) {
    x     <- seq_to_x(s)
    last2 <- tail(x, 2)              # last two tosses, whatever the sequence length
    imbalance_val    <- mean(x)
    terminal_run_val <- if (last2[1] == last2[2]) last2[2] else 0L
    list(
      seq              = s,
      imbalance        = imbalance_val, 
      imbalance_abs    = abs(imbalance_val),
      terminal_run     = terminal_run_val,       # signed: 1=HH, -1=TT, 0=none
      has_terminal_run = if (terminal_run_val == 0L) 0L else 1L  # unsigned dummy
    )
  }))
  master <- merge(master, seq_lookup, by = "seq", all.x = TRUE)
  msg("Sequence-level rules computed: imbalance, imbalance_abs, terminal_run, has_terminal_run | unique seq: ", length(uniq_seq))
  
  # Diagnostics
  n_missing_par <- master[
    is.na(age) | is.na(sex) | is.na(session) | is.na(treat) | is.na(label),
    uniqueN(pid)
  ]
  if (n_missing_par > 0) {
    warning("Missing participant covariates for ", n_missing_par, " pid(s).")
  }
  
  master[, sex := as.character(sex)]
  setnames(master, "label", "seat")
  
  outfile <- file.path(path_src, "master_sequences.csv")
  fwrite(master, outfile)
  
  msg("Master sequences saved: ", outfile,
      " | rows: ", nrow(master),
      " | unique pid: ", master[, uniqueN(pid)])
  
  invisible(master)
}