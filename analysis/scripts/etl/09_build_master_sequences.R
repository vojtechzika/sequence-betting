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