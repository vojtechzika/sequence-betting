# scripts/analysis/09_build_master_sequences.R
source(here::here("scripts", "00_setup.R"))

library(data.table)

build_master_sequences <- function(cfg) {
  

  dataset <- as.character(cfg$run$dataset)
  
  # ----------------------------
  # Inputs (BASE only)
  # ----------------------------
  f_seq <- file.path(path_clean_ds(dataset), "sequences.csv")
  f_par <- file.path(path_clean_ds(dataset), "participants.csv")
  
  stopifnot(file.exists(f_seq), file.exists(f_par))
  
  seq <- fread(f_seq, encoding = "UTF-8")
  par <- fread(f_par, encoding = "UTF-8")
  
  stopifnot("pid" %in% names(seq), "pid" %in% names(par))
  
  # ----------------------------
  # Keep only requested columns (pid-level)
  # ----------------------------
  par_keep <- c("pid", "age", "sex", "session", "treat", "label")
  missing_par <- setdiff(par_keep, names(par))
  if (length(missing_par) > 0) stop("participants.csv missing: ", paste(missing_par, collapse = ", "))
  
  par <- par[, ..par_keep]
  
  # ----------------------------
  # Enforce 1 row per pid in participants
  # ----------------------------
  dup <- par[, .N, by = pid][N > 1]
  if (nrow(dup) > 0) {
    stop(
      "participants.csv has duplicate pid rows. Example pids: ",
      paste(head(dup$pid, 5), collapse = ", "),
      " (n_dups=", nrow(dup), ")"
    )
  }
  
  # ----------------------------
  # Merge: sequences (many rows) + participants (one row)
  # ----------------------------
  master <- merge(seq, par, by = "pid", all.x = TRUE)
  
  # ----------------------------
  # Basic join diagnostics
  # ----------------------------
  n_missing_par <- master[
    is.na(age) | is.na(sex) | is.na(session) | is.na(treat) | is.na(label),
    uniqueN(pid)
  ]
  if (n_missing_par > 0) warning("Master (BASE): missing participant covariates for ", n_missing_par, " pid(s).")
  
  master[, sex := as.character(sex)]
  setnames(master, "label", "seat")
  
  # ----------------------------
  # Save
  # ----------------------------
  outfile <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  fwrite(master, outfile)
  
  msg("Master sequences (BASE) saved: ", outfile,
      " | rows: ", nrow(master),
      " | unique pid: ", master[, uniqueN(pid)])
  
  invisible(master)
}