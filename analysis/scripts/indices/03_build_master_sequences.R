# scripts/analysis/03_build_master_sequences.R
source(here::here("scripts", "00_setup.R"))

library(data.table)

build_master_sequences <- function(dataset = "pilot") {
  
  # ----------------------------
  # Inputs
  # ----------------------------
  f_seq  <- file.path(path_clean_ds(dataset), "sequences.csv")
  f_par  <- file.path(path_clean_ds(dataset), "participants.csv")
  f_lotr <- file.path(path_clean_ds(dataset), "lotr_scored.csv")
  f_mpl  <- file.path(path_clean_ds(dataset), "mpl_scored.csv")
  
  stopifnot(file.exists(f_seq), file.exists(f_par), file.exists(f_lotr), file.exists(f_mpl))
  
  seq  <- fread(f_seq,  encoding = "UTF-8")
  par  <- fread(f_par,  encoding = "UTF-8")
  lotr <- fread(f_lotr, encoding = "UTF-8")
  mpl  <- fread(f_mpl,  encoding = "UTF-8")
  
  stopifnot("pid" %in% names(seq), "pid" %in% names(par), "pid" %in% names(lotr), "pid" %in% names(mpl))
  
  # ----------------------------
  # Keep only the requested columns
  # ----------------------------
  par_keep  <- c("pid", "age", "sex", "session", "treat", "label")
  lotr_keep <- c("pid", "lotr_score", "lotr_z")
  mpl_keep  <- c("pid", "r_mean", "r_median", "inconsistent")
  
  missing_par  <- setdiff(par_keep,  names(par))
  missing_lotr <- setdiff(lotr_keep, names(lotr))
  missing_mpl  <- setdiff(mpl_keep,  names(mpl))
  
  if (length(missing_par)  > 0) stop("participants.csv missing: ", paste(missing_par, collapse = ", "))
  if (length(missing_lotr) > 0) stop("lotr_scored.csv missing: ", paste(missing_lotr, collapse = ", "))
  if (length(missing_mpl)  > 0) stop("mpl_scored.csv missing: ", paste(missing_mpl, collapse = ", "))
  
  par  <- par[,  ..par_keep]
  lotr <- lotr[, ..lotr_keep]
  mpl  <- mpl[,  ..mpl_keep]
  
  # ----------------------------
  # Enforce 1 row per pid in pid-level tables
  # ----------------------------
  chk_unique <- function(dt, name) {
    dup <- dt[, .N, by = pid][N > 1]
    if (nrow(dup) > 0) {
      stop(name, " has duplicate pid rows. Example pids: ",
           paste(head(dup$pid, 5), collapse = ", "),
           " (n_dups=", nrow(dup), ")")
    }
  }
  chk_unique(par,  "participants.csv")
  chk_unique(lotr, "lotr_scored.csv")
  chk_unique(mpl,  "mpl_scored.csv")
  
  # ----------------------------
  # Merge: sequences (many rows) + pid-level covariates (one row)
  # ----------------------------
  master <- merge(seq, par,  by = "pid", all.x = TRUE)
  master <- merge(master, lotr, by = "pid", all.x = TRUE)
  master <- merge(master, mpl,  by = "pid", all.x = TRUE)
  
  # ----------------------------
  # Basic join diagnostics
  # ----------------------------
  n_missing_par  <- master[is.na(age) | is.na(sex) | is.na(session) | is.na(treat) | is.na(label), uniqueN(pid)]
  n_missing_lotr <- master[is.na(lotr_score) | is.na(lotr_z), uniqueN(pid)]
  n_missing_mpl  <- master[is.na(r_mean) | is.na(r_median), uniqueN(pid)]
  
  if (n_missing_par  > 0) warning("Master: missing participant covariates for ", n_missing_par,  " pid(s).")
  if (n_missing_lotr > 0) warning("Master: missing LOT-R scores for ",          n_missing_lotr, " pid(s).")
  if (n_missing_mpl  > 0) warning("Master: missing MPL summaries for ",         n_missing_mpl,  " pid(s).")
  
  # Keep sex as-is (can include "other"). Optionally ensure it is character.
  master[, sex := as.character(sex)]
  
  # Rename some columns for clarity and consistency
  data.table::setnames(master, "label", "seat")
  data.table::setnames(master, "inconsistent", "r_is_inconsistent")
  
  # ----------------------------
  # Save
  # ----------------------------
  outfile <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  fwrite(master, outfile)
  
  msg("Master sequences saved:", outfile, "| rows:", nrow(master), "| unique pid:", master[, uniqueN(pid)])
  invisible(master)
}

# Example:
# build_master_sequences("pilot")