# ============================================================
# 01_read_merge_filter.R
#
# PURPOSE
#   Reads all oTree CSV exports from the raw data folder, filters to
#   completed participants (page = "finished"), de-duplicates across
#   backup files, and saves a single merged CSV.
#
# INPUT
#   path_raw/          -- one or more oTree CSV exports (.csv)
#
# OUTPUT
#   path_src/merged.csv
#
# DE-DUPLICATION
#   If a participant.code appears in more than one source file, the
#   script stops with an error listing the duplicates. This forces
#   explicit resolution rather than silently picking a version.
#   To use backup files: remove older exports from path_raw before running.
#
# NOTES
#   - Filters to rows where participant._current_page_name == "finished"
#   - Adds source_file and source_mtime provenance columns
#   - Handles variant participant code column names
# ============================================================

run_read_merge_filter <- function(cfg) {
  
  if (!dir.exists(path_raw)) stop("Raw data folder does not exist: ", path_raw)
  
  files <- list.files(
    path_raw,
    pattern     = "\\.csv$",
    full.names  = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0) stop("No CSV files found in: ", path_raw)
  
  msg("Data folder:", cfg$run$data_folder)
  msg("Raw dir:", path_raw)
  msg("Files detected:", length(files))
  
  total_rows_before <- 0L
  total_rows_after  <- 0L
  
  dt_list <- lapply(files, function(f) {
    x <- data.table::fread(f, encoding = "UTF-8")
    rows_before <- nrow(x)
    total_rows_before <<- total_rows_before + rows_before
    
    if (!"participant._current_page_name" %in% names(x)) {
      stop("Column 'participant._current_page_name' not found in ", basename(f))
    }
    
    x <- x[tolower(participant._current_page_name) == "finished"]
    rows_after <- nrow(x)
    total_rows_after <<- total_rows_after + rows_after
    
    msg("  File:", basename(f),
        "| before:", rows_before,
        "| kept:", rows_after,
        "| dropped:", rows_before - rows_after)
    
    x[, source_file  := basename(f)]
    x[, source_mtime := as.POSIXct(file.info(f)$mtime)]
    x
  })
  
  dt <- data.table::rbindlist(dt_list, fill = TRUE)
  
  # ---- Resolve participant code column ----
  code_candidates <- c("participant.code", "participant_code",
                       "participant_code_", "participantcode")
  code_col <- code_candidates[code_candidates %in% names(dt)][1]
  
  if (is.na(code_col) || !nzchar(code_col)) {
    stop(
      "Could not find participant code column. Tried: ",
      paste(code_candidates, collapse = ", "),
      "\nAvailable columns include: ",
      paste(head(names(dt), 50), collapse = ", ")
    )
  }
  
  if (code_col != "participant.code") {
    data.table::setnames(dt, code_col, "participant.code")
  }
  
  dt <- dt[!is.na(participant.code) & trimws(as.character(participant.code)) != ""]
  dt[, participant.code := as.character(participant.code)]
  
  # ---- De-duplication: stop if any participant appears in more than one file ----
  dup_check <- dt[, .(n_files = uniqueN(source_file)), by = participant.code]
  dups      <- dup_check[n_files > 1]
  
  if (nrow(dups) > 0) {
    stop(
      "\nThe following participant codes appear in more than one source file:\n",
      paste0("  ", dups$participant.code, collapse = "\n"),
      "\nRemove older export files from ", path_raw,
      " so that each participant appears in exactly one file.\n"
    )
  }
  
  n_unique <- uniqueN(dt$participant.code)
  
  outfile <- file.path(path_src, "merged.csv")
  data.table::fwrite(dt, outfile)
  
  msg("----------------------------------------")
  msg("Total rows before filter:", total_rows_before)
  msg("Total rows kept (finished):", total_rows_after)
  msg("Unique participants:", n_unique)
  msg("Saved to:", outfile)
  
  invisible(dt)
}