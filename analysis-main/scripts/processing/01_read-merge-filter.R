# ============================================================
# 01_read-merge-filter.R
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
#   path_src/raw-data.csv
#   path_csv/processing_merge-report.csv
#                             -- per raw file: rows read, rows finished,
#                                rows unfinished, participants kept;
#                                last row (TOTAL) sums over all files
#
# DE-DUPLICATION
#   If a participant.code appears more than once (within one file or
#   across files), the script stops with an error listing the duplicates.
#   This forces explicit resolution rather than silently picking a version
#   and guarantees one row per participant for all downstream tables.
#   To use backup files: remove older exports from path_raw before running.
#
# NOTES
#   - Filters to rows where participant._current_page_name == "finished"
#   - Adds the source_file provenance column
#   - Rows without a participant.code are dropped
# ============================================================

read_merge_filter <- function(cfg) {
  
  if (!dir.exists(path_raw)) stop("Raw data folder does not exist: ", path_raw)
  
  files <- list.files(
    path_raw,
    pattern     = "\\.csv$",
    full.names  = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0) stop("No CSV files found in: ", path_raw)
  
  # Row counts per raw file, collected for the sample-flow report
  flow_list <- vector("list", length(files))
  
  dt_list <- lapply(seq_along(files), function(i) {
    f <- files[i]
    x <- data.table::fread(f, encoding = "UTF-8")
    
    if (!"participant._current_page_name" %in% names(x)) {
      stop("Column 'participant._current_page_name' not found in ", basename(f))
    }
    
    n_read <- nrow(x)
    x <- x[tolower(participant._current_page_name) == "finished"]
    
    flow_list[[i]] <<- data.table::data.table(
      source_file   = basename(f),
      rows_read     = n_read,
      rows_finished = nrow(x)
    )
    
    x[, source_file := basename(f)]
    x
  })
  
  dt <- data.table::rbindlist(dt_list, fill = TRUE)
  
  # ---- Participant code ----
  if (!"participant.code" %in% names(dt)) {
    stop("Column 'participant.code' not found in the raw exports.")
  }
  
  dt <- dt[!is.na(participant.code) & trimws(as.character(participant.code)) != ""]
  dt[, participant.code := as.character(participant.code)]
  
  # ---- De-duplication: stop if any participant code appears more than once ----
  dups <- dt[, .N, by = participant.code][N > 1]
  
  if (nrow(dups) > 0) {
    stop(
      "\nThe following participant codes appear more than once ",
      "(within one file or across files):\n",
      paste0("  ", dups$participant.code, collapse = "\n"),
      "\nRemove older export files from ", path_raw,
      " or resolve the duplicates, so that each participant appears in exactly one row.\n"
    )
  }
  
  outfile <- file.path(path_src, "raw-data.csv")
  data.table::fwrite(dt, outfile)
  
  msg("Saved: ", outfile)
  
  # ---- Sample-flow report: rows and participants per raw file ----
  flow <- data.table::rbindlist(flow_list)
  flow[, rows_unfinished := rows_read - rows_finished]
  n_part <- dt[, .(participants = data.table::uniqueN(participant.code)),
               by = source_file]
  flow <- merge(flow, n_part, by = "source_file", all.x = TRUE, sort = FALSE)
  flow[is.na(participants), participants := 0L]
  flow <- rbind(
    flow,
    data.table::data.table(
      source_file     = "TOTAL",
      rows_read       = sum(flow$rows_read),
      rows_finished   = sum(flow$rows_finished),
      rows_unfinished = sum(flow$rows_unfinished),
      participants    = data.table::uniqueN(dt$participant.code)
    ),
    use.names = TRUE
  )
  data.table::setcolorder(
    flow, c("source_file", "rows_read", "rows_finished",
            "rows_unfinished", "participants")
  )
  
  outfile_flow <- file.path(path_csv, "processing_merge-report.csv")
  data.table::fwrite(flow, outfile_flow)
  msg("Saved: ", outfile_flow)
  
  invisible(dt)
}