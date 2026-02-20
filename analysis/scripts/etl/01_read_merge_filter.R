source(here::here("scripts", "00_setup.R"))

run_read_merge_filter <- function(dataset = "pilot") {
  
  # --------------------------------------------
  # Use dataset-specific raw folder
  # --------------------------------------------
  raw_dir <- path_raw_ds(dataset)
  if (!dir.exists(raw_dir)) stop("Raw dataset folder does not exist: ", raw_dir)
  
  files <- list.files(
    raw_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  stopifnot(length(files) > 0)
  
  msg("Dataset:", dataset)
  msg("Raw dir:", raw_dir)
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
    
    # provenance
    x[, source_file := basename(f)]
    x[, source_mtime := as.POSIXct(file.info(f)$mtime)]
    
    x
  })
  
  dt <- data.table::rbindlist(dt_list, fill = TRUE)
  
  # --------------------------------------------
  # De-duplicate by participant.code across backups
  # Keep the newest row per participant.code (by file mtime)
  # --------------------------------------------
  code_candidates <- c("participant.code", "participant_code", "participant_code_", "participantcode")
  code_col <- code_candidates[code_candidates %in% names(dt)][1]
  
  if (is.na(code_col) || !nzchar(code_col)) {
    stop(
      "Could not find participant code column. Tried: ",
      paste(code_candidates, collapse = ", "),
      "\nAvailable columns include: ",
      paste(head(names(dt), 50), collapse = ", ")
    )
  }
  
  # Normalize to a single column name for internal handling
  if (code_col != "participant.code") {
    data.table::setnames(dt, code_col, "participant.code")
    code_col <- "participant.code"
  }
  
  # Drop empty codes (should not happen, but protects the dedup logic)
  dt <- dt[!is.na(get(code_col)) & trimws(as.character(get(code_col))) != ""]
  dt[, (code_col) := as.character(get(code_col))]
  
  # Keep newest file version per participant.code
  data.table::setorder(dt, participant.code, source_mtime)
  n_before_dedup <- uniqueN(dt$participant.code)
  
  dt_unique <- dt[, .SD[.N], by = participant.code]  # last row within each code after sorting
  n_after_dedup <- nrow(dt_unique)
  
  msg("----------------------------------------")
  msg("Unique participant.code before dedup (from merged finished rows):", n_before_dedup)
  msg("Unique participant.code after dedup:", n_after_dedup)
  msg("Duplicates removed:", nrow(dt) - nrow(dt_unique))
  
  outfile <- file.path(path_clean_ds(dataset), "merged.csv")
  data.table::fwrite(dt_unique, outfile)
  
  msg("----------------------------------------")
  msg("Total rows before filter:", total_rows_before)
  msg("Total rows kept (finished, across all files):", total_rows_after)
  msg("Final merged unique rows:", nrow(dt_unique))
  msg("Saved to:", outfile)
  
  invisible(dt_unique)
}