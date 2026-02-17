source("scripts/00_setup.R")


run_read_merge_filter <- function(dataset) {
  
  pattern <- if (dataset == "pilot") {
    "pilot"
  } else {
    "main"
  }
  
  files <- list.files(
    path_raw,
    pattern = paste0(pattern, ".*\\.csv$"),
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  stopifnot(length(files) > 0)
  
  msg("Dataset:", dataset)
  msg("Files detected:", length(files))
  
  total_rows_before <- 0
  total_rows_after  <- 0
  
  dt_list <- lapply(files, function(f) {
    
    x <- data.table::fread(f)
    
    rows_before <- nrow(x)
    total_rows_before <<- total_rows_before + rows_before
    
    if (!"participant._current_page_name" %in% names(x)) {
      stop("Column 'participant._current_page_name' not found in ", basename(f))
    }
    
    x <- x[
      tolower(participant._current_page_name) == "finished"
    ]
    
    rows_after <- nrow(x)
    total_rows_after <<- total_rows_after + rows_after
    
    msg("  File:", basename(f),
        "| before:", rows_before,
        "| kept:", rows_after,
        "| dropped:", rows_before - rows_after)
    
    x[, source_file := basename(f)]
    x
  })
  
  dt <- data.table::rbindlist(dt_list, fill = TRUE)
  
  outfile <- file.path(path_clean_ds(dataset), paste0("merged.csv"))
  data.table::fwrite(dt, outfile)
  
  msg("----------------------------------------")
  msg("Total rows before filter:", total_rows_before)
  msg("Total rows kept:", total_rows_after)
  msg("Total rows dropped:", total_rows_before - total_rows_after)
  msg("Final merged rows:", nrow(dt))
  msg("Saved to:", outfile)
}