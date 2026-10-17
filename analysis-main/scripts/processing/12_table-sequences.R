# ============================================================
# 12_table-sequences.R
#
# PURPOSE
#   Extracts sequence-level (trial-level) variables from the raw
#   oTree export, melts from wide to long format, applies controlled
#   coercions, and saves raw-sequences.csv.
#
# INPUT
#   path_src/raw-data.csv
#
# OUTPUT
#   path_src/raw-sequences.csv
#
# NOTES
#   - One row per participant x sequence trial
#   - treat is the treatment recorded by oTree for each trial; it is
#     expected to agree with the participant-level treatment
#   - Stops if a non-empty value of a numeric column cannot be read as a
#     number, or if realized contains a label not defined in
#     cfg$design$seq$label_recode (no value is changed silently)
# ============================================================

table_sequences <- function(cfg) {
  
  infile <- file.path(path_src, "raw-data.csv")
  stopifnot(file.exists(infile))
  dt <- data.table::fread(infile, encoding = "UTF-8")
  stopifnot("participant.code" %in% names(dt))
  
  # Sequence-level variables of the raw export (right) and the short pipeline
  # names used in raw-sequences.csv (left). The order defines the column order of
  # the output.
  #
  # Variables of the export intentionally not carried over:
  #   oTree meta data:       id_in_group, role, payoff, t
  #   eye-tracking / device: aoi_boxes_px, viewport_w_px, viewport_h_px, dpr
  #   comprehension checks:  cq_keep_endowment, cq_multiplier, cq_toss,
  #                          cq_payment, cq_recency
  vars <- c(
    trial_id  = "trial_id",
    seq       = "seq",
    uid       = "uid",
    block     = "block",
    pos       = "pos_in_block",
    multi     = "m_used",
    treat     = "treatment",
    realized  = "realized",
    win       = "win",
    earn      = "round_earnings",
    side      = "side",
    stake     = "stake",
    screen_ms = "screen_time_ms",
    btn_order = "button_order"
  )
  
  # ---- Identify sequence indices present (expect 1..64) ----
  seq_ids <- sort(unique(as.integer(sub(
    "^sequences\\.(\\d+)\\.player\\..*$", "\\1",
    grep("^sequences\\.(\\d+)\\.player\\.", names(dt), value = TRUE)
  ))))
  seq_ids <- seq_ids[!is.na(seq_ids)]
  stopifnot(length(seq_ids) > 0)
  
  # ---- Build measure.vars ONCE, filtered to existing columns ----
  mv <- lapply(unname(vars), function(v) {
    keep <- paste0("sequences.", seq_ids, ".player.", v)
    keep[keep %in% names(dt)]
  })
  names(mv) <- unname(vars)
  
  cols_to_melt <- unlist(mv, use.names = FALSE)
  stopifnot(length(cols_to_melt) > 0)
  
  # ---- Keep only pid + the exact sequences.* columns implied by vars ----
  id_cols <- "participant.code"
  dt_small <- dt[, c(id_cols, cols_to_melt), with = FALSE]
  
  # ---- Coercion helpers ----
  null_strings <- c("", "NA", "NaN", "None", "null", "NULL")
  
  to_num <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% null_strings] <- NA
    suppressWarnings(as.numeric(x))
  }
  
  # Stops if a non-empty value cannot be read as a number
  check_numeric <- function(x, name) {
    x <- trimws(as.character(x))
    x[x %in% null_strings] <- NA
    bad <- unique(x[!is.na(x) & is.na(suppressWarnings(as.numeric(x)))])
    if (length(bad) > 0) {
      stop("Column '", name, "' has values that are not numeric: ",
           paste(head(bad, 10), collapse = ", "))
    }
  }
  
  # Normalize melt cols to character to avoid melt type warnings
  dt_small[, (cols_to_melt) := lapply(.SD, as.character), .SDcols = cols_to_melt]
  
  # ---- Melt wide -> long (ONLY vars) ----
  long <- data.table::melt(
    dt_small,
    id.vars      = id_cols,
    measure.vars = mv,
    variable.name = "seq_no"
  )
  
  # Map seq_no (1..K) to actual index
  long[, seq_no := seq_ids[seq_no]]
  
  # Rename to pipeline names (participant.code -> pid; raw -> short names)
  data.table::setnames(long, c("participant.code", unname(vars)), c("pid", names(vars)))
  
  # Drop empty rows
  long <- long[!is.na(uid) & uid != ""]
  
  # ---- Controlled coercion (post-melt) ----
  num_cols <- intersect(
    c("multi", "win", "earn", "stake", "screen_ms", "block", "pos"),
    names(long)
  )
  for (cc in c(num_cols, intersect("trial_id", names(long)))) {
    check_numeric(long[[cc]], cc)
  }
  for (cc in num_cols) long[, (cc) := to_num(get(cc))]
  
  if ("trial_id" %in% names(long)) long[, trial_id := as.integer(to_num(trial_id))]
  if ("realized" %in% names(long)) {
    long[, realized := toupper(trimws(as.character(realized)))]
    valid_labels <- toupper(names(cfg$design$seq$label_recode))
    unknown <- setdiff(unique(long$realized[!is.na(long$realized) & long$realized != ""]),
                       valid_labels)
    if (length(unknown) > 0) {
      stop("Column 'realized' has labels not in label_recode: ",
           paste(unknown, collapse = ", "),
           ". Expected: ", paste(valid_labels, collapse = ", "))
    }
    long[!realized %in% valid_labels, realized := NA_character_]  # empty -> NA
  }
  
  # ---- Keep final columns (pid, seq_no, then vars in the defined order) ----
  long <- long[, c("pid", "seq_no", names(vars)), with = FALSE]
  
  outfile <- file.path(path_src, "raw-sequences.csv")
  data.table::fwrite(long, outfile)
  msg("Saved: ", outfile)
  invisible(long)
}