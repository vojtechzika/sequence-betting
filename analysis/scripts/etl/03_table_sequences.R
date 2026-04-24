# ============================================================
# 03_table_sequences.R
#
# PURPOSE
#   Extracts sequence-level (trial-level) variables from the merged
#   oTree export, melts from wide to long format, applies controlled
#   coercions, and saves sequences.csv.
#
# INPUT
#   path_src/merged.csv
#
# OUTPUT
#   path_src/sequences.csv
#
# NOTES
#   - One row per participant x sequence trial
#   - Commented variables in vars are intentionally excluded
# ============================================================

make_sequences_table <- function(cfg) {
  
  infile <- file.path(path_src, "merged.csv")
  stopifnot(file.exists(infile))
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # All sequence-level variables available in raw oTree export.
  # Active variables define the schema of sequences.csv.
  # Commented variables are intentionally excluded from the long table.
  vars <- c(
    
    # --- Structural / oTree meta ---
    # "id_in_group",
    # "role",
    # "payoff",
    # "t",
    
    # --- Core sequence task variables ---
    "trial_id",
    "seq",
    "uid",
    "block",
    "pos_in_block",
    "m_used",
    # "treatment",
    "realized",
    "win",
    "round_earnings",
    "side",
    "stake",
    "screen_time_ms",
    "button_order"
    
    # --- Eye-tracking / device diagnostics ---
    # "aoi_boxes_px",
    # "viewport_w_px",
    # "viewport_h_px",
    # "dpr",
    
    # --- Comprehension checks ---
    # "cq_keep_endowment",
    # "cq_multiplier",
    # "cq_toss",
    # "cq_payment",
    # "cq_recency"
  )
  
  # ---- Identify sequence indices present (expect 1..64) ----
  seq_ids <- sort(unique(as.integer(sub(
    "^sequences\\.(\\d+)\\.player\\..*$", "\\1",
    grep("^sequences\\.(\\d+)\\.player\\.", names(dt), value = TRUE)
  ))))
  seq_ids <- seq_ids[!is.na(seq_ids)]
  stopifnot(length(seq_ids) > 0)
  
  # ---- Build measure.vars ONCE, filtered to existing columns ----
  mv <- lapply(vars, function(v) {
    keep <- paste0("sequences.", seq_ids, ".player.", v)
    keep[keep %in% names(dt)]
  })
  names(mv) <- vars
  
  cols_to_melt <- unlist(mv, use.names = FALSE)
  stopifnot(length(cols_to_melt) > 0)
  
  # ---- Keep only pid + the exact sequences.* columns implied by vars ----
  id_cols <- "participant.code"
  id_cols <- id_cols[id_cols %in% names(dt)]
  
  dt_small <- dt[, c(id_cols, cols_to_melt), with = FALSE]
  
  # Safety: ensure nothing else sequences.* slipped in
  other_seq_cols <- setdiff(grep("^sequences\\.", names(dt_small), value = TRUE), cols_to_melt)
  if (length(other_seq_cols) > 0) {
    stop("Unexpected sequences.* columns present: ", paste(other_seq_cols, collapse = ", "))
  }
  
  # ---- Coercion helper ----
  to_num <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "None", "null", "NULL")] <- NA
    suppressWarnings(as.numeric(x))
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
  
  # Drop empty rows
  long <- long[!is.na(uid) & uid != ""]
  
  # Rename participant.code -> pid
  if ("participant.code" %in% names(long)) {
    data.table::setnames(long, "participant.code", "pid")
  }
  
  # ---- Controlled coercion (post-melt) ----
  num_cols <- intersect(
    c("m_used", "win", "round_earnings", "stake", "screen_time_ms", "block", "pos_in_block"),
    names(long)
  )
  for (cc in num_cols) long[, (cc) := to_num(get(cc))]
  
  if ("trial_id" %in% names(long)) long[, trial_id := as.integer(to_num(trial_id))]
  if ("realized" %in% names(long)) {
    long[, realized := toupper(trimws(as.character(realized)))]
    valid_labels <- toupper(names(cfg$design$seq$label_recode))
    long[!realized %in% valid_labels, realized := NA_character_]
  }
  
  # ---- Rename columns ----
  rename_map <- c(
    pos_in_block   = "pos",
    treatment      = "treat",
    round_earnings = "earn",
    screen_time_ms = "screen_ms",
    button_order   = "btn_order",
    m_used         = "multi"
  )
  old <- names(rename_map)
  old <- old[old %in% names(long)]
  data.table::setnames(long, old, rename_map[old])
  
  # ---- Keep final columns derived ONLY from vars (+ pid, seq_no) ----
  final_vars <- vars
  final_vars <- ifelse(final_vars %in% names(rename_map), rename_map[final_vars], final_vars)
  keep_cols  <- c("pid", "seq_no", final_vars)
  keep_cols  <- keep_cols[keep_cols %in% names(long)]
  long       <- long[, ..keep_cols]
  
  outfile <- file.path(path_src, "sequences.csv")
  data.table::fwrite(long, outfile)
  msg("Sequences table saved:", outfile, "| rows:", nrow(long))
  invisible(long)
}