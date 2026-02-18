# scripts/03_table_sequences.R
source(here::here("scripts", "00_setup.R"))

make_sequences_table <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "merged.csv")
  stopifnot(file.exists(infile))
  dt <- data.table::fread(infile, encoding = "UTF-8")
  
  # All sequence-level variables available in raw oTree export.
  # Active variables define the schema of sequences.csv.
  # Commented variables are intentionally excluded from the long table.
  vars <- c(
    
    # --- Structural / oTree meta ---
    # "id_in_group",     # oTree internal group position (not analytically relevant in 1-player task)
    # "role",            # oTree role label (unused in single-role design)
    # "payoff",          # cumulative or round payoff assigned by oTree (redundant if round_earnings kept)
    # "t",               # internal trial counter (may duplicate trial_id or pos_in_block)
    
    # --- Core sequence task variables ---
    "trial_id",          # unique trial identifier within participant
    "seq",               # 6-toss sequence shown (e.g., HHTOHO)
    "uid",               # unique sequence identifier (stable across participants)
    "block",             # block number
    "pos_in_block",      # position within block
    "m_used",            # multiplier used in this trial
    # "treatment",       # treatment condition (sequence-level; constant per participant)
    "realized",          # actual realized outcome (H or O)
    "win",               # indicator whether participant won (1/0 or TRUE/FALSE)
    "round_earnings",    # earnings from this trial
    "side",              # side (H or O) a participant chose to bet on
    "stake",             # amount wagered
    "screen_time_ms",    # decision screen time in milliseconds
    "button_order"      # randomized order of response buttons (HO / OH)
    
    # --- Eye-tracking / device diagnostics ---
    # "aoi_boxes_px",    # AOI bounding boxes (pixel coordinates)
    # "viewport_w_px",   # browser viewport width (px)
    # "viewport_h_px",   # browser viewport height (px)
    # "dpr",             # device pixel ratio
    
    # --- Comprehension checks ---
    # --- Only first sequence has comprehension checks, so these will be mostly NA in the long format, plus participant MUST complete checks to continue ---
    # "cq_keep_endowment",  # comprehension: keep endowment question
    # "cq_multiplier",      # comprehension: multiplier question
    # "cq_toss",            # comprehension: coin toss probability question
    # "cq_payment",         # comprehension: payment rule question
    # "cq_recency"          # comprehension: recency understanding question
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
  
  # ---- Coercion helpers ----
  to_num <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "None", "null", "NULL")] <- NA
    suppressWarnings(as.numeric(x))
  }
  
  to_bool01 <- function(x) {
    x <- trimws(tolower(as.character(x)))
    x[x %in% c("", "na", "nan", "none", "null", "NULL")] <- NA
    data.table::fifelse(
      x %in% c("1", "true", "t", "yes", "y"),
      1,
      data.table::fifelse(
        x %in% c("0", "false", "f", "no", "n"),
        0,
        suppressWarnings(as.numeric(x))
      )
    )
  }
  
  # Normalize melt cols to character to avoid melt type warnings
  dt_small[, (cols_to_melt) := lapply(.SD, as.character), .SDcols = cols_to_melt]
  
  # ---- Melt wide -> long (ONLY vars) ----
  long <- data.table::melt(
    dt_small,
    id.vars = id_cols,
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
    long[!realized %in% c("H", "O"), realized := NA_character_]
  }
  
  # ---- Rename columns (this does NOT define what we keep) ----
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
  # Translate vars through rename_map where applicable
  final_vars <- vars
  final_vars <- ifelse(final_vars %in% names(rename_map), rename_map[final_vars], final_vars)
  keep_cols <- c("pid", "seq_no", final_vars)
  keep_cols <- keep_cols[keep_cols %in% names(long)]
  long <- long[, ..keep_cols]
  
  # ---- Save ----
  outfile <- file.path(path_clean_ds(dataset), "sequences.csv")
  data.table::fwrite(long, outfile)
  msg("Sequences table saved:", outfile, "| rows:", nrow(long))
  invisible(long)
}