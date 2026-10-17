# ============================================================
# 19_master-sequences.R
#
# PURPOSE
#   Builds the trial-level analysis table. Takes the sequence data
#   extracted from oTree, validates it against the participant table
#   and the design configuration, recodes oTree coin labels to the
#   labels used in the analysis, and adds deterministic sequence-level
#   rules.
#
# INPUT
#   path_src/raw-sequences.csv     -- trial-level data with oTree labels
#   path_src/participants.csv      -- used only to validate treatment
#
# OUTPUT
#   path_der/master-sequences.csv  -- one row per participant x sequence
#
# NOTES
#   - Label recode (cfg$design$seq$label_recode), e.g. O -> T, is
#     applied to seq, realized, side and btn_order; any label outside
#     the expected set stops the script
#   - Validation: trial treatment equals the participant's treatment,
#     and the trial multiplier equals the multiplier configured for
#     that treatment
#   - Two sequence-level rules are computed from `seq` (the string of
#     prior outcomes), independent of the participant's bet:
#       imbalance        : mean of the tosses coded H = +1, T = -1
#                          (-1..+1; 0 when heads and tails are equal)
#       imbalance_abs    : absolute value of imbalance
#       terminal_run     : sign of the terminal run of length >= 2:
#                          +1 if the sequence ends "HH", -1 if it ends
#                          "TT", 0 otherwise (coded 0 rather than NA so
#                          that it enters linear models without
#                          dropping rows)
#       has_terminal_run : 1 if terminal_run is not 0, else 0
#   - Scope: this table holds only deterministic functions of the trial
#     data and the configuration. Quantities that depend on estimates or
#     model selection are added in later stages
#   - Rows are ordered by pid and seq_no
# ============================================================

master_sequences <- function(cfg) {
  
  f_seq <- file.path(path_src, "raw-sequences.csv")
  f_par <- file.path(path_src, "participants.csv")
  stopifnot(file.exists(f_seq), file.exists(f_par))
  
  trials <- fread(f_seq, encoding = "UTF-8")
  par    <- fread(f_par, encoding = "UTF-8")
  
  need <- c("pid", "seq_no", "seq", "multi", "treat", "realized", "side", "btn_order")
  if (!all(need %in% names(trials))) {
    stop("raw-sequences.csv missing: ", paste(setdiff(need, names(trials)), collapse = ", "))
  }
  stopifnot(all(c("pid", "treat") %in% names(par)))
  
  trials[, pid   := as.character(pid)]
  trials[, treat := as.character(treat)]
  par[,    pid   := as.character(pid)]
  par[,    treat := as.character(treat)]
  
  # --------------------------------------------------------
  # Validation against participants and design
  # --------------------------------------------------------
  if (anyDuplicated(par$pid)) stop("participants.csv has duplicate pid rows.")
  if (!all(trials$pid %in% par$pid)) stop("Trials with pid missing from participants.csv.")
  
  par_treat <- par$treat[match(trials$pid, par$pid)]
  if (any(trials$treat != par_treat)) {
    stop("Trial treatment differs from participants.csv for some pid(s).")
  }
  
  mult_cfg <- unlist(cfg$design$seq$treatments)
  if (!all(trials$treat %in% names(mult_cfg))) stop("Trial treatment not in cfg$design$seq$treatments.")
  if (any(abs(as.numeric(trials$multi) - mult_cfg[trials$treat]) > 1e-9)) {
    stop("Trial multiplier differs from the configured treatment multiplier.")
  }
  
  # --------------------------------------------------------
  # Label recode (oTree -> analysis labels)
  # --------------------------------------------------------
  recode    <- unlist(cfg$design$seq$label_recode)      # names = oTree label
  char_map  <- recode[nchar(names(recode)) == 1L]       # single-character labels
  
  # Whole-label recode (side, realized): every value must be a known label
  recode_whole <- function(x, what) {
    out <- unname(recode[as.character(x)])
    if (anyNA(out)) stop("Unknown label in ", what, ".")
    out
  }
  # Character-wise recode (seq, btn_order): every character must be known
  recode_chars <- function(x, what) {
    u   <- unique(as.character(x))
    new <- vapply(u, function(s) {
      ch <- strsplit(s, "")[[1]]
      m  <- unname(char_map[ch])
      if (anyNA(m)) stop("Unknown character in ", what, ": ", s)
      paste(m, collapse = "")
    }, character(1))
    unname(new[match(as.character(x), u)])
  }
  
  trials[, seq       := recode_chars(seq,       "seq")]
  trials[, btn_order := recode_chars(btn_order, "btn_order")]
  trials[, realized  := recode_whole(realized,  "realized")]
  trials[, side      := recode_whole(side,      "side")]
  
  # --------------------------------------------------------
  # Sequence-level rules (computed once per unique seq)
  # --------------------------------------------------------
  heads <- as.character(cfg$design$seq$side_labels$heads)
  tails <- as.character(cfg$design$seq$side_labels$tails)
  
  rules <- function(s) {
    toss <- strsplit(s, "")[[1]]
    if (length(toss) < 2L || !all(toss %in% c(heads, tails))) {
      stop("Sequence not a string of ", heads, "/", tails, " tosses: ", s)
    }
    x     <- ifelse(toss == heads, 1L, -1L)
    last2 <- tail(x, 2L)
    run   <- if (last2[1] == last2[2]) last2[2] else 0L
    list(seq = s, imbalance = mean(x), imbalance_abs = abs(mean(x)),
         terminal_run = run, has_terminal_run = as.integer(run != 0L))
  }
  lookup <- rbindlist(lapply(unique(trials$seq), rules))
  
  # Join by position so that trial order is preserved
  idx <- match(trials$seq, lookup$seq)
  trials[, c("imbalance", "imbalance_abs", "terminal_run", "has_terminal_run") :=
           lookup[idx, .(imbalance, imbalance_abs, terminal_run, has_terminal_run)]]
  
  setorder(trials, pid, seq_no)
  
  outfile <- file.path(path_der, "master-sequences.csv")
  fwrite(trials, outfile)
  msg("Saved: ", outfile)
  
  invisible(trials)
}
