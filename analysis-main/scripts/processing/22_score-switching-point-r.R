# ============================================================
# 22_score-switching-point-r.R
#
# PURPOSE
#   Classical Holt-Laury switching-point estimation of CRRA r_i
#   following Holt and Laury (2002). For each consistent participant
#   finds the switching row and computes the implied r interval
#   and midpoint estimate.
#
# INPUT
#   path_src/holt-laury.csv        -- MPL choices c1..cK (A/B)
#   path_src/participants.csv      -- treatment assignment and sex
#
# OUTPUT
#   path_der/switching-point-r.csv  -- participant-level results
#   path_csv/processing_switching-point-r-summary.csv -- summary by group
#
# NOTES
#   - Inconsistent participants (more than one switch, or a switch from
#     the risky back to the safe lottery) are flagged and excluded
#   - Indifference condition EU(A) = EU(B) is solved with uniroot
#   - Midpoint of the implied r interval is the point estimate
#   - Participants who always choose the safe (risky) option have an
#     open-ended interval; their point estimate is the boundary
#     shifted by +0.5 (-0.5) and they are excluded from the summary
#   - Summary groups are the treatments in cfg$design$seq$treatments,
#     optionally split by sex
#   - Unrecognised choices, unknown treatments or sex codes, and
#     failed root finding stop the script
# ============================================================

score_switching_point_r <- function(cfg) {
  
  design <- cfg$design
  
  label_safe  <- toupper(trimws(as.character(design$mpl$label_safe)))
  label_risky <- toupper(trimws(as.character(design$mpl$label_risky)))
  
  K_cfg  <- as.integer(design$mpl$K)
  A_high <- as.numeric(design$mpl$A_high)
  A_low  <- as.numeric(design$mpl$A_low)
  B_high <- as.numeric(design$mpl$B_high)
  B_low  <- as.numeric(design$mpl$B_low)
  
  p_sched <- as.numeric(design$mpl$p)
  stopifnot(length(p_sched) == K_cfg, all(p_sched > 0), all(p_sched < 1))
  
  treatments <- names(design$seq$treatments)
  
  f_mpl          <- file.path(path_src, "holt-laury.csv")
  f_participants <- file.path(path_src, "participants.csv")
  stopifnot(file.exists(f_mpl), file.exists(f_participants))
  
  f_out_pid <- file.path(path_der, "switching-point-r.csv")
  f_out_sum <- file.path(path_csv, "processing_switching-point-r-summary.csv")
  
  mpl          <- fread(f_mpl,          encoding = "UTF-8")
  participants <- fread(f_participants, encoding = "UTF-8")
  
  mpl[,          pid   := as.character(pid)]
  participants[, pid   := as.character(pid)]
  participants[, treat := as.character(treat)]
  
  # Treatment and sex come from the participant table (one row per pid)
  if (anyDuplicated(participants$pid)) stop("Repeated pid in participants.csv.")
  if (!all(mpl$pid %in% participants$pid)) {
    stop("MPL participants missing from participants.csv.")
  }
  if (!all(participants[pid %in% mpl$pid, treat] %in% treatments)) {
    stop("Unknown treatment in participants.csv.")
  }
  
  sex_map <- c(F = "Women", M = "Men", O = "Other")
  participants[, sex_label := unname(sex_map[toupper(trimws(sex))])]
  if (anyNA(participants[pid %in% mpl$pid, sex_label])) {
    stop("Unknown sex code in participants.csv.")
  }
  
  mpl <- merge(mpl, participants[, .(pid, treat, sex_label)], by = "pid")
  
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  stopifnot(length(choice_cols) == K_cfg)
  
  # Long format: one row per participant and MPL row
  long <- melt(
    mpl,
    id.vars       = c("pid", "treat", "sex_label"),
    measure.vars  = choice_cols,
    variable.name = "row",
    value.name    = "choice"
  )
  
  long[, row    := as.integer(sub("^c", "", row))]
  long[, choice := toupper(trimws(as.character(choice)))]
  if (!all(long$choice %in% c(label_safe, label_risky))) {
    stop("MPL choices contain labels other than the configured safe/risky labels.")
  }
  long[, y := fifelse(choice == label_safe, 1L, 0L)]
  
  # CRRA utility and the indifference condition between lotteries A and B
  crra_u <- function(x, r) {
    if (abs(r - 1) < 1e-10) return(log(x))
    x^(1 - r) / (1 - r)
  }
  
  eu_diff <- function(r, p) {
    eu_a <- p * crra_u(A_high, r) + (1 - p) * crra_u(A_low, r)
    eu_b <- p * crra_u(B_high, r) + (1 - p) * crra_u(B_low, r)
    eu_a - eu_b
  }
  
  indiff_r <- function(p, r_lo = -10, r_hi = 10) {
    uniroot(function(r) eu_diff(r, p), lower = r_lo, upper = r_hi,
            tol = 1e-10)$root
  }
  
  # r at which a participant is indifferent in each MPL row
  indiff_at_row <- sapply(seq_len(K_cfg), function(k) indiff_r(p_sched[k]))
  
  pids <- sort(unique(long$pid))
  
  results <- rbindlist(lapply(pids, function(pid_i) {
    
    d  <- long[pid == pid_i][order(row)]
    yy <- d$y
    
    # Safe choices (1) must precede risky choices (0): inconsistent if there
    # is more than one switch or any switch from risky back to safe
    n_switch     <- sum(yy[-1] != yy[-length(yy)])
    reversal     <- any(diff(yy) > 0)
    inconsistent <- as.integer(n_switch > 1 || reversal)
    
    base <- list(pid = pid_i, treat = d$treat[1], sex_label = d$sex_label[1])
    
    if (inconsistent == 1L) {
      return(as.data.table(c(base, list(
        inconsistent = 1L, switch_row = NA_integer_,
        r_lo = NA_real_, r_hi = NA_real_, r_hl = NA_real_,
        note = "inconsistent"
      ))))
    }
    
    switch_pos <- which(yy == 0)
    
    if (length(switch_pos) == 0) {
      return(as.data.table(c(base, list(
        inconsistent = 0L, switch_row = K_cfg + 1L,
        r_lo = indiff_at_row[K_cfg], r_hi = Inf,
        r_hl = indiff_at_row[K_cfg] + 0.5, note = "always safe"
      ))))
    }
    
    if (switch_pos[1] == 1L) {
      return(as.data.table(c(base, list(
        inconsistent = 0L, switch_row = 0L,
        r_lo = -Inf, r_hi = indiff_at_row[1],
        r_hl = indiff_at_row[1] - 0.5, note = "always risky"
      ))))
    }
    
    k    <- switch_pos[1]
    r_lo <- indiff_at_row[k - 1]
    r_hi <- indiff_at_row[k]
    
    as.data.table(c(base, list(
      inconsistent = 0L, switch_row = as.integer(k),
      r_lo = r_lo, r_hi = r_hi, r_hl = (r_lo + r_hi) / 2, note = "ok"
    )))
  }))
  
  setcolorder(results, c("pid", "treat", "inconsistent", "switch_row",
                         "r_lo", "r_hi", "r_hl", "note", "sex_label"))
  
  fwrite(results, f_out_pid)
  msg("Saved: ", f_out_pid)
  
  # Summary statistics of the point estimate among interval-identified participants
  summarise_r <- function(d, group_label) {
    d <- d[note == "ok"]
    if (nrow(d) == 0) return(data.table(
      group = group_label, n = 0L,
      r_mean = NA_real_, r_sd = NA_real_, r_median = NA_real_,
      r_q25 = NA_real_,  r_q75 = NA_real_,
      r_min = NA_real_,  r_max = NA_real_
    ))
    data.table(
      group    = group_label,
      n        = nrow(d),
      r_mean   = mean(d$r_hl),
      r_sd     = sd(d$r_hl),
      r_median = median(d$r_hl),
      r_q25    = quantile(d$r_hl, 0.25),
      r_q75    = quantile(d$r_hl, 0.75),
      r_min    = min(d$r_hl),
      r_max    = max(d$r_hl)
    )
  }
  
  # Pooled, by sex, by treatment, and by treatment x sex
  sexes  <- c("Women", "Men")
  groups <- list(summarise_r(results, "pooled"))
  for (s in sexes) {
    groups[[length(groups) + 1L]] <-
      summarise_r(results[sex_label == s], paste("pooled /", s))
  }
  for (tr in treatments) {
    groups[[length(groups) + 1L]] <- summarise_r(results[treat == tr], tr)
    for (s in sexes) {
      groups[[length(groups) + 1L]] <-
        summarise_r(results[treat == tr & sex_label == s], paste(tr, "/", s))
    }
  }
  
  fwrite(rbindlist(groups), f_out_sum)
  msg("Saved: ", f_out_sum)
  
  invisible(results)
}
