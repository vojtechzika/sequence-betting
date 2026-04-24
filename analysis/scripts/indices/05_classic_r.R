# ============================================================
# 05_classic_r.R
#
# PURPOSE
#   Classical Holt-Laury switching-point estimation of CRRA r_i
#   following Holt and Laury (2002). For each consistent participant
#   finds the switching row and computes the implied r interval
#   and midpoint estimate.
#
# INPUT
#   path_src/mpl.csv
#   path_src/master_sequences.csv  -- for treatment assignment
#   path_src/participants.csv      -- for sex variable in summary
#
# OUTPUT
#   path_out/mpl_switching.csv         -- participant-level results
#   path_out/mpl_switching_summary.csv -- summary statistics by group
#
# NOTES
#   - Inconsistent participants (multiple switches or reversals)
#     are flagged and excluded from r estimates
#   - Indifference condition: EU(A) = EU(B) solved via uniroot
#   - Midpoint of implied r interval used as point estimate
# ============================================================

classic_r <- function(cfg) {
  
  design <- cfg$design
  
  label_safe  <- toupper(trimws(as.character(design$mpl$label_safe)))
  label_risky <- toupper(trimws(as.character(design$mpl$label_risky)))
  
  K_cfg  <- as.integer(design$mpl$K)
  A_high <- as.numeric(design$mpl$A_high)
  A_low  <- as.numeric(design$mpl$A_low)
  B_high <- as.numeric(design$mpl$B_high)
  B_low  <- as.numeric(design$mpl$B_low)
  
  p_sched <- (1:K_cfg) / K_cfg
  
  f_mpl          <- file.path(path_src, "mpl.csv")
  f_master       <- file.path(path_src, "master_sequences.csv")
  f_participants <- file.path(path_src, "participants.csv")
  
  stopifnot(file.exists(f_mpl), file.exists(f_master))
  
  f_out_pid <- file.path(path_out, "mpl_switching.csv")
  f_out_sum <- file.path(path_out, "mpl_switching_summary.csv")
  
  skip_pid <- should_skip(f_out_pid, cfg, "output", "HL switching r (participant)")
  skip_sum <- should_skip(f_out_sum, cfg, "output", "HL switching r (summary)")
  
  if (skip_pid && skip_sum) return(invisible(NULL))
  
  mpl    <- fread(f_mpl,    encoding = "UTF-8")
  master <- fread(f_master, encoding = "UTF-8")
  
  mpl[,    pid   := as.character(pid)]
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  
  pid_treat <- unique(master[, .(pid, treat)])
  mpl <- merge(mpl, pid_treat, by = "pid", all.x = TRUE)
  
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  
  long <- melt(
    mpl,
    id.vars       = c("pid", "treat"),
    measure.vars  = choice_cols,
    variable.name = "row",
    value.name    = "choice"
  )
  
  long[, row    := as.integer(sub("^c", "", row))]
  long[, choice := toupper(trimws(as.character(choice)))]
  long[, y      := fifelse(
    choice == label_safe,  1L,
    fifelse(choice == label_risky, 0L, NA_integer_)
  )]
  long[, p := p_sched[row]]
  
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
    tryCatch(
      uniroot(function(r) eu_diff(r, p), lower = r_lo, upper = r_hi,
              tol = 1e-10)$root,
      error = function(e) NA_real_
    )
  }
  
  indiff_at_row <- sapply(seq_len(K_cfg), function(k) indiff_r(p_sched[k]))
  
  pids <- sort(unique(long$pid))
  
  results <- rbindlist(lapply(pids, function(pid_i) {
    
    d <- long[pid == pid_i & !is.na(y)][order(row)]
    
    if (nrow(d) == 0) {
      return(data.table(pid = pid_i, treat = NA_character_,
                        inconsistent = NA_integer_,
                        switch_row = NA_integer_,
                        r_lo = NA_real_, r_hi = NA_real_,
                        r_hl = NA_real_, note = "no data"))
    }
    
    yy           <- d$y
    n_switch     <- sum(yy[-1] != yy[-length(yy)])
    reversal     <- any(diff(yy) > 0)
    inconsistent <- as.integer(n_switch > 1 || reversal)
    
    if (inconsistent == 1L) {
      return(data.table(
        pid = pid_i, treat = d$treat[1], inconsistent = 1L,
        switch_row = NA_integer_, r_lo = NA_real_,
        r_hi = NA_real_, r_hl = NA_real_, note = "inconsistent"
      ))
    }
    
    switch_pos <- which(yy == 0)
    
    if (length(switch_pos) == 0) {
      return(data.table(
        pid = pid_i, treat = d$treat[1], inconsistent = 0L,
        switch_row = K_cfg + 1L,
        r_lo = indiff_at_row[K_cfg], r_hi = Inf,
        r_hl = indiff_at_row[K_cfg] + 0.5, note = "always safe"
      ))
    }
    
    if (switch_pos[1] == 1L) {
      return(data.table(
        pid = pid_i, treat = d$treat[1], inconsistent = 0L,
        switch_row = 0L,
        r_lo = -Inf, r_hi = indiff_at_row[1],
        r_hl = indiff_at_row[1] - 0.5, note = "always risky"
      ))
    }
    
    k    <- switch_pos[1]
    r_lo <- indiff_at_row[k - 1]
    r_hi <- indiff_at_row[k]
    
    if (is.na(r_lo)) r_lo <- -Inf
    if (is.na(r_hi)) r_hi <-  Inf
    
    r_hl <- if (is.finite(r_lo) && is.finite(r_hi)) (r_lo + r_hi) / 2 else NA_real_
    
    data.table(
      pid = pid_i, treat = d$treat[1], inconsistent = 0L,
      switch_row = as.integer(k),
      r_lo = r_lo, r_hi = r_hi, r_hl = r_hl, note = "ok"
    )
  }))
  
  if (file.exists(f_participants)) {
    participants <- fread(f_participants, encoding = "UTF-8")
    participants[, pid := as.character(pid)]
    participants[, sex_label := fcase(
      toupper(trimws(sex)) == "F", "Women",
      toupper(trimws(sex)) == "M", "Men",
      default = "Other"
    )]
    results <- merge(results, participants[, .(pid, sex_label)],
                     by = "pid", all.x = TRUE)
  } else {
    results[, sex_label := NA_character_]
  }
  
  if (!skip_pid) {
    fwrite(results, f_out_pid)
    msg("Saved: ", f_out_pid)
  }
  
  summarise_r <- function(d, group_label) {
    d <- d[inconsistent == 0L & !is.na(r_hl) & is.finite(r_hl) & note == "ok"]
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
  
  if (!skip_sum) {
    summary_groups <- rbindlist(list(
      summarise_r(results,                                        "pooled"),
      summarise_r(results[treat == "m25"],                        "FN"),
      summarise_r(results[treat == "m19"],                        "FP"),
      summarise_r(results[sex_label == "Women"],                  "pooled / Women"),
      summarise_r(results[sex_label == "Men"],                    "pooled / Men"),
      summarise_r(results[treat == "m25" & sex_label == "Women"], "FN / Women"),
      summarise_r(results[treat == "m25" & sex_label == "Men"],   "FN / Men"),
      summarise_r(results[treat == "m19" & sex_label == "Women"], "FP / Women"),
      summarise_r(results[treat == "m19" & sex_label == "Men"],   "FP / Men")
    ))
    fwrite(summary_groups, f_out_sum)
    msg("Saved: ", f_out_sum)
  }
  
  msg("HL switching-point r estimation complete.")
  msg("  Total:        ", nrow(results),                        " participants")
  msg("  Consistent:   ", nrow(results[inconsistent == 0L]),    " participants")
  msg("  Inconsistent: ", nrow(results[inconsistent == 1L]),    " participants")
  msg("  Always safe:  ", nrow(results[note == "always safe"]), " participants")
  msg("  Always risky: ", nrow(results[note == "always risky"]), " participants")
  
  invisible(results)
}