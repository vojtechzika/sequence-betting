# ============================================================
# scripts/analysis/52_ex1_1_sequences.R
#
# EX1.I (part 2): Sequence-level directional tendencies (chi_s)
#   PREREG-LITERAL implementation.
#
# Implements prereg:
# 1) Uses anchor-based similarity *probabilities* d_a(s) computed in 51_ex1_1_anchors.R:
#      d_a(s) = P(|theta_s - theta_a| < eps | data)
#    and the corresponding *scalar* weights:
#      w_a(s) = (d_a(s)+eta) / (dH(s)+dT(s)+d0(s)+3*eta)
#
# 2) Directional score per betting trial (i,s) with side choice h_is in {0,1}:
#      z_is = P(hot|h_is,s) - P(gambler|h_is,s)  in [-1,1]
#    Under the prereg weight rule, this simplifies to:
#      z_is = (2*h_is - 1) * ( wH(s) - wT(s) )
#    because the w0(s)*0.5 term cancels in the difference.
#
# 3) Sequence-level directional tendency:
#      chi_s = E_i[z_is]
#    Operationalized here as the mean of z_is over betting trials on sequence s.
#
# IMPORTANT (prereg consequence):
# - Because wH(s), wT(s), w0(s) are defined as posterior *probabilities* (scalars),
#   chi_s is a plug-in posterior functional (no per-draw chi_s^(k) is constructed here).
#   Therefore the prereg probabilities H(.)=P(chi_s>delta), T(.)=P(chi_s<-delta),
#   and P(|chi_s|>delta) collapse to 0/1 indicators under this literal implementation.
#
# Outputs (per treatment):
#   data/clean/<ds>/output/ex1_1_sequences_<tr>.csv
#
# Requires upstream:
#   - RQ4 fit artifacts
#   - 51_ex1_1_anchors.R producing: models/ex1_1_anchors_<tr>.rds
# ============================================================

library(data.table)

ex1_1_sequences <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds     <- as.character(cfg$run$dataset)
  design <- cfg$design
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Required design labels
  # ----------------------------
  stopifnot(!is.null(design$seq), !is.null(design$seq$side_labels))
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads))
  stopifnot(length(lab_tails) == 1L, nzchar(lab_tails))
  stopifnot(lab_heads != lab_tails)
  
  # delta grid: main + sensitivities (prereg default)
  delta_vec <- c(0.05, 0.03, 0.08)
  if (!is.null(design$rhos) && !is.null(design$rhos$ex1_delta)) {
    delta_vec <- as.numeric(design$rhos$ex1_delta)
  }
  stopifnot(length(delta_vec) >= 1L, all(is.finite(delta_vec)), all(delta_vec > 0))
  
  # eps main and grid must match the anchors file
  eps_main <- 0.05
  
  # ----------------------------
  # IO paths
  # ----------------------------
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  # load master once
  dt <- fread(infile, encoding = "UTF-8")
  required <- c("pid", "treat", "seq", "stake", "side")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) stop("master_sequences.csv missing columns: ", paste(missing, collapse = ", "))
  
  dt[, pid := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ----------------------------
    # Load anchor weights (from 51_)
    # ----------------------------
    f_anc <- file.path(mod_dir, paste0("ex1_1_anchors_", tr, ".rds"))
    if (!file.exists(f_anc)) {
      stop("EX1.1 sequences: missing anchors RDS: ", f_anc, "\nRun 51_ex1_1_anchors(cfg) first.")
    }
    anc <- readRDS(f_anc)
    
    if (is.null(anc$seq_levels) || is.null(anc$table)) {
      stop("EX1.1 sequences: anchors RDS must contain seq_levels and table.")
    }
    seq_levels <- as.character(anc$seq_levels)
    tab_anc <- as.data.table(anc$table)
    
    # Sanity: same sequences
    if (!identical(tab_anc$sequence, seq_levels)) {
      # allow reordering if needed
      setkey(tab_anc, sequence)
      tab_anc <- tab_anc[.(seq_levels)]
      if (anyNA(tab_anc$sequence)) stop("EX1.1 sequences: anchors table could not be aligned to seq_levels.")
    }
    
    # Determine eps grid from anchors RDS (preferred)
    eps_vec <- anc$eps_vec
    if (is.null(eps_vec)) eps_vec <- c(0.05, 0.03, 0.08)
    eps_vec <- as.numeric(eps_vec)
    stopifnot(length(eps_vec) >= 1L, all(is.finite(eps_vec)), all(eps_vec > 0))
    
    # ----------------------------
    # Prepare observed data for this treatment (betting trials only)
    # ----------------------------
    d <- dt[treat == tr & is.finite(stake) & stake > 0]
    if (nrow(d) == 0) {
      warning("EX1.1 sequences: no betting trials for treatment='", tr, "'. Skipping.")
      next
    }
    
    # Validate side values among betting trials
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0) {
      stop(
        "EX1.1 sequences: side values outside {heads,tails} after stake>0 (treatment='", tr, "'). ",
        "Expected {'", lab_heads, "','", lab_tails, "'}. Examples: ",
        paste(unique(head(bad$side, 10)), collapse = ", ")
      )
    }
    
    d[, h := as.integer(side == lab_heads)]          # 1=Heads, 0=Tails
    d[, sign := 2L * h - 1L]                         # in {-1,+1}
    d[, sid := match(seq, seq_levels)]
    if (anyNA(d$sid)) {
      stop("EX1.1 sequences: some seq values not found in RQ4 seq_levels for treatment='", tr, "'.")
    }
    
    # Base output table
    tbl <- data.table(sequence = seq_levels)
    
    # ----------------------------
    # For each eps: compute scalar weights and then chi_s
    # ----------------------------
    for (eps in eps_vec) {
      
      suf <- gsub("\\.", "", sprintf("%.2f", eps))
      
      col_wH <- paste0("wH_eps_", suf)
      col_wT <- paste0("wT_eps_", suf)
      col_w0 <- paste0("w0_eps_", suf)
      
      if (!all(c(col_wH, col_wT, col_w0) %in% names(tab_anc))) {
        stop(
          "EX1.1 sequences: anchors table missing weight columns for eps=", eps, " in treatment='", tr, "'.\n",
          "Expected: ", col_wH, ", ", col_wT, ", ", col_w0
        )
      }
      
      wH <- as.numeric(tab_anc[[col_wH]])
      wT <- as.numeric(tab_anc[[col_wT]])
      # w0 exists but cancels in z_is; still keep for completeness if you ever need it:
      # w0 <- as.numeric(tab_anc[[col_w0]])
      
      if (any(!is.finite(wH)) || any(!is.finite(wT))) stop("Non-finite weights in anchors table (tr=", tr, ", eps=", eps, ").")
      
      diff_w <- wH - wT  # length S
      
      # trial-level z_is (prereg literal)
      d[, z := sign * diff_w[sid]]
      
      # chi_s = mean(z_is) over betting trials on s
      chi_by_seq <- d[, .(chi = mean(z, na.rm = TRUE), n_bets = .N), by = sid]
      chi_vec <- rep(NA_real_, length(seq_levels))
      n_vec   <- rep(0L, length(seq_levels))
      chi_vec[chi_by_seq$sid] <- chi_by_seq$chi
      n_vec[chi_by_seq$sid]   <- chi_by_seq$n_bets
      
      tbl[, (paste0("n_bets_eps_", suf)) := n_vec]
      tbl[, (paste0("chi_eps_", suf)) := chi_vec]
      
      # prereg probabilities H(.) and T(.) collapse under scalar chi (0/1)
      for (delta in delta_vec) {
        dsu <- gsub("\\.", "", sprintf("%.2f", delta))
        tbl[, (paste0("H_eps_", suf, "_delta_", dsu)) := as.numeric(chi_vec >  delta)]
        tbl[, (paste0("T_eps_", suf, "_delta_", dsu)) := as.numeric(chi_vec < -delta)]
        tbl[, (paste0("Pabs_eps_", suf, "_delta_", dsu)) := as.numeric(abs(chi_vec) > delta)]
      }
    }
    
    # ----------------------------
    # Direction label (7-point) for MAIN eps + MAIN delta
    # (collapses to deterministic under prereg-literal scalar chi)
    # ----------------------------
    suf_eps_main <- gsub("\\.", "", sprintf("%.2f", eps_main))
    dmain <- delta_vec[1]
    dsu_main <- gsub("\\.", "", sprintf("%.2f", dmain))
    
    colH <- paste0("H_eps_", suf_eps_main, "_delta_", dsu_main)
    colT <- paste0("T_eps_", suf_eps_main, "_delta_", dsu_main)
    
    if (!all(c(colH, colT) %in% names(tbl))) {
      stop("EX1.1 sequences: missing main H/T columns for labeling. Check eps_main/delta config.")
    }
    
    tbl[, direction_label :=
          fifelse(get(colT) >= 0.95, "strong_Tbar",
                  fifelse(get(colT) >= 0.80, "moderate_Tbar",
                          fifelse(get(colT) >= 0.50, "weak_Tbar",
                                  fifelse(get(colH) >= 0.50, "weak_Hbar",
                                          fifelse(get(colH) >= 0.80, "moderate_Hbar",
                                                  fifelse(get(colH) >= 0.95, "strong_Hbar", "neutral"))))))]
    
    setorder(tbl, sequence)
    
    f_out <- file.path(out_dir, paste0("ex1_1_sequences_", tr, ".csv"))
    fwrite(tbl, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}

# Example:
# ex1_1_sequences(cfg)