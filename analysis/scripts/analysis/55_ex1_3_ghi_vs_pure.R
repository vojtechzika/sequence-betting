# ============================================================
# scripts/analysis/55_ex1_3_pure_approximation.R
# (patched: safe correlation draws)
# ============================================================

library(data.table)

ex1_3_ghi_vs_pure <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds     <- as.character(cfg$run$dataset)
  design <- cfg$design
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  stopifnot(!is.null(design$seq), !is.null(design$seq$side_labels), !is.null(design$seq$anchor_labels))
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads))
  stopifnot(length(lab_tails) == 1L, nzchar(lab_tails))
  stopifnot(lab_heads != lab_tails)
  
  pure_heads <- as.character(design$seq$anchor_labels$pure_heads)
  pure_tails <- as.character(design$seq$anchor_labels$pure_tails)
  stopifnot(length(pure_heads) == 1L, nzchar(pure_heads))
  stopifnot(length(pure_tails) == 1L, nzchar(pure_tails))
  stopifnot(pure_heads != pure_tails)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  master <- fread(f_master, encoding = "UTF-8")
  req <- c("pid", "treat", "seq", "stake", "side")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0) stop("master_sequences.csv missing: ", paste(miss, collapse = ", "))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[, side := as.character(side)]
  master[is.na(stake), stake := 0]
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # --- load chi draws ---
    f_chi <- file.path(mod_dir, paste0("ex1_1_participants_", tr, ".rds"))
    if (!file.exists(f_chi)) stop("EX1.3: Missing EX1.1 participants RDS for treatment='", tr, "': ", f_chi)
    
    chi_obj <- readRDS(f_chi)
    
    pid_levels <- NULL
    if (!is.null(chi_obj$pid_levels)) pid_levels <- as.character(chi_obj$pid_levels)
    if (is.null(pid_levels) && !is.null(chi_obj$pid)) pid_levels <- as.character(chi_obj$pid)
    
    chi_draws <- NULL
    if (!is.null(chi_obj$chi_draws)) chi_draws <- chi_obj$chi_draws
    if (is.null(chi_draws) && !is.null(chi_obj$chi_rep)) chi_draws <- chi_obj$chi_rep
    
    if (is.null(pid_levels) || is.null(chi_draws)) {
      stop("EX1.3: ex1_1_participants_<tr>.rds must contain pid_levels (or pid) and chi_draws (or chi_rep).")
    }
    stopifnot(is.matrix(chi_draws), ncol(chi_draws) == length(pid_levels))
    K <- nrow(chi_draws)
    stopifnot(K >= 2L)
    
    # --- compute chi_pure from betting trials on both pure sequences ---
    d <- master[treat == tr & pid %in% pid_levels & seq %in% c(pure_heads, pure_tails)]
    if (nrow(d) == 0) stop("EX1.3: No master rows for pure sequences after filters (tr='", tr, "').")
    
    d <- d[is.finite(stake) & stake > 0]
    d <- d[side %in% c(lab_heads, lab_tails)]
    d[, h := as.integer(side == lab_heads)]
    
    chk <- d[, .N, by = .(pid, seq)]
    bad_chk <- chk[N != 1]
    if (nrow(bad_chk) > 0) {
      stop(
        "EX1.3: expected exactly 1 BET (stake>0) per participant on each pure sequence.\n",
        "Violations (first rows):\n",
        paste(capture.output(print(head(bad_chk, 10))), collapse = "\n")
      )
    }
    
    pure_wide <- dcast(d, pid ~ seq, value.var = "h")
    
    need_cols <- c("pid", pure_heads, pure_tails)
    miss_cols <- setdiff(need_cols, names(pure_wide))
    if (length(miss_cols) > 0) stop("EX1.3: dcast failed; missing columns: ", paste(miss_cols, collapse = ", "))
    
    pure_wide <- pure_wide[is.finite(get(pure_heads)) & is.finite(get(pure_tails))]
    pid_keep <- intersect(pid_levels, pure_wide$pid)
    
    if (length(pid_keep) == 0) stop("EX1.3: No participants have betting trials on BOTH pure sequences (tr='", tr, "').")
    
    if (length(pid_keep) < length(pid_levels)) {
      dropped <- setdiff(pid_levels, pid_keep)
      warning(
        "EX1.3: Dropping ", length(dropped),
        " participant(s) missing a pure-sequence BET (stake>0) in treatment='", tr, "'. Example: ",
        paste(head(dropped, 10), collapse = ", ")
      )
    }
    
    idx_keep <- match(pid_keep, pid_levels)
    stopifnot(!anyNA(idx_keep))
    chi_draws_keep <- chi_draws[, idx_keep, drop = FALSE]
    
    setkey(pure_wide, pid)
    chi_pure <- pure_wide[.(pid_keep),
                          as.integer(get(pure_heads) == 1L) - as.integer(get(pure_tails) == 1L)]
    if (anyNA(chi_pure)) stop("EX1.3: unexpected NA chi_pure after complete-case filtering.")
    
    if (!is.finite(stats::sd(chi_pure)) || stats::sd(chi_pure) <= 0) {
      stop("EX1.3: chi_pure has zero variance after filtering; correlation undefined (tr='", tr, "').")
    }
    
    # ----------------------------
    # SAFE correlation draws
    # ----------------------------
    safe_cor <- function(x, y) {
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3) return(NA_real_)
      sx <- stats::sd(x[ok])
      sy <- stats::sd(y[ok])
      if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) return(NA_real_)
      suppressWarnings(stats::cor(x[ok], y[ok]))
    }
    
    rho_draws <- vapply(seq_len(K), function(k) safe_cor(chi_pure, chi_draws_keep[k, ]), numeric(1))
    
    n_bad <- sum(!is.finite(rho_draws))
    if (n_bad > 0) {
      warning(
        "EX1.3: ", n_bad, " / ", K, " correlation draw(s) were non-finite (typically sd=0 or missing chi draws). ",
        "Summaries will use only finite rho draws."
      )
    }
    
    rho_draws_fin <- rho_draws[is.finite(rho_draws)]
    if (length(rho_draws_fin) < 10) {
      stop("EX1.3: too few finite rho draws (", length(rho_draws_fin), "/", K, ") for treatment='", tr, "'.")
    }
    
    rho_tbl <- data.table(
      dataset = ds,
      treatment = tr,
      N = length(pid_keep),
      K = K,
      K_finite = length(rho_draws_fin),
      rho_median = stats::median(rho_draws_fin),
      rho_mean   = mean(rho_draws_fin),
      rho_q025   = stats::quantile(rho_draws_fin, 0.025),
      rho_q975   = stats::quantile(rho_draws_fin, 0.975),
      p_rho_gt0  = mean(rho_draws_fin > 0)
    )
    
    f_rho <- file.path(out_dir, paste0("ex1_3_rho_", tr, ".csv"))
    fwrite(rho_tbl, f_rho)
    msg("Saved: ", f_rho)
    
    # ----------------------------
    # Grouped chi summaries by chi_pure in {-1,0,1}
    # ----------------------------
    groups <- sort(unique(chi_pure))
    grp_rows <- list()
    
    for (g in groups) {
      idxg <- which(chi_pure == g)
      if (length(idxg) == 0) next
      
      chi_g_draws <- rowMeans(chi_draws_keep[, idxg, drop = FALSE], na.rm = TRUE)
      chi_g_draws <- chi_g_draws[is.finite(chi_g_draws)]
      if (length(chi_g_draws) < 10) next
      
      grp_rows[[as.character(g)]] <- data.table(
        dataset = ds,
        treatment = tr,
        chi_pure = g,
        n = length(idxg),
        chi_median = stats::median(chi_g_draws),
        chi_mean   = mean(chi_g_draws),
        chi_q025   = stats::quantile(chi_g_draws, 0.025),
        chi_q975   = stats::quantile(chi_g_draws, 0.975)
      )
    }
    
    grp_tbl <- rbindlist(grp_rows, use.names = TRUE, fill = TRUE)
    setorder(grp_tbl, chi_pure)
    
    f_grp <- file.path(out_dir, paste0("ex1_3_groups_", tr, ".csv"))
    fwrite(grp_tbl, f_grp)
    msg("Saved: ", f_grp)
    
    outputs[[tr]] <- list(rho = rho_tbl, groups = grp_tbl)
  }
  
  invisible(outputs)
}

# Example:
# ex1_3_ghi_vs_pure(cfg)