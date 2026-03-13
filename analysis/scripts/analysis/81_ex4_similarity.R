# ============================================================
# scripts/analysis/81_ex4_similarity.R
#
# EX4: Similarity between anchor-based and canonical rule-based
# classifications
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

ex4_similarity <- function(cfg) {
  
  ds     <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  
  clean_dir <- path_clean_ds(ds)
  mod_dir   <- path_mod_ds(ds)
  out_dir   <- path_out_ds(ds)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(clean_dir, "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  # ----------------------------------------------------------
  # helpers
  # ----------------------------------------------------------
  
  parse_seq_vec <- function(s) {
    
    ch <- strsplit(as.character(s), "", fixed = TRUE)[[1]]
    
    if (length(ch) != 6L) {
      stop("Sequence must have length 6: ", s)
    }
    
    out <- ifelse(
      ch == lab_heads, 1L,
      ifelse(ch == lab_tails, -1L, NA_integer_)
    )
    
    if (anyNA(out)) {
      stop("Sequence contains unknown symbols: ", s)
    }
    
    out
  }
  
  get_choice_pm1 <- function(dt) {
    
    if ("side_ist" %in% names(dt)) {
      x <- as.integer(dt$side_ist)
      return(ifelse(
        is.na(x), NA_integer_,
        ifelse(x == 1L, 1L,
               ifelse(x == 0L, -1L, NA_integer_))
      ))
    }
    
    if ("side" %in% names(dt)) {
      x <- as.character(dt$side)
      return(ifelse(
        x == lab_heads, 1L,
        ifelse(x == lab_tails, -1L, NA_integer_)
      ))
    }
    
    stop("master_sequences.csv must contain either 'side_ist' or 'side'.")
  }
  
  summ_post <- function(x, prob = TRUE) {
    x <- x[is.finite(x)]
    if (!length(x)) {
      return(c(
        median = NA_real_,
        mean   = NA_real_,
        q025   = NA_real_,
        q975   = NA_real_,
        p_gt0  = NA_real_
      ))
    }
    q <- quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    c(
      median = median(x, na.rm = TRUE),
      mean   = mean(x, na.rm = TRUE),
      q025   = q[1],
      q975   = q[2],
      p_gt0  = if (prob) mean(x > 0, na.rm = TRUE) else NA_real_
    )
  }
  
  summ_scalar <- function(x) {
    c(
      median = x,
      mean   = x,
      q025   = NA_real_,
      q975   = NA_real_,
      p_gt0  = if (is.finite(x)) as.numeric(x > 0) else NA_real_
    )
  }
  
  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (!is.finite(sd(x[ok])) || !is.finite(sd(y[ok])) ||
        sd(x[ok]) <= 0 || sd(y[ok]) <= 0) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok]))
  }
  
  safe_agree <- function(x, y) {
    ok <- is.finite(x) & is.finite(y) & (x != 0) & (y != 0)
    if (!any(ok)) return(NA_real_)
    mean(sign(x[ok]) == sign(y[ok]))
  }
  
  safe_r2_1 <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (!is.finite(sd(x[ok])) || !is.finite(sd(y[ok])) ||
        sd(x[ok]) <= 0 || sd(y[ok]) <= 0) return(NA_real_)
    fit <- lm(y[ok] ~ x[ok])
    unname(summary(fit)$r.squared)
  }
  
  safe_r2_multi <- function(X, y) {
    ok <- is.finite(y) & rowSums(!is.finite(X)) == 0L
    if (sum(ok) < 3L) return(NA_real_)
    if (!is.finite(sd(y[ok])) || sd(y[ok]) <= 0) return(NA_real_)
    
    Xok <- X[ok, , drop = FALSE]
    if (ncol(Xok) == 0L) return(NA_real_)
    if (any(!is.finite(apply(Xok, 2, sd))) || any(apply(Xok, 2, sd) <= 0)) {
      return(NA_real_)
    }
    
    fit <- tryCatch(
      lm(y[ok] ~ Xok),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    
    unname(summary(fit)$r.squared)
  }
  
  # ----------------------------------------------------------
  # load master once
  # ----------------------------------------------------------
  
  master <- fread(f_master, encoding = "UTF-8")
  
  req <- c("pid", "treat", "seq", "stake")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0) {
    stop("master_sequences.csv missing columns: ", paste(miss, collapse = ", "))
  }
  
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq   := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  
  master[, choice_pm1 := get_choice_pm1(.SD)]
  
  if ("screen_ms" %in% names(master)) {
    master[, screen_ms := as.numeric(screen_ms)]
  }
  if ("lotr_score" %in% names(master)) {
    master[, lotr_score := as.numeric(lotr_score)]
  }
  
  dt0 <- master[stake > 0 & is.finite(choice_pm1)]
  
  outputs <- list()
  
  rules <- c("q_last", "q_run2", "q_run3", "q_imb", "q_alt")
  rule_labels <- sub("^q_", "", rules)
  
  for (tr in tr_vec) {
    
    msg("Running EX4 similarity for treatment: ", tr)
    
    f_ex1_seq <- file.path(mod_dir, paste0("ex1_1_", tr, "_sequences.rds"))
    f_ex1_pid <- file.path(mod_dir, paste0("ex1_1_", tr, "_participants.rds"))
    
    f_seq_sum  <- file.path(out_dir, paste0("ex4_sequence_summary_", tr, ".csv"))
    f_seq_diag <- file.path(out_dir, paste0("ex4_sequence_diagnostics_", tr, ".csv"))
    f_pid_sum  <- file.path(out_dir, paste0("ex4_participant_summary_", tr, ".csv"))
    f_pred     <- file.path(out_dir, paste0("ex4_rule_predictors_", tr, ".csv"))
    f_r2       <- file.path(out_dir, paste0("ex4_sequence_r2_", tr, ".csv"))
    
    skip_all <- should_skip(
      paths = c(f_seq_sum, f_seq_diag, f_pid_sum, f_pred, f_r2),
      cfg   = cfg,
      type  = "output",
      label = paste0("EX4 similarity (", ds, "/", tr, ")")
    )
    if (skip_all) next
    
    stopifnot(file.exists(f_ex1_seq), file.exists(f_ex1_pid))
    
    ex1_seq <- readRDS(f_ex1_seq)
    ex1_pid <- readRDS(f_ex1_pid)
    
    chi_s_draws <- ex1_seq$chi_draws
    chi_i_draws <- ex1_pid$chi_draws
    
    seq_levels <- as.character(ex1_seq$seq_levels)
    pid_levels <- as.character(ex1_pid$pid_levels)
    
    stopifnot(is.matrix(chi_s_draws), is.matrix(chi_i_draws))
    stopifnot(ncol(chi_s_draws) == length(seq_levels))
    stopifnot(ncol(chi_i_draws) == length(pid_levels))
    
    dt <- copy(dt0[treat == tr])
    if (nrow(dt) == 0L) next
    
    # --------------------------------------------------------
    # sequence features
    # --------------------------------------------------------
    
    seq_map <- data.table(seq = seq_levels)
    seq_mat <- t(vapply(seq_levels, parse_seq_vec, integer(6L)))
    
    seq_map[, `:=`(
      x1 = seq_mat[, 1],
      x2 = seq_mat[, 2],
      x3 = seq_mat[, 3],
      x4 = seq_mat[, 4],
      x5 = seq_mat[, 5],
      x6 = seq_mat[, 6]
    )]
    
    seq_map[, B := x1 + x2 + x3 + x4 + x5 + x6]
    seq_map[, A := (x2 != x1) + (x3 != x2) + (x4 != x3) + (x5 != x4) + (x6 != x5)]
    seq_map[, run_len := {
      rr <- rle(c(x1, x2, x3, x4, x5, x6))
      rr$lengths[length(rr$lengths)]
    }, by = seq]
    
    dt <- merge(dt, seq_map, by = "seq", all.x = TRUE)
    
    # --------------------------------------------------------
    # canonical rules
    # --------------------------------------------------------
    
    dt[, q_last := choice_pm1 * x6]
    dt[, q_run2 := ifelse(run_len >= 2L, choice_pm1 * x6, NA_real_)]
    dt[, q_run3 := ifelse(run_len >= 3L, choice_pm1 * x6, NA_real_)]
    dt[, q_imb  := ifelse(B == 0L, NA_real_, -choice_pm1 * sign(B))]
    dt[, q_alt  := ifelse(
      A <= 2L, choice_pm1 * x6,
      ifelse(A >= 4L, -choice_pm1 * x6, NA_real_)
    )]
    
    # --------------------------------------------------------
    # sequence rule scores Q_s^(k)
    # --------------------------------------------------------
    
    q_seq_list <- lapply(seq_along(rules), function(i) {
      r <- rules[i]
      out <- dt[!is.na(get(r)),
                .(Q = mean(get(r)), N = .N),
                by = seq]
      out[, rule := rule_labels[i]]
      out
    })
    
    q_seq <- rbindlist(q_seq_list, fill = TRUE)
    
    q_seq <- merge(
      CJ(seq = seq_levels, rule = rule_labels, unique = TRUE),
      q_seq,
      by = c("seq", "rule"),
      all.x = TRUE,
      sort = FALSE
    )
    
    q_seq <- merge(
      data.table(seq = seq_levels, seq_id = seq_along(seq_levels)),
      q_seq,
      by = "seq",
      all.x = TRUE,
      sort = FALSE
    )
    
    # --------------------------------------------------------
    # participant rule scores Q_i^(k)
    # --------------------------------------------------------
    
    q_pid_list <- lapply(seq_along(rules), function(i) {
      r <- rules[i]
      out <- dt[!is.na(get(r)),
                .(Q = mean(get(r)), N = .N),
                by = pid]
      out[, rule := rule_labels[i]]
      out
    })
    
    q_pid <- rbindlist(q_pid_list, fill = TRUE)
    
    q_pid <- merge(
      CJ(pid = pid_levels, rule = rule_labels, unique = TRUE),
      q_pid,
      by = c("pid", "rule"),
      all.x = TRUE,
      sort = FALSE
    )
    
    # --------------------------------------------------------
    # sequence similarity
    # --------------------------------------------------------
    
    seq_sum_rows  <- list()
    seq_diag_rows <- list()
    
    for (rule_nm in rule_labels) {
      
      qs <- q_seq[rule == rule_nm][match(seq_levels, seq), Q]
      
      rho_draws <- apply(
        chi_s_draws, 1L,
        function(chi) safe_cor(qs, chi)
      )
      
      agr_draws <- apply(
        chi_s_draws, 1L,
        function(chi) safe_agree(qs, chi)
      )
      
      seq_sum_rows[[paste0(rule_nm, "_rho")]] <- data.table(
        level  = "sequence",
        rule   = rule_nm,
        metric = "correlation",
        t(summ_post(rho_draws, TRUE))
      )
      
      seq_sum_rows[[paste0(rule_nm, "_agree")]] <- data.table(
        level  = "sequence",
        rule   = rule_nm,
        metric = "agreement",
        t(summ_post(agr_draws, FALSE))
      )
      
      p_match <- vapply(
        seq_along(seq_levels),
        function(j) {
          if (!is.finite(qs[j]) || qs[j] == 0) return(NA_real_)
          mean(sign(qs[j]) == sign(chi_s_draws[, j]) &
                 chi_s_draws[, j] != 0,
               na.rm = TRUE)
        },
        numeric(1L)
      )
      
      seq_diag_rows[[rule_nm]] <- data.table(
        rule     = rule_nm,
        sequence = seq_levels,
        seq_id   = seq_along(seq_levels),
        Q_seq    = qs,
        p_match  = p_match
      )
    }
    
    tbl_seq_sum  <- rbindlist(seq_sum_rows, fill = TRUE)
    tbl_seq_diag <- rbindlist(seq_diag_rows, fill = TRUE)
    
    tbl_seq_sum[, `:=`(dataset = ds, treatment = tr)]
    setcolorder(tbl_seq_sum, c(
      "dataset", "treatment", "level", "rule", "metric",
      "median", "mean", "q025", "q975", "p_gt0"
    ))
    
    tbl_seq_diag[, `:=`(dataset = ds, treatment = tr)]
    setcolorder(tbl_seq_diag, c(
      "dataset", "treatment", "rule", "sequence", "seq_id", "Q_seq", "p_match"
    ))
    
    # --------------------------------------------------------
    # participant similarity
    # --------------------------------------------------------
    
    pid_sum_rows <- list()
    
    for (rule_nm in rule_labels) {
      
      qi <- q_pid[rule == rule_nm][match(pid_levels, pid), Q]
      
      rho_draws <- apply(
        chi_i_draws, 1L,
        function(chi) safe_cor(qi, chi)
      )
      
      agr_draws <- apply(
        chi_i_draws, 1L,
        function(chi) safe_agree(qi, chi)
      )
      
      pid_sum_rows[[paste0(rule_nm, "_rho")]] <- data.table(
        level  = "participant",
        rule   = rule_nm,
        metric = "correlation",
        t(summ_post(rho_draws, TRUE))
      )
      
      pid_sum_rows[[paste0(rule_nm, "_agree")]] <- data.table(
        level  = "participant",
        rule   = rule_nm,
        metric = "agreement",
        t(summ_post(agr_draws, FALSE))
      )
    }
    
    tbl_pid_sum <- rbindlist(pid_sum_rows, fill = TRUE)
    
    tbl_pid_sum[, `:=`(dataset = ds, treatment = tr)]
    setcolorder(tbl_pid_sum, c(
      "dataset", "treatment", "level", "rule", "metric",
      "median", "mean", "q025", "q975", "p_gt0"
    ))
    
    # --------------------------------------------------------
    # exploratory predictors of canonical rules
    # non-posterior, descriptive only
    # --------------------------------------------------------
    
    dtr <- master[treat == tr & pid %in% pid_levels]
    
    d_rt <- dtr[stake > 0 & is.finite(screen_ms) & screen_ms > 0]
    rt_pid <- d_rt[, .(rt_log_mean = mean(log(screen_ms))), by = pid]
    opt_pid <- dtr[, .(lotr_score = lotr_score[1]), by = pid]
    
    pred_pid <- merge(
      data.table(pid = pid_levels),
      opt_pid,
      by = "pid",
      all.x = TRUE
    )
    pred_pid <- merge(pred_pid, rt_pid, by = "pid", all.x = TRUE)
    
    chi_i_med <- apply(chi_i_draws, 2, median, na.rm = TRUE)
    
    q_pid_wide <- dcast(
      q_pid[, .(pid, rule, Q)],
      pid ~ rule,
      value.var = "Q"
    )
    
    pred_tbl <- merge(pred_pid, q_pid_wide, by = "pid", all.x = TRUE)
    pred_tbl[, chi_i_median := chi_i_med[match(pid, pid_levels)]]
    
    q_rule_names <- sort(setdiff(names(q_pid_wide), "pid"))
    
    pred_rows <- list()
    
    for (rule_nm in q_rule_names) {
      
      qv <- pred_tbl[[rule_nm]]
      
      pred_rows[[paste0(rule_nm, "_chi")]] <- data.table(
        rule = rule_nm,
        predictor = "chi_i_median",
        t(summ_scalar(safe_cor(qv, pred_tbl$chi_i_median)))
      )
      
      pred_rows[[paste0(rule_nm, "_opt")]] <- data.table(
        rule = rule_nm,
        predictor = "lotr_score",
        t(summ_scalar(safe_cor(qv, pred_tbl$lotr_score)))
      )
      
      pred_rows[[paste0(rule_nm, "_rt")]] <- data.table(
        rule = rule_nm,
        predictor = "rt_log_mean",
        t(summ_scalar(safe_cor(qv, pred_tbl$rt_log_mean)))
      )
    }
    
    tbl_pred <- rbindlist(pred_rows, fill = TRUE)
    tbl_pred[, `:=`(dataset = ds, treatment = tr)]
    setcolorder(tbl_pred, c(
      "dataset", "treatment", "rule", "predictor",
      "median", "mean", "q025", "q975", "p_gt0"
    ))
    
    # --------------------------------------------------------
    # sequence-level variance explained
    # --------------------------------------------------------
    
    seq_r2_rows <- list()
    
    for (rule_nm in rule_labels) {
      
      qs <- q_seq[rule == rule_nm][match(seq_levels, seq), Q]
      
      r2_draws <- apply(
        chi_s_draws, 1L,
        function(chi) safe_r2_1(qs, chi)
      )
      
      seq_r2_rows[[paste0(rule_nm, "_r2")]] <- data.table(
        level  = "sequence",
        model  = "single_rule",
        rule   = rule_nm,
        metric = "r_squared",
        t(summ_post(r2_draws, FALSE))
      )
    }
    
    q_seq_wide <- dcast(
      q_seq[, .(seq, rule, Q)],
      seq ~ rule,
      value.var = "Q"
    )
    
    X_joint <- as.matrix(
      q_seq_wide[match(seq_levels, seq), setdiff(names(q_seq_wide), "seq"), with = FALSE]
    )
    
    r2_joint_draws <- apply(
      chi_s_draws, 1L,
      function(chi) safe_r2_multi(X_joint, chi)
    )
    
    seq_r2_rows[["joint"]] <- data.table(
      level  = "sequence",
      model  = "joint_rules",
      rule   = "all",
      metric = "r_squared",
      t(summ_post(r2_joint_draws, FALSE))
    )
    
    tbl_r2 <- rbindlist(seq_r2_rows, fill = TRUE)
    tbl_r2[, `:=`(dataset = ds, treatment = tr)]
    setcolorder(tbl_r2, c(
      "dataset", "treatment", "level", "model", "rule", "metric",
      "median", "mean", "q025", "q975", "p_gt0"
    ))
    
    # --------------------------------------------------------
    # save
    # --------------------------------------------------------
    
    fwrite(tbl_seq_sum,  f_seq_sum)
    fwrite(tbl_seq_diag, f_seq_diag)
    fwrite(tbl_pid_sum,  f_pid_sum)
    fwrite(tbl_pred,     f_pred)
    fwrite(tbl_r2,       f_r2)
    
    msg("Saved: ", f_seq_sum)
    msg("Saved: ", f_seq_diag)
    msg("Saved: ", f_pid_sum)
    msg("Saved: ", f_pred)
    msg("Saved: ", f_r2)
    
    outputs[[tr]] <- list(
      sequence_summary     = tbl_seq_sum,
      sequence_diagnostics = tbl_seq_diag,
      participant_summary  = tbl_pid_sum,
      rule_predictors      = tbl_pred,
      sequence_r2          = tbl_r2
    )
  }
  
  invisible(outputs)
}